#!/usr/bin/env Rscript
## Purpose: Derive degree/role/hub structure from the per-cell-line signed module
## networks produced by 3-sign_export_humannet.R / 3-sign_export_reactome.R /
## 3-sign_export_huri.R, after they've been loaded into Cytoscape, analysed, and
## exported back out (4-export_cyto_humannet.R / 4-export_cyto_reactome.R /
## 4-export_cyto_huri.R). Branch-agnostic: it walks all three branches in one run
## and tags each network's output by branch, so it isn't duplicated per branch
## the way the earlier steps are.
##
## For each network: compute in/out degree and the out/in ratio, classify each
## node as Source/Sink/Relay, and flag potential hubs. A hub needs all three of:
## out-degree fan-out AND betweenness (bottleneck) in this network's own top 10%,
## AND net out-flow (out/in ratio > 1) - the cutoffs adapt per network rather
## than being fixed constants. Role logic follows stem-gabi's sort_modules,
## adapted to the signed-network columns (rho/sign instead of weight).
##
## Input: any "<CL>_..._cyto.graphml" under output/gabi/<branch>/<CL>/ (branch =
## humannet/reactome/huri), i.e. whatever you've exported from Cytoscape so far,
## from any or all branches. Branch, cell line, 'directed' vs 'backbone', and the
## threshold are all read out of the file's path/name - keep those words in your
## Cytoscape export file names so this script can tell networks apart.
##
## Usage: Rscript 5-sort_modules.R   (run from gabi-maboss/scripts/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/gabi-maboss/scripts")
getwd()   # sanity check (kept for posterity)

## Graph I/O; banners hidden.
suppressMessages(library(igraph))

# 1. INPUTS --------------------------------------------------------------------

## Gabi output root, covering all three branches (humannet/reactome/huri); every
## matching file found anywhere under here is picked up (however many cell
## lines/branches/thresholds/types exported so far).
gabi_dir <- "../output/gabi/"

## Where per-network node/edge tables (and the cross-network summary) land.
## Per-network tables go in a branch subfolder (humannet/reactome/huri) so the
## folder doesn't fill up with ~50 flat CSVs; the cross-network summary stays
## at the top level since it spans all three branches.
output_dir <- "../output/module_sorting/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## Every Cytoscape re-export, regardless of cell line/threshold/network type.
graph_files <- list.files(gabi_dir, pattern = "_cyto\\.graphml$",
                          recursive = TRUE, full.names = TRUE)

if (length(graph_files) == 0)
  stop("No *_cyto.graphml files found under ", gabi_dir,
       " - export your analysed network(s) from Cytoscape first.")

# 2. FUNCTION: IDENTIFY ONE FILE -----------------------------------------------

## Cell line = the file's parent folder name; branch = the grandparent folder
## name (humannet/reactome/huri, per gabi_dir's layout above); 'directed' (full
## signed net) vs 'backbone' and the threshold suffix (d0, d0.05, d0.1, d0.2, ...)
## are read out of the file name itself, since both can be pushed to Cytoscape
## per cell line. Branch is folded into 'code' so identically-thresholded
## networks from different branches (e.g. Reactome and HuRI both at d0) don't
## collide.
describe_network_file <- function(file_path) {
  file_name <- basename(file_path)
  cl        <- basename(dirname(file_path))
  branch    <- basename(dirname(dirname(file_path)))

  net_type <- if (grepl("backbone", file_name, ignore.case = TRUE)) "backbone" else "directed"

  thr_hit <- regmatches(file_name, regexpr("(?<=_d)[0-9]+(\\.[0-9]+)?", file_name, perl = TRUE))
  thr     <- if (length(thr_hit) && nzchar(thr_hit)) thr_hit else NA_character_

  code <- paste(c(cl, branch, net_type, if (!is.na(thr)) paste0("d", thr)), collapse = "_")

  list(file_path = file_path, file_name = file_name, cl = cl, branch = branch,
       net_type = net_type, thr = thr, code = code)
}

# 3. FUNCTION: PROCESS ONE NETWORK ---------------------------------------------

## Per-node tally of a node vector against the full node set, e.g. how many
## activation edges leave each node. Returns one integer per node in 'levels' order.
count_by_node <- function(nodes, levels) as.integer(table(factor(nodes, levels = levels)))

## Reads one Cytoscape-exported graphml, derives degree/role/hub info per node,
## writes a node table, and returns role/hub counts for the cross-network summary.
process_network <- function(info, output_dir) {
  cat(sprintf("[%s] Reading %s\n", info$code, info$file_name))
  g <- read_graph(info$file_path, format = "graphml")

  ## Edges carry whatever the signing pipeline + Cytoscape attached (score,
  ## bidirected, rho, sign, sign_num, plus Cytoscape extras). 'sign' states
  ## activation/inhibition directly, so no derived +/- column is needed. Pulled
  ## up front so the per-node sign tallies below can reuse it.
  edge_df <- igraph::as_data_frame(g, what = "edges")

  ## Out/in degree and their ratio. In-degree 0 makes the ratio undefined: a
  ## pure source (out > 0) is Inf; a node with no edges at all is 0 (so neither
  ## trips the ratio > 1 hub test below).
  out_deg <- degree(g, mode = "out")
  in_deg  <- degree(g, mode = "in")
  out_in_ratio <- ifelse(in_deg == 0 & out_deg > 0, Inf,
                  ifelse(in_deg == 0 & out_deg == 0, 0,
                         out_deg / in_deg))

  ## Roles: Source (only out), Sink (only in), Relay (everything else, incl. isolated).
  role <- ifelse(in_deg == 0 & out_deg > 0, "Source",
          ifelse(in_deg > 0 & out_deg == 0, "Sink",
                 "Relay"))

  ## Hub = top-10% out-degree (fan-out) AND top-10% betweenness (bottleneck) AND
  ## net out-flow. Cutoffs are this network's own 90th percentile, so the rule
  ## adapts per cell line/threshold export.
  btw         <- betweenness(g, directed = TRUE)
  out_deg_cut <- quantile(out_deg, 0.90, na.rm = TRUE)
  btw_cut     <- quantile(btw, 0.90, na.rm = TRUE)
  hub         <- ifelse(out_deg >= out_deg_cut & btw >= btw_cut & out_in_ratio > 1, "Y", "N")

  ## Biological cluster number (first part of the module name, e.g. "6.3_12" -> "6").
  node_names      <- names(in_deg)
  cluster_numbers <- sub("\\..*", "", node_names)

  ## Per-node counts of outgoing/incoming edges by sign.
  act <- edge_df$sign == "activation"
  inh <- edge_df$sign == "inhibition"

  node_info <- data.frame(
    name           = node_names,
    Cluster_Number = cluster_numbers,
    indegree       = in_deg,
    outdegree      = out_deg,
    out_in_ratio   = out_in_ratio,
    betweenness    = btw,
    EdgeCount      = in_deg + out_deg,   # total degree, for the gabi_style.xml NODE_FILL_COLOR mapping
    Out_Activation = count_by_node(edge_df$from[act], node_names),
    Out_Inhibition = count_by_node(edge_df$from[inh], node_names),
    In_Activation  = count_by_node(edge_df$to[act], node_names),
    In_Inhibition  = count_by_node(edge_df$to[inh], node_names),
    role           = role,
    Potential_Hub  = hub,
    stringsAsFactors = FALSE
  )

  branch_dir <- file.path(output_dir, info$branch)
  if (!dir.exists(branch_dir)) dir.create(branch_dir, recursive = TRUE)
  write.csv(node_info, file = file.path(branch_dir, paste0(info$code, "_node_info.csv")), row.names = FALSE)

  ## For the cross-network summary.
  role_counts <- table(factor(role, levels = c("Relay", "Sink", "Source")))
  hub_count   <- sum(node_info$Potential_Hub == "Y")

  list(role_counts = as.numeric(role_counts), hub_count = hub_count)
}

# 4. RUN: every Cytoscape-exported network -------------------------------------

## Lists keyed by 'code' (e.g. "AH_humannet_directed_d0.2", "AH_reactome_backbone_d0").
role_summary <- list()
hub_summary  <- list()

for (file_path in graph_files) {
  info <- describe_network_file(file_path)

  ## Two exports resolving to the same code would silently overwrite each other's
  ## output - suffix and warn instead so it's obvious the export names need fixing.
  if (info$code %in% names(role_summary)) {
    dup_n <- sum(grepl(paste0("^", info$code, "(_dup[0-9]+)?$"), names(role_summary))) + 1
    warning(sprintf(
      "[%s] Code collides with an earlier file - rename Cytoscape exports to include 'directed'/'backbone' (and the threshold) so they don't overwrite each other. Suffixing this one as _dup%d.",
      info$file_path, dup_n))   # branch is already folded into 'code'; a collision means two exports share cell line/branch/type/threshold
    info$code <- paste0(info$code, "_dup", dup_n)
  }

  result <- process_network(info, output_dir)
  role_summary[[info$code]] <- result$role_counts
  hub_summary[[info$code]]  <- result$hub_count
}

## Hub counts as a one-row data frame.
hub_summary_df <- as.data.frame(do.call(cbind, hub_summary))
hub_summary_df <- cbind(role = "Hubs", hub_summary_df)

## Role counts as a three-row data frame.
role_summary_df <- as.data.frame(do.call(cbind, role_summary))
role_summary_df <- cbind(role = c("Relay", "Sink", "Source"), role_summary_df)

## Combine and write the cross-network summary.
combined_summary <- rbind(role_summary_df, hub_summary_df)
write.csv(combined_summary, file = paste0(output_dir, "all_net_role_count.csv"), row.names = FALSE)

cat("Done.\n")

# END SCRIPT -------------------------------------------------------------------
