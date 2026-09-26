#!/usr/bin/env Rscript
## Purpose: Build undirected meta-node scaffolds for Gabi (1.0.1) from the
## Reactome Functional Interaction network, in batch. Runs over every cell line,
## writing one GraphML scaffold per cell line at the d0 (no-threshold) cutoff for
## Gabi to orient.
##
## Reactome/HGNC counterpart of 1-undir_scaff_humannet.R, differing only in two
## inputs: gene-level edges come from Reactome FI (HGNC symbols) instead of
## HumanNet-FN (Entrez), and annotation is joined on the HGNC column. The shared
## contraction + density logic lives in ../scaffold_common.R.
##
## Pipeline (per cell line): gene-level Reactome FI edges (read once) -> keep genes
## present in the annotation (matched by HGNC) -> label each with its
## Cluster_Profile -> contract genes into modules, counting gene-gene edges
## between each module pair -> keep a module edge A-B only if
## n_edges(A,B) / (|A|*|B|) >= threshold -> write one plain undirected GraphML.
##
## Usage: Rscript 1-undir_scaff_reactome.R   (run from gabi-maboss/scripts/reactome/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/gabi-maboss/scripts/reactome")
getwd()   # sanity check (kept for posterity)

## Graph + GraphML I/O; banners hidden.
suppressMessages(library(igraph))

## Shared build_metagraph() + logmsg().
source("../scaffold_common.R")

# 1. INPUTS --------------------------------------------------------------------

## 1a. Reactome FI gene-gene edge list: tab-separated WITH a header (Gene1, Gene2,
##     Annotation, Direction, Score); columns 1-2 are the two HGNC symbols of an
##     undirected edge, other columns (incl. direction) ignored. Read as character
##     so symbols match exactly. Shared, so read once (section 2).
reactome_file       <- "../../data/FIsInGene_04142025_with_annotations.txt"
reactome_has_header <- TRUE

## 1b. Cell lines to process (each needs a matching annotation file below).
cell_lines <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")

## 1c. Per-cell-line module table (cyto_module HGNC output); %s = cell-line code.
##     Columns: HGNC (join key to Reactome FI), Cluster_Number, Multi_Profiles,
##     Cluster_Profile (module label).
annot_tmpl   <- "../../data/cyto_module_output/hgnc/%s_genes_in_modules.csv"
col_gene     <- "HGNC"             # HGNC symbol column (join key to Reactome FI)
col_metanode <- "Cluster_Profile"  # Module label column, read verbatim

## 1d. Edge-density threshold: keep a scaffold edge A-B only if the fraction of
##     the |A|*|B| possible gene-gene pairs that are actually linked reaches thr.
##     Reactome starts with thr = 0 (no threshold).
thr     <- 0
thr_str <- formatC(thr, format = "g")   # e.g. "0"; avoids scientific notation

## 1e. Output. One flat dir; the cell-line code + threshold are carried in the
##     filename. Kept separate from the HumanNet scaffolds so neither clobbers
##     the other. %s slots in out_tmpl are the cell-line code then the threshold.
out_dir  <- "../../output/scaffold/reactome"
out_tmpl <- "%s_reactome_metanode_scaffold_d%s.graphml"

## 1f. Logs. One plain-text log per cell line (gene/module/edge counts, density
##     quantiles, per-threshold edge/isolated counts). %s = cell-line code.
log_dir  <- "../../logs/scaffold/reactome"
log_tmpl <- "%s_scaffold.log"

# 2. SHARED SETUP --------------------------------------------------------------

## Gene-level undirected graph from Reactome FI, read once and reused for every
## cell line (only this read is cell-line-invariant; annotation varies per line).
## Read the edge list as strings so gene symbols match the annotation exactly;
## quote/comment disabled so stray "/# chars can't break parsing.
fi <- read.table(reactome_file, sep = "\t", header = reactome_has_header,
                 colClasses = "character", quote = "", comment.char = "")

## Undirected graph from the two symbol columns only; simplify() drops duplicate
## edges and self-loops for a clean gene network.
g_gene_full <- simplify(graph_from_data_frame(fi[, 1:2], directed = FALSE))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# 3. RUN: every cell line ------------------------------------------------------

for (cl in cell_lines) {
  annot_file <- sprintf(annot_tmpl, cl)
  log_con    <- file(file.path(log_dir, sprintf(log_tmpl, cl)), open = "wt")

  ## Skip a missing cell line without aborting the batch.
  if (!file.exists(annot_file)) {
    logmsg(log_con, "[%s] SKIP - annotation not found: %s\n", cl, annot_file)
    close(log_con)
    next
  }

  built  <- build_metagraph(g_gene_full, annot_file, col_gene, col_metanode)
  g_meta <- built$g_meta

  logmsg(log_con, "[%s] genes used: %d | meta-nodes: %d | candidate edges: %d\n",
         cl, built$n_genes, vcount(g_meta), ecount(g_meta))

  ## Density spread helps choose thresholds: a big gap between the 75th/90th
  ## percentiles flags a long tail of near-saturated pairs above a mass of weak
  ## ones, so a fixed thr cuts roughly half the candidate edges while isolating
  ## few/no extra modules.
  if (ecount(g_meta))
    logmsg(log_con, "[%s] density quantiles (50/75/90/max): %.4f / %.4f / %.4f / %.4f\n",
           cl,
           quantile(E(g_meta)$density, 0.50),
           quantile(E(g_meta)$density, 0.75),
           quantile(E(g_meta)$density, 0.90),
           max(E(g_meta)$density))

  ## Keep edges meeting the threshold; keep all modules (even if now edgeless).
  g_out <- subgraph_from_edges(g_meta, E(g_meta)[E(g_meta)$density >= thr],
                               delete.vertices = FALSE)
  ## Output stays a plain scaffold (no weights, no density) so Gabi can read it.
  g_out <- delete_edge_attr(g_out, "weight")
  g_out <- delete_edge_attr(g_out, "density")

  out_path <- file.path(out_dir, sprintf(out_tmpl, cl, thr_str))
  write_graph(g_out, out_path, format = "graphml")

  ## Modules left edgeless after the filter: still present, ignored by Gabi.
  iso <- sum(degree(g_out) == 0L)
  logmsg(log_con, "[%s]   thr=%-5s edges=%-6d isolated=%-4d -> %s\n",
         cl, thr_str, ecount(g_out), iso, out_path)

  close(log_con)
}

cat("Done.\n")

## Downstream note: 2-gabi_reactome.R wants expression as SAMPLES-in-rows,
## NODES-in-columns; the '*_cluster_profile_mean_exp.csv' file is the transpose,
## and that reshape now happens inside 2-gabi_reactome.R. Its row labels are the
## same Cluster_Profile values used here, so node names line up by construction.

# END SCRIPT -------------------------------------------------------------------
