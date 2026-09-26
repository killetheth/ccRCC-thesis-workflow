#!/usr/bin/env Rscript
## Purpose: Annotate every cell line's directed Gabi network in one pass, then
## export each to GraphML (and optionally Cytoscape).
##
## Pipeline (per cell line): load the directed net (result$directed.net) +
## expression -> join 'score' and 'bidirected' from result$scores onto each edge,
## matched on the undirected node pair -> add an activation/inhibition 'sign' +
## the underlying 'rho' from endpoint co-expression (Spearman; rho >= 0 ->
## activation) -> write GraphML for the full signed net AND the resolved backbone
## (bidirected == FALSE) -> optionally push the styled signed net to Cytoscape.
##
## Result: every edge carries direction (Gabi) + confidence (score) + sign (rho)
## - everything the Boolean rules need.
##
## Why rho is computed HERE rather than read off the Gabi result: not because
## ClNetScOn cannot sign - it can, and normally would. gabi.directed() ends by
## giving each scaffold-only edge a Spearman rho as its `weight`. That block never
## fires in this pipeline, because 2-gabi_<branch>.R hands the scaffold over as a
## 2-column TSV, leaving the graph with no `weight` attribute for its sentinel
## test to find; it then skips silently, warning included. See that script's
## header for the full account. What is computed below is the same quantity over
## the same 18 samples, so these are the signs ClNetScOn would have assigned.
##
## Caveat for methods: 'sign' is contemporaneous co-expression across the 18
## timepoint x replicate samples, not a lagged/mechanistic sign, and it is
## thresholded at zero, so rho = 0.001 is recorded as activation on the same
## footing as rho = 0.90. Both properties belong to the statistic, not to this
## implementation - ClNetScOn's own sign is the same Spearman coefficient.
## Constant-expression modules give rho = NA (reported) and need a
## manual/default call.
##
## Threshold is a variable (default 0 for reactome); change 'thr' in section 1
## to sign a different scaffold density. Cell lines not yet oriented at 'thr' are
## skipped.
##
## Reactome counterpart of 3-sign_export_humannet.R (thr 0 vs 0.2, reactome/ vs
## humannet/ output subtree).
##
## Usage: Rscript 3-sign_export_reactome.R   (run from gabi-maboss/scripts/reactome/;
## if push_to_cytoscape is TRUE, Cytoscape must already be open with the target
## session loaded)

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

# 1. INPUTS --------------------------------------------------------------------

## Cell lines to process; each needs its own ../../output/gabi/reactome/<CL>/ folder.
cell_lines <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")

## Scaffold density threshold to sign+export; matches 1-undir_scaff_reactome.R's
## d<thr> suffix and 2-gabi_reactome.R's orientation output.
thr     <- 0
thr_str <- formatC(thr, format = "g")   # string form used in every file name below

## Per-cell-line Gabi output folder; %s = cell-line code.
gabi_dir <- "../../output/gabi/reactome/%s"

## Gabi result object (list with $directed.net and $scores); written by 2-gabi_reactome.R.
result_tmpl <- "%s_reactome_gabi_result_d%s.rds"

## Expression matrix, samples (rows) x meta-nodes (columns), used for 'sign'.
expr_tmpl <- "%s_reactome_expr_gabi.tsv"

## Output names: full signed network (every edge) and resolved backbone (bidirected == FALSE).
graphml_full_tmpl     <- "%s_reactome_directed_signed_d%s.graphml"
graphml_backbone_tmpl <- "%s_reactome_backbone_signed_d%s.graphml"

## Whether to also push each cell line's signed networks to the Cytoscape session
## that's already open. Session is left unsaved - save it yourself in Cytoscape
## once you've checked the result.
push_to_cytoscape <- TRUE

## Collection (within that session) to push each cell line's pair of networks
## into - one per cell line, named to match the existing convention in the
## session. Must match the Network-panel name exactly (case/spaces/parentheses)
## or a re-run creates a sibling collection instead of adding to this one.
cytoscape_collection <- c(
  AH = "A498 - Hypoxia (AH)",
  AN = "A498 - Normoxia (AN)",
  CH = "CAKI1 - Hypoxia (CH)",
  CN = "CAKI1 - Normoxia (CN)",
  HH = "HUVECS - Hypoxia (HH)",
  HN = "HUVECS - Normoxia (HN)",
  UH = "UMRC2 - Hypoxia (UH)",
  UN = "UMRC2 - Normoxia (UN)"
)

# 2. FUNCTION: SIGN + EXPORT ONE CELL LINE -------------------------------------

## Full annotate-then-export pipeline for a single cell line 'cl' at 'thr_str';
## returns invisibly after skipping if that cell line hasn't been oriented at
## this threshold yet, or if the expression matrix is missing a needed column.
sign_export_cell_line <- function(cl, thr_str, gabi_dir, result_tmpl, expr_tmpl,
                                  graphml_full_tmpl, graphml_backbone_tmpl,
                                  push_to_cytoscape, cytoscape_ok, cytoscape_collection) {

  cl_dir     <- sprintf(gabi_dir, cl)
  result_rds <- file.path(cl_dir, sprintf(result_tmpl, cl, thr_str))

  ## Not every cell line has reached this threshold yet - skip, don't abort.
  if (!file.exists(result_rds)) {
    cat(sprintf("[%s] SKIP - no d%s result yet: %s\n", cl, thr_str, result_rds))
    return(invisible(NULL))
  }

  expr_tsv <- file.path(cl_dir, sprintf(expr_tmpl, cl))
  res      <- readRDS(result_rds)   # carries $directed.net and $scores
  g        <- res$directed.net      # the directed meta-node network to annotate

  ## Expression matrix; rows = sample, columns = meta-node, names kept verbatim.
  expr <- as.matrix(read.table(expr_tsv, sep = "\t", header = TRUE,
                               row.names = 1, check.names = FALSE))

  cat(sprintf("[%s] Directed net: %d nodes, %d edges\n", cl, vcount(g), ecount(g)))

  el   <- as_edgelist(g)                                # from/to meta-node names
  miss <- setdiff(unique(as.vector(el)), colnames(expr))  # endpoints without expression

  ## Skip this cell line rather than aborting the rest of the batch.
  if (length(miss)) {
    cat(sprintf("[%s] SKIP - no expression column for: %s\n",
                cl, paste(head(miss), collapse = ", ")))
    return(invisible(NULL))
  }

  ## Per-edge key, endpoints sorted so A->B and B->A collapse to one score row.
  key_el <- paste(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]), sep = "|")
  key_sc <- paste(pmin(res$scores$from, res$scores$to),
                  pmax(res$scores$from, res$scores$to), sep = "|")
  idx    <- match(key_el, key_sc)

  E(g)$score      <- res$scores$score[idx]        # confidence from the scaffold scoring
  E(g)$bidirected <- res$scores$bidirected[idx]   # TRUE where Gabi couldn't resolve one direction

  ## Spearman correlation between each edge's two endpoint expression columns;
  ## suppressWarnings() silences the tie warning for near-constant modules.
  rho <- mapply(function(a, b)
                  suppressWarnings(cor(expr[, a], expr[, b], method = "spearman")),
                el[, 1], el[, 2])
  E(g)$rho      <- as.numeric(rho)                          # kept for downstream QC
  E(g)$sign     <- ifelse(rho >= 0, "activation", "inhibition")  # rho >= 0 -> activation
  E(g)$sign_num <- ifelse(rho >= 0, 1L, -1L)                # +1/-1 for Boolean-rule arithmetic

  cat(sprintf("[%s] Edges: %d resolved, %d bidirected | %d activation, %d inhibition, %d unsignable\n",
              cl, sum(!E(g)$bidirected, na.rm = TRUE), sum(E(g)$bidirected, na.rm = TRUE),
              sum(rho >= 0, na.rm = TRUE), sum(rho < 0, na.rm = TRUE), sum(is.na(rho))))

  ## Drop the still-bidirected edges; keep every meta-node even if now isolated.
  backbone <- subgraph_from_edges(g, which(!E(g)$bidirected), delete.vertices = FALSE)

  graphml_full     <- file.path(cl_dir, sprintf(graphml_full_tmpl, cl, thr_str))
  graphml_backbone <- file.path(cl_dir, sprintf(graphml_backbone_tmpl, cl, thr_str))
  write_graph(g, graphml_full, format = "graphml")            # full signed net (QC)
  write_graph(backbone, graphml_backbone, format = "graphml") # resolved backbone
  cat(sprintf("[%s] Wrote:\n  %s\n  %s\n", cl, graphml_full, graphml_backbone))

  ## Skipped unless push_to_cytoscape is TRUE *and* Cytoscape opened successfully.
  if (push_to_cytoscape && cytoscape_ok) {
    cl_collection <- cytoscape_collection[[cl]]   # this cell line's own collection

    tryCatch({
      ## Push 'g' as 'title' into 'collection' - but only if no network of that
      ## title is already there. Additive only: a re-run must never clobber a
      ## network you've since laid out/analysed. To refresh one, delete it in
      ## Cytoscape first. getNetworkList() has no 'collection' filter, so the
      ## collection's own network names are looked up and matched instead.
      push_network <- function(g, title, collection) {
        collection_suid <- tryCatch(
          Find(function(s) identical(getCollectionName(s), collection),
               cyrestGET("collections")),
          error = function(e) NULL)

        if (!is.null(collection_suid)) {
          net_names <- vapply(getCollectionNetworks(collection_suid), getNetworkName, character(1))
          if (title %in% net_names) {
            cat(sprintf("  SKIP - '%s' already in Cytoscape collection '%s'\n", title, collection))
            return(invisible(NULL))
          }
        }
        createNetworkFromIgraph(g, title = title, collection = collection)
      }

      push_network(g, sprintf("%s_reactome_directed_signed", cl), cl_collection)
      push_network(backbone, sprintf("%s_reactome_backbone_signed", cl), cl_collection)

      ## Map 'sign' to an arrow shape so activation/inhibition reads visually.
      setEdgeTargetArrowShapeMapping(
        "sign",
        table.column.values = c("activation", "inhibition"),
        shapes              = c("ARROW", "T"),
        style.name          = "default")
    }, error = function(e)
         message(sprintf("[%s] RCy3 push skipped (%s) - import the GraphML files and map 'sign' to arrow shape.",
                         cl, conditionMessage(e))))
  }

  invisible(NULL)
}

# 3. RUN: every cell line ------------------------------------------------------

## Confirm Cytoscape is reachable once, before the loop. cytoscape_ok gates every
## push; FALSE just means "no Cytoscape, skip pushing". Networks go into whichever
## session is already open - this script never opens or saves a session file.
cytoscape_ok <- FALSE

if (push_to_cytoscape) {
  tryCatch({
    suppressMessages(library(RCy3))
    cytoscapePing()
    cytoscape_ok <- TRUE
  }, error = function(e)
       message(sprintf("Cytoscape push skipped for this run (%s) - import the GraphML files instead.",
                       conditionMessage(e))))
}

for (cl in cell_lines) {
  sign_export_cell_line(cl, thr_str, gabi_dir, result_tmpl, expr_tmpl,
                        graphml_full_tmpl, graphml_backbone_tmpl,
                        push_to_cytoscape, cytoscape_ok, cytoscape_collection)
}

cat("Done. Session left unsaved - save it yourself once checked.\n")

# END SCRIPT -------------------------------------------------------------------
