#!/usr/bin/env Rscript
## Purpose: Build undirected meta-node scaffolds for Gabi (1.0.1) from the Human
## Reference Interactome (HuRI), in batch. Runs over every cell line, writing one
## GraphML scaffold per cell line at the d0 (no-threshold) cutoff for Gabi to
## orient.
##
## HuRI/HGNC counterpart of ../reactome/1-undir_scaff_reactome.R: gene-level
## edges come from HI-union + Lit-BM (unioned for the largest interaction set)
## instead of Reactome FI, and those files ship as Ensembl gene ID (ENSG) pairs,
## so this script maps them to HGNC symbols before joining -annotation is then
## joined on HGNC, same as Reactome. The shared contraction + density logic
## lives in ../scaffold_common.R.
##
## Pipeline (per cell line): HI-union + Lit-BM edges (read once, unioned, mapped
## Ensembl -> HGNC) -> keep genes present in the annotation (matched by HGNC) ->
## label each with its Cluster_Profile -> contract genes into modules, counting
## gene-gene edges between each module pair -> keep a module edge A-B only if
## n_edges(A,B) / (|A|*|B|) >= threshold -> write one plain undirected GraphML.
##
## Usage: Rscript 1-undir_scaff_huri.R   (run from gabi-maboss/scripts/huri/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/gabi-maboss/scripts/huri")
getwd()   # sanity check (kept for posterity)

## Graph + GraphML I/O; ID mapping; banners hidden.
suppressMessages({ library(igraph); library(AnnotationDbi); library(org.Hs.eg.db) })

## Shared build_metagraph() + logmsg().
source("../scaffold_common.R")

# 1. INPUTS --------------------------------------------------------------------

## 1a. HuRI gene-gene edge lists, unioned for the largest interaction set:
##     header-less, tab-separated; columns 1-2 are the two Ensembl gene IDs
##     (ENSG) of an undirected edge. Read as character so IDs match exactly.
##     Shared, so read + unioned once (section 2).
hi_union_file <- "../../data/HI-union.tsv"
lit_bm_file   <- "../../data/Lit-BM.tsv"

## 1b. Cell lines to process (each needs a matching annotation file below).
cell_lines <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")

## 1c. Per-cell-line module table (cyto_module HGNC output); %s = cell-line code.
##     Columns: HGNC (join key, post Ensembl->HGNC mapping), Cluster_Number,
##     Multi_Profiles, Cluster_Profile (module label).
annot_tmpl   <- "../../data/cyto_module_output/hgnc/%s_genes_in_modules.csv"
col_gene     <- "HGNC"             # HGNC symbol column (join key, post-mapping)
col_metanode <- "Cluster_Profile"  # Module label column, read verbatim

## 1d. Edge-density threshold: keep a scaffold edge A-B only if the fraction of
##     the |A|*|B| possible gene-gene pairs that are actually linked reaches thr.
##     HuRI starts with thr = 0 (no threshold), same as Reactome: HI-union +
##     Lit-BM is already a sparse, curated binary interaction set.
thr     <- 0
thr_str <- formatC(thr, format = "g")   # e.g. "0"; avoids scientific notation

## 1e. Output. One flat dir; the cell-line code + threshold are carried in the
##     filename. Kept separate from the HumanNet/Reactome scaffolds so none
##     clobbers another. %s slots in out_tmpl are the cell-line code then the
##     threshold.
out_dir  <- "../../output/scaffold/huri"
out_tmpl <- "%s_huri_metanode_scaffold_d%s.graphml"

## 1f. Logs. One plain-text log per cell line (gene/module/edge counts, density
##     quantiles, per-threshold edge/isolated counts). %s = cell-line code.
log_dir  <- "../../logs/scaffold/huri"
log_tmpl <- "%s_scaffold.log"

# 2. SHARED SETUP --------------------------------------------------------------

if (!file.exists(hi_union_file))
  stop("HI-union file not found: ", hi_union_file,
       " -download it from the HuRI download page and place it there.")
if (!file.exists(lit_bm_file))
  stop("Lit-BM file not found: ", lit_bm_file,
       " -download it from the HuRI download page and place it there.")

## Read both edge lists as strings so Ensembl IDs match exactly; quote/comment
## disabled so stray "/# chars can't break parsing.
hi <- read.table(hi_union_file, sep = "\t", header = FALSE,
                 colClasses = "character", quote = "", comment.char = "")
lb <- read.table(lit_bm_file, sep = "\t", header = FALSE,
                 colClasses = "character", quote = "", comment.char = "")

## Union of both sources for the largest interaction set; simplify() dedupes
## edges appearing in both files and drops any self-loops.
g_ensg <- simplify(graph_from_data_frame(rbind(hi[, 1:2], lb[, 1:2]), directed = FALSE))
cat(sprintf("HI-union: %d edges | Lit-BM: %d edges | union (deduped): %d genes, %d edges\n",
            nrow(hi), nrow(lb), vcount(g_ensg), ecount(g_ensg)))

## Map Ensembl gene IDs -> HGNC symbols so HuRI joins onto the same hgnc/ module
## annotation the Reactome branch uses. Genes with no symbol are dropped;
## distinct Ensembl IDs that collapse onto the same symbol are merged by the
## simplify() below (parallel edges -> one, any resulting self-loops removed).
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = V(g_ensg)$name, column = "SYMBOL",
                               keytype = "ENSEMBL", multiVals = "first"))
mapped      <- !is.na(sym)
g_gene_full <- induced_subgraph(g_ensg, which(mapped))
V(g_gene_full)$name <- sym[mapped]
g_gene_full <- simplify(g_gene_full, remove.multiple = TRUE, remove.loops = TRUE)
cat(sprintf("Ensembl -> HGNC: %d/%d genes mapped -> %d genes, %d edges after merge\n",
            sum(mapped), vcount(g_ensg), vcount(g_gene_full), ecount(g_gene_full)))

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

## Downstream note: 2-gabi_huri.R wants expression as SAMPLES-in-rows,
## NODES-in-columns; the '*_cluster_profile_mean_exp.csv' file is the transpose,
## and that reshape now happens inside 2-gabi_huri.R. Its row labels are the
## same Cluster_Profile values used here, so node names line up by construction.

# END SCRIPT -------------------------------------------------------------------
