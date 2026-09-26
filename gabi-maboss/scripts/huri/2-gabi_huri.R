#!/usr/bin/env Rscript
## Purpose: Orient one cell line's HuRI meta-node scaffold with ClNetScOn
## (GABI 1.0.1), following demo_directed_from_scaffold.R. Cell line and density
## threshold are taken as arguments, so this one script replaces per-cell-line
## copies; run_all_gabi_huri.sh drives it over all eight.
##
## HuRI/HGNC counterpart of ../reactome/2-gabi_reactome.R. It differs only in
## which scaffold it orients (the HuRI scaffold instead of the Reactome one),
## and where outputs land (a huri/ subtree, so no branch clobbers another). The
## default threshold (d0, no-threshold) matches Reactome's. The expression input
## is shared: the same '*_cluster_profile_mean_exp.csv' labels line up with the
## HuRI scaffold node names by construction (both are Cluster_Profile labels).
##
## NB - A SILENT CONSEQUENCE OF THE 2-COLUMN TSV EDGE LIST WRITTEN IN SECTION 2.
## It has no third column, so ClNetScOn's load.scaffold() builds the graph from
## `el[, 1:2]` and the scaffold igraph ends up with NO `weight` edge attribute at
## all. That matters because ClNetScOn would otherwise SIGN these edges for us:
## gabi.directed() finishes by giving every scaffold-only edge a Spearman rho as
## its `weight` (tagging it test_stat_name = "Rho"). But it FINDS those edges with
##
##     scaffold.edges <- which(is.na(E(net.directed)$weight) |
##                             E(net.directed)$weight == -2)
##
## i.e. by looking for a sentinel value in an attribute it assumes exists. With no
## attribute at all, E(...)$weight is NULL, is.na(NULL) is logical(0), and which()
## returns integer(0) - so the block sees zero scaffold edges and skips. Its own
## "score.matrix is NULL" warning sits inside the same length(scaffold.edges) > 0
## test, so nothing is printed either: it fails silently, and the directed .gml
## comes back with no weight and no sign.
##
## That is why 3-sign_export_<branch>.R computes rho itself. The two are the same
## quantity - Spearman between the same two modules over the same 18 samples
## (tma() does no replicate aggregation here, because the module names are unique
## so hasreplicates = 0) - so the signs assigned downstream are the ones ClNetScOn
## would have assigned; nothing is wrong with the numbers. But if the scaffold is
## ever handed over WITH a weight column (or with NA/-2 sentinels), that block
## will fire and edges will carry a rho from both routes: check for double-signing
## before trusting either.
##
## Usage: Rscript 2-gabi_huri.R <CELL_LINE> [THRESHOLD]   (runs on the Linux
##        box; THRESHOLD defaults to 0)

# 0. SETUP ---------------------------------------------------------------------

## Graph I/O + ClNetScOn; banners hidden.
suppressMessages({ library(ClNetScOn); library(igraph) })

## Pin the working directory to this script's own folder, derived from its own
## location, so the heavy orient runs on any machine/OS without a hardcoded
## path. Rscript exposes the path via --file; no RStudio dependency.
.file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(.file)) stop("Run this script with Rscript (needs --file to locate its own folder).")
setwd(normalizePath(dirname(.file[1])))
getwd()   # sanity check (kept for posterity)

# 1. INPUTS --------------------------------------------------------------------

## Cell-line code (required) and density threshold (optional; HuRI uses 0).
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript 2-gabi_huri.R <CELL_LINE> [THRESHOLD]")
cl      <- args[1]
thr     <- if (length(args) >= 2) as.numeric(args[2]) else 0
thr_str <- formatC(thr, format = "g")   # matches 1-undir_scaff_huri.R's d<thr> suffix

meanexp_file  <- sprintf("../../data/cyto_module_output/exp_data/%s_cluster_profile_mean_exp.csv", cl)
scaffold_file <- sprintf("../../output/scaffold/huri/%s_huri_metanode_scaffold_d%s.graphml", cl, thr_str)

out_dir       <- sprintf("../../output/gabi/huri/%s", cl)
expr_out      <- file.path(out_dir, sprintf("%s_huri_expr_gabi.tsv", cl))
scaffold_tsv  <- file.path(out_dir, sprintf("%s_huri_metanode_scaffold_d%s.tsv", cl, thr_str))
directed_out  <- file.path(out_dir, sprintf("%s_huri_directed_net_d%s.gml", cl, thr_str))
result_out    <- file.path(out_dir, sprintf("%s_huri_gabi_result_d%s.rds", cl, thr_str))

mwq_solver    <- "qclique"   # "qclique" (fast) or "cliquer" (exact)
cores         <- 12

## Log file: mirrors everything printed below (including gabi.scaffold.master's
## own verbose = TRUE trace) to a per-cell-line log, so a long run's progress
## survives even if the console scrolls away or the session is detached.
log_dir  <- sprintf("../../logs/gabi/huri/%s", cl)
log_file <- file.path(log_dir, sprintf("%s_gabi_d%s.log", cl, thr_str))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

## split = TRUE keeps echoing cat()/print() to the console as well as the log;
## message()/warning() (type = "message") can't split, so those go to the log
## only. Reset via on.exit so a later error still leaves stdout/stderr sane.
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

# 2. SCAFFOLD -> 2-COLUMN TSV EDGE LIST ----------------------------------------

## (+ node set for the coverage check.)
g_sc     <- read_graph(scaffold_file, format = "graphml")
sc_nodes <- V(g_sc)$name
write.table(as_edgelist(g_sc, names = TRUE), scaffold_tsv,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("[%s d%s] Scaffold edge list: %d edges -> %s\n", cl, thr_str, ecount(g_sc), scaffold_tsv))

# 3. EXPRESSION -> GABI LAYOUT -------------------------------------------------

## Samples in rows, meta-nodes in columns.
me   <- read.csv(meanexp_file, check.names = FALSE, stringsAsFactors = FALSE)
labs <- me[["Cluster_Profile"]]
expr <- me[, setdiff(names(me), "Cluster_Profile"), drop = FALSE]   # nodes x samples
rownames(expr) <- labs
mat  <- t(as.matrix(expr))                                          # samples x nodes

missing <- setdiff(sc_nodes, colnames(mat))
if (length(missing))
  stop(length(missing), " scaffold node(s) lack an expression column, e.g. ",
       paste(head(missing), collapse = ", "))
mat <- mat[, intersect(colnames(mat), sc_nodes), drop = FALSE]

write.table(mat, expr_out, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("[%s] Expression: %d samples x %d meta-nodes -> %s\n",
            cl, nrow(mat), ncol(mat), expr_out))

# 4. ORIENT --------------------------------------------------------------------

## The demo's call, now with a TSV scaffold.
result <- gabi.scaffold.master(
  data       = expr_out,
  scaffold   = scaffold_tsv,
  mwq.solver = mwq_solver,
  cores      = cores,
  verbose    = TRUE               # watch where the time goes
)

directed <- result$directed.net

# 5. INSPECT + SAVE ------------------------------------------------------------

cat(sprintf("[%s d%s] Directed network: %d nodes, %d directed edges\n",
            cl, thr_str, vcount(directed), ecount(directed)))
cat("Edge attributes on directed.net: ",
    paste(edge_attr_names(directed), collapse = ", "), "\n")

save.graph(directed, directed_out)
saveRDS(result, result_out)
cat("Directed net -> ", directed_out, "\nFull result -> ", result_out, "\n", sep = "")

# END SCRIPT -------------------------------------------------------------------
