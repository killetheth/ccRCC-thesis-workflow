#!/usr/bin/env Rscript
## Purpose: Pull every cell line's analysed network back out of the running
## Cytoscape session into the *_cyto.graphml files 5-sort_modules.R expects.
##
## Other half of 3-sign_export_reactome.R's push step: that script creates
## "<CL>_reactome_directed_signed" / "<CL>_reactome_backbone_signed" networks via
## RCy3; once you've laid out/clustered/annotated them in Cytoscape, run this to
## export each to output/gabi/reactome/<CL>/<CL>_<type>_signed_d<thr>_cyto.graphml.
##
## Usage: Rscript 4-export_cyto_reactome.R   (run from gabi-maboss/scripts/reactome/,
## with Cytoscape open and the analysed session loaded)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/gabi-maboss/scripts/reactome")
getwd()   # sanity check (kept for posterity)

## Cytoscape integration; banners hidden.
suppressMessages(library(RCy3))

cytoscapePing()

# 1. INPUTS --------------------------------------------------------------------

## Threshold this session holds; matches 3-sign_export_reactome.R's 'thr'.
thr_str <- "0"

## Cell lines expected in the session; anything else found is skipped.
cell_lines <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")

## Per-cell-line Gabi output folder; %s = cell-line code.
gabi_dir <- "../../output/gabi/reactome/%s"

## Matches the titles 3-sign_export_reactome.R's push_network() creates, e.g.
## "AH_reactome_directed_signed" / "AH_reactome_backbone_signed". Restricted to
## the reactome branch since this script's gabi_dir (above) only routes there --
## a humannet/huri push into the same session would otherwise be misfiled into
## the reactome cell-line folder.
name_pattern <- "^([A-Z0-9]+)_reactome_(directed|backbone)_signed$"

# 2. RUN: every matching network in the session --------------------------------

networks <- getNetworkList()
exported <- 0L

for (net_name in networks) {
  m <- regmatches(net_name, regexec(name_pattern, net_name))[[1]]

  ## Not one of this pipeline's networks - leave it alone.
  if (length(m) == 0) next

  cl       <- m[2]
  net_type <- m[3]

  if (!cl %in% cell_lines) {
    warning(sprintf("[%s] Unrecognised cell-line code parsed from '%s' - skipping.", cl, net_name))
    next
  }

  cl_dir <- sprintf(gabi_dir, cl)

  if (!dir.exists(cl_dir)) {
    warning(sprintf("[%s] Expected output folder missing: %s - skipping '%s'.", cl, cl_dir, net_name))
    next
  }

  ## 5-sort_modules.R reads cell line from the parent folder and
  ## 'directed'/'backbone' + threshold from the file name itself.
  out_file <- file.path(cl_dir, sprintf("%s_reactome_%s_signed_d%s_cyto", cl, net_type, thr_str))

  tryCatch({
    ## exportNetwork()'s 'network' arg is not reliably honoured by Cytoscape's
    ## "network export" command - without first switching the active network,
    ## every call silently exports whatever was last active in the UI instead.
    setCurrentNetwork(net_name)
    exportNetwork(filename = out_file, type = "graphML", network = net_name, overwriteFile = TRUE)
    cat(sprintf("[%s] Exported '%s' -> %s.graphml\n", cl, net_name, out_file))
    exported <- exported + 1L
  }, error = function(e)
       message(sprintf("[%s] Export failed for '%s' (%s)", cl, net_name, conditionMessage(e))))
}

cat(sprintf("Done. %d network(s) exported.\n", exported))

# END SCRIPT -------------------------------------------------------------------
