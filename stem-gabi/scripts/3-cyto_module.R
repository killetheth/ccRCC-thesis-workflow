#!/usr/bin/env Rscript
## Purpose: Build the Gabi-ready module tables + mean-expression files from the
## Cytoscape exports and replicate expression data, keying the module gene tables
## by both Entrez ID (feeds the HumanNet branch) and HGNC symbol (feeds the
## Reactome/HuRI branches). Runs both ID spaces in one invocation, both reading
## the same for_cytomod/ export (one row per gene, with both a Entrez `name`
## and a current, non-duplicated `HGNC` column).
##
## Shared pipeline lives in cyto_module_common.R.
##
## Usage: Rscript 3-cyto_module.R   (run from stem-gabi/scripts/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/stem-gabi/scripts")
getwd()   # sanity check (kept for posterity)

## Shared run_cyto_module() (also attaches the tidyverse libraries).
source("cyto_module_common.R")

# 1. RUN: both ID spaces --------------------------------------------------------

run_cyto_module(
  cytoscape_data_in    = "../output/cytoscape_exp/for_cytomod/",
  genes_in_modules_col = "name",   # Entrez ID (join key for the HumanNet branch)
  genes_in_modules_out = "../output/cyto_module_output/entrez/"
)

run_cyto_module(
  cytoscape_data_in    = "../output/cytoscape_exp/for_cytomod/",
  genes_in_modules_col = "HGNC",   # HGNC symbol (join key for the Reactome/HuRI branches)
  genes_in_modules_out = "../output/cyto_module_output/hgnc/"
)

# END SCRIPT -------------------------------------------------------------------
