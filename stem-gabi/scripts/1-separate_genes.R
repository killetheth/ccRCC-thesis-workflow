#!/usr/bin/env Rscript
## Purpose: Separate STEM's per-cell-line gene tables into profile-specific .csv
## files, one file per STEM profile per cell line.
##
## Usage: Rscript 1-separate_genes.R   (run from stem-gabi/scripts/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/stem-gabi/scripts")
getwd()   # sanity check (kept for posterity)

## Libraries.
suppressMessages({
  library(dplyr)
  library(tidyr)
})

# 1. INPUTS --------------------------------------------------------------------

cell_lines    <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")
stem_in_dir   <- "../output/stem_outputs/stem_main_genes"
genes_out_dir <- "../output/genes_by_profile"

# 2. FUNCTION: PROCESS AND SAVE PROFILES ---------------------------------------

process_profiles <- function(cell_line) {
  file_path <- file.path(stem_in_dir, paste0(cell_line, "_table.txt"))
  gene_data <- read.delim(file_path, sep = "\t", header = TRUE)

  ## One row per profile (split any multi-profile genes), profile as numeric.
  split_profiles <- gene_data %>%
    separate_rows(Profile, sep = ";") %>%
    mutate(Profile = as.numeric(Profile))

  ## Split into a list of data frames, one per profile.
  profiles_list <- split(split_profiles, split_profiles$Profile)

  cell_genes_out_dir <- file.path(genes_out_dir, cell_line)
  dir.create(cell_genes_out_dir, recursive = TRUE, showWarnings = FALSE)

  ## Save each profile to its own .csv.
  lapply(names(profiles_list), function(profile) {
    output_file <- file.path(cell_genes_out_dir, paste0(cell_line, "_profile_", profile, ".csv"))
    write.csv(profiles_list[[profile]], output_file, row.names = FALSE)
  })

  message("Processed: ", cell_line)
}

# 3. RUN: every cell line ------------------------------------------------------

lapply(cell_lines, process_profiles)

# END SCRIPT -------------------------------------------------------------------
