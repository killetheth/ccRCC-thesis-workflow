#!/usr/bin/env Rscript
## Purpose: Count significant genes across all cell line/conditions, build a
## Cytoscape-compatible gene-profile file per cell line, and count unique genes.
##
## Note: this is reliant on the profiles written by 1-separate_genes.R to be
## first manually sorted into significant/ and not_significant/ subfolders per
## cell line.
##
## Reads the profile CSVs written by 1-separate_genes.R after they've been sorted
## into significant/ and not_significant/ subfolders per cell line.
##
## Usage: Rscript 2-genes_in_profile.R   (run from stem-gabi/scripts/)

# 0. SETUP ---------------------------------------------------------------------

## Clean slate.
rm(list = ls())

## Pin the working directory to this script's own folder; every path below is
## relative to it, so the script runs the same from anywhere on the Windows
## workstation.
setwd("C:/ccRCC-causal-networks/stem-gabi/scripts")
getwd()   # sanity check (kept for posterity)

## Libraries.
suppressMessages(library(dplyr))

# 1. INPUTS --------------------------------------------------------------------

cell_lines <- c("AH", "AN", "CH", "CN", "HH", "HN", "UH", "UN")

## Significant / not-significant profile folders (one per cell line).
genes_by_profile <- "../output/genes_by_profile"
sig_folders    <- setNames(file.path(genes_by_profile, cell_lines, "significant"), cell_lines)
notsig_folders <- setNames(file.path(genes_by_profile, cell_lines, "not_significant"), cell_lines)

## STEM output gene tables (the original per-cell-line gene lists).
gene_files <- setNames(
  file.path("../output/stem_outputs/stem_main_genes", paste0(cell_lines, "_table.txt")),
  cell_lines)

## Output folders.
cyto_profile_dir <- "../output/genes-profile_for_cytoscape"
gene_count_dir   <- "../output/gene_count"

# 2. STATISTICS AFTER STEM PROFILING -------------------------------------------

## All genes originally fed in, across the 8 cell lines.
original_gene_list <- unique(unlist(lapply(gene_files, function(f) read.delim(f, sep = "\t")$Gene)))

## All significant genes from the profile CSVs.
all_sig_genes <- unlist(lapply(sig_folders, function(folder) {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  unlist(lapply(files, function(f) read.csv(f)$Gene))
}))

unique_sig_genes <- unique(all_sig_genes)
sig_in_original  <- intersect(unique_sig_genes, original_gene_list)

cat("Total genes fed in:", length(original_gene_list), "\n")
cat("Significant genes found:", length(sig_in_original), "\n")
cat("Percentage covered:", round(100 * length(sig_in_original) / length(original_gene_list), 2), "%\n\n")

## Same breakdown per cell line.
for (cell_line in cell_lines) {
  original_genes <- unique(read.delim(gene_files[[cell_line]], sep = "\t")$Gene)

  profile_files <- list.files(sig_folders[[cell_line]], pattern = "\\.csv$", full.names = TRUE)
  sig_genes     <- unique(unlist(lapply(profile_files, function(f) read.csv(f)$Gene)))
  overlap_genes <- intersect(original_genes, sig_genes)

  cat("Cell line:", cell_line, "\n")
  cat("  Total genes originally fed in:", length(original_genes), "\n")
  cat("  Sig. genes in profiles:", length(sig_genes), "\n")
  cat("  Sig. genes in original set:", length(overlap_genes), "\n")
  cat("  Sig. genes %:", round(100 * length(overlap_genes) / length(original_genes), 2), "%\n\n")
}

# 3. CYTOSCAPE-COMPATIBLE PROFILE FILE PER CELL LINE ---------------------------

dir.create(cyto_profile_dir, showWarnings = FALSE, recursive = TRUE)

for (cell_line in cell_lines) {
  profile_files <- list.files(sig_folders[[cell_line]], pattern = "profile_.*\\.csv", full.names = TRUE)

  ## Keep the gene, its profile number, and the cell-line code.
  profile_data <- lapply(profile_files, function(file_path) {
    profile_number <- gsub(".*_profile_(\\d+)\\.csv", "\\1", basename(file_path))
    read.csv(file_path) %>%
      select(Gene) %>%
      mutate(Profile = profile_number, NEWCellLine = cell_line)
  })

  ## Collapse to one row per gene, listing all its profiles.
  combined_profiles <- bind_rows(profile_data) %>%
    group_by(Gene) %>%
    summarise(
      Profiles    = paste(sort(unique(Profile)), collapse = ", "),
      NEWCellLine = first(NEWCellLine),
      .groups     = "drop"
    )

  write.csv(combined_profiles,
            file.path(cyto_profile_dir, paste0("Combined_", cell_line, "_Profiles.csv")),
            row.names = FALSE)
}

# 4. COUNT UNIQUE GENES (significant + not significant) ------------------------

## Files contain the gene list for each case; the count is carried in the name.
gene_count <- function(folder_list, folder_type) {
  for (cell_line in names(folder_list)) {
    gene_files    <- list.files(folder_list[[cell_line]], pattern = "\\.csv$", full.names = TRUE)
    unique_genes  <- unique(unlist(lapply(gene_files, function(f) read.csv(f)$Gene)))

    out_dir <- file.path(gene_count_dir, folder_type)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_file <- file.path(out_dir, paste0(cell_line, "_", folder_type, "_unique_genes_", length(unique_genes), ".csv"))
    write.table(unique_genes, out_file, row.names = FALSE, col.names = FALSE, sep = ",")

    message(length(unique_genes), " unique genes in ", cell_line, " (", folder_type, ")")
  }
}

gene_count(sig_folders, "significant")
gene_count(notsig_folders, "not_significant")

# END SCRIPT -------------------------------------------------------------------
