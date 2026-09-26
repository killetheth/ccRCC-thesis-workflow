## Shared pipeline for the step-3 module builder, sourced by 3-cyto_module.R,
## which calls run_cyto_module() once per ID space (Entrez, HGNC) — the two
## calls differ only in their Cytoscape input folder and the gene-ID column
## emitted in *_genes_in_modules.csv. No side effects at source time.
##
## run_cyto_module() extracts Cluster_Number/Profiles from the Cytoscape tables,
## integrates the replicate expression data, decomposes each network by
## Cluster_Profile (taking the mean of the expression data), and writes the
## Gabi-ready module tables + mean-expression files. Nodes with an NA
## Cluster_Number keep their Entrez ID as the cluster.

suppressMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

## cytoscape_data_in     - folder of Cytoscape .csv exports (Cluster_Number, Profiles)
## genes_in_modules_col  - gene-ID column to emit in *_genes_in_modules.csv ("name" or "HGNC")
## genes_in_modules_out  - folder for the *_genes_in_modules.csv tables (ID-space specific)
run_cyto_module <- function(cytoscape_data_in, genes_in_modules_col, genes_in_modules_out) {

  ## Shared input/output locations (identical for both ID spaces).
  exp_data_in    <- "../exp_data/replicates/name_added/"
  expanded_out   <- "../output/cyto_module_output/gene_lists_expanded/"
  gene_lists_out <- "../output/cyto_module_output/gene_lists/"
  gabi_input_out <- "../output/cyto_module_output/gabi_input/"
  counts_out     <- "../output/cyto_module_output/cluster_profile_gene_counts.csv"

  for (d in c(genes_in_modules_out, expanded_out, gene_lists_out, gabi_input_out, dirname(counts_out)))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # 1. THE CYTOSCAPE DATA ------------------------------------------------------

  process_cyt <- function(fp) {
    read_csv(fp, show_col_types = FALSE) %>%
      select(any_of(c("name", "Cell_Line", "Multi_Profiles", "Cluster_Number", "HGNC")))
  }

  cyt_files    <- list.files(cytoscape_data_in, pattern = "\\.csv$", full.names = TRUE)
  cyt_combined <- bind_rows(lapply(cyt_files, process_cyt))

  cat(sum(is.na(cyt_combined$Cluster_Number)), "missing Cluster_Number\n")

  ## Fill a missing Cluster_Number with the gene's Entrez ID.
  cyt_combined <- cyt_combined %>%
    mutate(Cluster_Number = ifelse(is.na(Cluster_Number),
                                   as.character(name), as.character(Cluster_Number)))

  # 2. THE EXPRESSION DATA -----------------------------------------------------

  process_exp <- function(fp) {
    exp <- read_delim(fp, delim = "\t", show_col_types = FALSE)
    exp$Cell_Line <- str_extract(basename(fp), "^[A-Z]+")   # cell line/condition file prefix
    exp
  }

  exp_files    <- list.files(exp_data_in, full.names = TRUE)
  exp_combined <- bind_rows(lapply(exp_files, process_exp))

  ## Merge Cytoscape + expression (keep only genes from the Cytoscape files),
  ## then fill missing expression with 100 and confirm nothing else is NA.
  cyt_exp_combined <- cyt_combined %>%
    left_join(exp_combined, by = c("name", "Cell_Line")) %>%
    mutate(across(where(is.numeric), ~coalesce(., 100)))
  stopifnot(sum(is.na(cyt_exp_combined)) == 0)

  # 3. MODULE GENE LISTS -------------------------------------------------------

  genes_by_net <- split(cyt_exp_combined, cyt_exp_combined$Cell_Line)

  ## 3a. Expanded: one row per gene, with Cluster_Profile + Entrez + HGNC.
  for (network in names(genes_by_net)) {
    df <- genes_by_net[[network]] %>%
      mutate(Cluster_Profile = paste(as.character(Cluster_Number), Multi_Profiles, sep = ".")) %>%
      select(Cluster_Profile, name, HGNC) %>%
      arrange(Cluster_Profile)
    write_csv(df, file.path(expanded_out, paste0(network, "_module_gene_lists_expanded.csv")))
  }

  ## 3b. Collapsed: one row per module, genes comma-joined (HGNC).
  for (network in names(genes_by_net)) {
    df <- genes_by_net[[network]] %>%
      mutate(Cluster_Profile = paste(as.character(Cluster_Number), Multi_Profiles, sep = ".")) %>%
      select(Cluster_Profile, name, HGNC) %>%
      group_by(Cluster_Profile) %>%
      summarise(
        name            = paste(name, collapse = ", "),
        Genes_in_Module = paste(HGNC, collapse = ", "),
        HGNC_List       = list(HGNC),
        .groups         = "drop"
      ) %>%
      rowwise() %>%
      mutate(Node_ID = ifelse(length(unlist(HGNC_List)) > 1,
                              Cluster_Profile, unlist(HGNC_List)[1])) %>%
      ungroup() %>%
      select(-HGNC_List) %>%
      arrange(Cluster_Profile)
    write_csv(df, file.path(gene_lists_out, paste0(network, "_module_gene_lists.csv")))
  }

  ## 3c. Gabi-ready module table: the ID column ('name' or 'HGNC') + Cluster_Profile.
  for (network in names(genes_by_net)) {
    df <- genes_by_net[[network]] %>%
      mutate(Cluster_Profile = paste(as.character(Cluster_Number), Multi_Profiles, sep = ".")) %>%
      select(all_of(genes_in_modules_col), Cluster_Number, Multi_Profiles, Cluster_Profile) %>%
      arrange(Cluster_Profile)
    write_csv(df, file.path(genes_in_modules_out, paste0(network, "_genes_in_modules.csv")))
  }

  ## Gene counts per Cluster_Profile (>= 5 genes), for interest.
  profile_gene_counts <- cyt_exp_combined %>%
    mutate(Cluster_Profile = paste(as.character(Cluster_Number), Multi_Profiles, sep = ".")) %>%
    group_by(Cell_Line, Cluster_Profile) %>%
    summarise(Gene_Count = n_distinct(name), .groups = "drop") %>%
    filter(Gene_Count >= 5)
  write_csv(profile_gene_counts, counts_out)

  # 4. AGGREGATE EXPRESSION (mean per module) ----------------------------------

  final_df <- cyt_exp_combined %>%
    mutate(Cluster_Profile = paste(as.character(Cluster_Number), Multi_Profiles, sep = ".")) %>%
    group_by(Cell_Line, Cluster_Profile) %>%
    summarise(across(starts_with(c("0h", "1h", "2h", "4h", "8h", "24h")),
                     ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    relocate(Cell_Line, Cluster_Profile)

  # 5. CLEAN UP ----------------------------------------------------------------

  ## Strip the 'h' from timepoint column names for Gabi, then confirm no NAs.
  names(final_df) <- gsub("([0-9]+)h\\.?", "\\1.", names(final_df))
  stopifnot(sum(is.na(final_df)) == 0)

  ## Order timepoint columns numerically.
  timepoint_cols     <- setdiff(names(final_df), c("Cell_Line", "Cluster_Profile"))
  ordered_timepoints <- timepoint_cols[order(as.numeric(sub("\\..*", "", timepoint_cols)))]
  final_df <- final_df %>% select(Cell_Line, Cluster_Profile, all_of(ordered_timepoints))

  # 6. SPLIT INTO 8 NETWORKS ---------------------------------------------------

  split_by_network <- split(final_df, final_df$Cell_Line)
  for (network in names(split_by_network)) {
    df_network <- split_by_network[[network]] %>% select(-Cell_Line)
    write_csv(df_network, file.path(gabi_input_out, paste0(network, "_cluster_profile_mean_exp.csv")))
  }

  cat("Done.\n")
}
