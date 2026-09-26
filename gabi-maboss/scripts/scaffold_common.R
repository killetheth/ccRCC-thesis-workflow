## Shared helpers for the step-1 scaffold builders, sourced by both
## 1-undir_scaff_humannet.R and 1-undir_scaff_reactome.R (which differ only in
## their gene-level edge source and join key). No side effects at source time.

## logmsg: print 'fmt' (sprintf-style) to the console AND append it to the
## current cell line's open log connection 'con'. Used in place of cat()
## throughout the per-cell-line loop.
logmsg <- function(con, fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(msg)
  cat(msg, file = con)
}

## build_metagraph: contract a cell line's genes into a weighted module graph
## whose edges carry 'density' = n_edges(A,B) / (|A| * |B|). The contraction is
## done once here; the threshold sweep in the caller just filters this graph's
## edges, so it isn't recomputed per cutoff. 'col_gene' is the join key column
## (Entrez for HumanNet, HGNC for Reactome); 'col_metanode' is the module label.
build_metagraph <- function(g_gene_full, annot_file, col_gene, col_metanode) {
  ## Read the annotation table (e.g. AH_genes_in_modules.csv); IDs/labels kept
  ## as strings so they match the gene network exactly, headers kept verbatim.
  annot <- read.csv(annot_file, check.names = FALSE,
                    colClasses = "character", stringsAsFactors = FALSE)
  gene     <- as.character(annot[[col_gene]])      # per-row gene ID/symbol
  metanode <- as.character(annot[[col_metanode]])  # per-row module label (aligned)

  ## Keep genes present in BOTH network and annotation - the join that links the
  ## gene network to the cell line's module labels - then prune to just those.
  keep   <- intersect(V(g_gene_full)$name, gene)
  g_gene <- induced_subgraph(g_gene_full, keep)

  ## Empty overlap = mismatched ID spaces -> fail loudly.
  if (vcount(g_gene) == 0L)
    stop("No overlap between the gene network IDs and the '", col_gene,
         "' column - check both use the same ID space.")

  ## Module label per surviving vertex, realigned from CSV row order to g_gene's
  ## vertex order ('gene'/'metanode' are parallel CSV columns).
  mn        <- metanode[match(V(g_gene)$name, gene)]
  ## Distinct labels fix the group numbering and the contracted vertex order.
  mn_levels <- unique(mn)
  ## contract() needs an integer group ID per vertex: each label's position.
  mapping   <- match(mn, mn_levels)
  ## |A| for every module A, one entry per group in mn_levels order.
  mn_size   <- tabulate(mapping, nbins = length(mn_levels))

  ## Collapse same-group vertices into one super-vertex per module; edges become
  ## parallel edges / self-loops, deduped by simplify() below.
  g_meta <- contract(g_gene, mapping)
  ## Reattach the Cluster_Profile label to each module (contract() numbers new
  ## vertices by group ID, and mn_levels[k] is group k's label).
  V(g_meta)$name <- mn_levels

  ## Weight 1 per edge, summed on collapse, so each surviving edge's weight =
  ## the gene-gene edge count between that module pair (loops removed first).
  E(g_meta)$weight <- 1
  g_meta <- simplify(g_meta, remove.multiple = TRUE, remove.loops = TRUE,
                     edge.attr.comb = list(weight = "sum"))

  ## density = actual cross-links / max possible (|A| * |B|) per module edge.
  ep <- ends(g_meta, E(g_meta), names = FALSE)
  E(g_meta)$density <- E(g_meta)$weight / (mn_size[ep[, 1]] * mn_size[ep[, 2]])

  list(g_meta = g_meta, n_genes = vcount(g_gene))
}
