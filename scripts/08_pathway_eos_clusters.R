#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## Per-cluster GO and KEGG enrichment on eosinophil cluster markers (from step 04)
## SLURM: --mem=64G --cpus-per-task=4 --time=2:00:00 -p batch

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
enrich_dir <- file.path(BASE, "enrichment")
fig_dir <- file.path(BASE, "figures")
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load cluster markers from step 04 ────────────────────────────────────────
cat("Loading eos_cluster_markers.csv...\n")
markers <- read.csv(file.path(BASE, "eos_cluster_markers.csv"), stringsAsFactors = FALSE)
cat(sprintf("  %d marker genes across %d clusters\n",
            nrow(markers), length(unique(markers$cluster))))

clusters <- sort(unique(markers$cluster))
cat(sprintf("  Clusters: %s\n", paste(clusters, collapse = ", ")))

# ── Per-cluster enrichment ────────────────────────────────────────────────────
all_go <- list()
all_kegg <- list()

for (clust in clusters) {
  cat(sprintf("\n=== Cluster %s ===\n", clust))

  deg_sub <- markers %>%
    filter(cluster == clust, p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    arrange(desc(avg_log2FC))
  cat(sprintf("  Significant markers: %d\n", nrow(deg_sub)))

  if (nrow(deg_sub) < 5) {
    cat("  Too few markers, skipping.\n")
    next
  }

  # Map gene symbols to Entrez IDs
  mapping <- bitr(deg_sub$gene, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  cat(sprintf("  Mapped %d / %d genes to Entrez IDs\n",
              nrow(mapping), nrow(deg_sub)))

  if (nrow(mapping) < 5) next

  gene_entrez <- unique(mapping$ENTREZID)

  # GO Biological Process
  cat("  Running GO:BP enrichment...\n")
  ego <- tryCatch({
    enrichGO(gene = gene_entrez,
             OrgDb = org.Mm.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2,
             readable = TRUE)
  }, error = function(e) { cat(sprintf("  GO error: %s\n", e$message)); NULL })

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    cat(sprintf("  GO:BP terms enriched: %d\n", nrow(as.data.frame(ego))))
    write.csv(as.data.frame(ego),
              file.path(enrich_dir, sprintf("go_bp_%s.csv", clust)),
              row.names = FALSE)
    all_go[[clust]] <- as.data.frame(ego) %>%
      mutate(cluster = clust) %>%
      head(10)
  } else {
    cat("  No significant GO:BP terms.\n")
  }

  # KEGG
  cat("  Running KEGG enrichment...\n")
  ekegg <- tryCatch({
    enrichKEGG(gene = gene_entrez,
               organism = "mmu",
               pvalueCutoff = 0.05)
  }, error = function(e) { cat(sprintf("  KEGG error: %s\n", e$message)); NULL })

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    cat(sprintf("  KEGG pathways enriched: %d\n", nrow(as.data.frame(ekegg))))
    write.csv(as.data.frame(ekegg),
              file.path(enrich_dir, sprintf("kegg_%s.csv", clust)),
              row.names = FALSE)
    all_kegg[[clust]] <- as.data.frame(ekegg) %>%
      mutate(cluster = clust) %>%
      head(10)
  } else {
    cat("  No significant KEGG pathways.\n")
  }

  # Per-cluster dotplot
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    pdf(file.path(fig_dir, sprintf("08_go_dotplot_%s.pdf", clust)),
        width = 10, height = 7)
    print(dotplot(ego, showCategory = 15) +
            ggtitle(sprintf("GO:BP Enrichment — %s", clust)) +
            theme(plot.title = element_text(face = "bold")))
    dev.off()
  }
}

# ── Combined GO heatmap (top 6 terms per cluster) ────────────────────────────
cat("\n=== Generating combined GO heatmap ===\n")

if (length(all_go) > 0) {
  combined_go <- bind_rows(all_go)

  # Select top 6 per cluster
  top_terms <- combined_go %>%
    group_by(cluster) %>%
    slice_min(order_by = p.adjust, n = 6) %>%
    ungroup()

  # Build -log10(padj) matrix for heatmap
  unique_terms <- unique(top_terms$Description)
  unique_clusters <- sort(unique(top_terms$cluster))

  heatmap_mat <- matrix(0, nrow = length(unique_terms), ncol = length(unique_clusters),
                        dimnames = list(unique_terms, unique_clusters))

  for (i in seq_len(nrow(top_terms))) {
    term <- top_terms$Description[i]
    clust <- top_terms$cluster[i]
    heatmap_mat[term, clust] <- -log10(top_terms$p.adjust[i])
  }

  # Truncate extreme values for visualization
  heatmap_mat[heatmap_mat > 20] <- 20

  # Order rows by cluster: C0 terms first, then C1, etc.
  row_order <- c()
  for (cl in unique_clusters) {
    cl_terms <- top_terms %>% filter(cluster == cl) %>%
      arrange(p.adjust) %>% pull(Description)
    # Only add terms not already in the list (avoid duplicates)
    cl_terms <- setdiff(cl_terms, row_order)
    row_order <- c(row_order, cl_terms)
  }
  # Add any remaining terms
  row_order <- c(row_order, setdiff(unique_terms, row_order))
  heatmap_mat <- heatmap_mat[row_order, , drop = FALSE]

  # Row annotation: which cluster each term belongs to
  term_cluster <- top_terms %>%
    group_by(Description) %>%
    slice_min(order_by = p.adjust, n = 1) %>%
    ungroup()
  row_anno <- data.frame(
    Cluster = term_cluster$cluster[match(row_order, term_cluster$Description)],
    row.names = row_order
  )
  anno_colors <- list(Cluster = eos_colors[unique_clusters])

  # Add gaps between cluster blocks
  gap_positions <- cumsum(table(row_anno$Cluster)[unique_clusters])
  gap_positions <- gap_positions[-length(gap_positions)]

  pdf(file.path(fig_dir, "08_go_heatmap_clusters.pdf"),
      width = 10, height = max(8, nrow(heatmap_mat) * 0.3))
  pheatmap(heatmap_mat,
           color = colorRampPalette(c("white", "#FDDBC7", "#D6604D", "#B2182B"))(50),
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_anno,
           annotation_colors = anno_colors,
           gaps_row = gap_positions,
           fontsize_row = 8, fontsize_col = 11,
           main = "GO:BP Enrichment by Eosinophil Cluster",
           angle_col = 0,
           border_color = "grey90")
  dev.off()
  cat("  Saved: 08_go_heatmap_clusters.pdf\n")

  # Save combined tables
  write.csv(combined_go,
            file.path(enrich_dir, "go_bp_all_clusters_top10.csv"),
            row.names = FALSE)
}

# Combined KEGG
if (length(all_kegg) > 0) {
  combined_kegg <- bind_rows(all_kegg)
  write.csv(combined_kegg,
            file.path(enrich_dir, "kegg_all_clusters_top10.csv"),
            row.names = FALSE)
  cat("  Saved: kegg_all_clusters_top10.csv\n")
}

cat("\nDone! 08_pathway_eos_clusters.R\n")
