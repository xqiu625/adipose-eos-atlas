#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 32_fig2J_pseudotime_heatmap.R
# Fig 2J: Gene expression dynamics along Monocle3 pseudotime trajectory
# Heatmap of top 1,000 HVGs ordered by pseudotime (BM + Blood + C0-C5)
#
# Input:  $BASE/monocle_threeway/monocle3_threeway_cds.rds
# Output: revision/figures/fig2/fig2J_pseudotime_heatmap.pdf
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(pheatmap)
  library(viridis)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
PNAS <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision/figures/fig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Group colors (for annotation) -------------------------------------------
group_colors <- c(
  "Gurtner BM"    = "#2E8B57",
  "C3"            = "#CBAACB",
  "Gurtner Blood" = "#DAA520",
  "C1"            = "#E84A5F",
  "C5"            = "#5BA4A4",
  "C4"            = "#B03A48",
  "C0"            = "#6CA6CD",
  "C2"            = "#FF6F1F"
)

# =============================================================================
# 1. Load Monocle3 three-way CDS (as SCE — no monocle3 package needed)
# =============================================================================
cat("Loading monocle3_threeway_cds.rds ...\n")
cds <- readRDS(file.path(BASE, "monocle_threeway", "monocle3_threeway_cds.rds"))
cat(sprintf("  %d cells x %d genes\n", ncol(cds), nrow(cds)))

# Extract pseudotime directly from colData (monocle3 stores it there)
pt_all <- colData(cds)$pseudotime
if (is.null(pt_all)) {
  # Fallback: try monocle3_pseudotime slot name
  pt_all <- cds@colData$monocle3_pseudotime
}
names(pt_all) <- colnames(cds)

# =============================================================================
# 2. Subset to main figure groups (BM + Blood + C0-C5)
# =============================================================================
cat("Subsetting to main groups ...\n")
main_groups <- c("Gurtner BM", "Gurtner Blood",
                 "C0", "C1", "C2", "C3", "C4", "C5")
grp <- colData(cds)$main_group
keep <- !is.na(grp) & grp %in% main_groups & is.finite(pt_all)
cds_main <- cds[, keep]
pt <- pt_all[keep]
cat(sprintf("  Main group cells with finite pseudotime: %d\n", ncol(cds_main)))

# =============================================================================
# 3. Identify top 1,000 HVGs
# =============================================================================
cat("Identifying top 1,000 HVGs ...\n")

# Get normalized counts (log-normalized)
counts_mat <- as.matrix(assay(cds_main))
# Normalize: library-size normalize + log1p
lib_sizes <- colSums(counts_mat)
norm_mat <- log1p(t(t(counts_mat) / lib_sizes * 1e4))

# Compute variance of each gene across cells
gene_vars <- apply(norm_mat, 1, var)
gene_means <- rowMeans(norm_mat)

# Filter: keep genes expressed in >= 1% of cells
pct_expressed <- rowMeans(counts_mat > 0)
expressed_genes <- names(pct_expressed[pct_expressed >= 0.01])
cat(sprintf("  Genes expressed in >= 1%% of cells: %d\n", length(expressed_genes)))

# Select top 1000 by variance among expressed genes
gene_vars_filt <- gene_vars[expressed_genes]
top1000 <- names(sort(gene_vars_filt, decreasing = TRUE))[1:min(1000, length(gene_vars_filt))]
cat(sprintf("  Selected %d HVGs\n", length(top1000)))

# =============================================================================
# 4. Order cells by pseudotime and bin for visualization
# =============================================================================
cat("Ordering cells and binning ...\n")

# Order cells by pseudotime
cell_order <- order(pt)
pt_ordered <- pt[cell_order]
groups_ordered <- as.character(colData(cds_main)$main_group)[cell_order]

# Bin cells into 200 bins along pseudotime for cleaner visualization
n_bins <- 200
bin_idx <- cut(seq_along(pt_ordered), breaks = n_bins, labels = FALSE)

# Average expression per bin
norm_ordered <- norm_mat[top1000, cell_order]
bin_mat <- sapply(1:n_bins, function(b) {
  cells_in_bin <- which(bin_idx == b)
  if (length(cells_in_bin) == 1) {
    norm_ordered[, cells_in_bin]
  } else {
    rowMeans(norm_ordered[, cells_in_bin, drop = FALSE])
  }
})
colnames(bin_mat) <- paste0("bin_", 1:n_bins)

# Bin pseudotime (mean per bin) and majority group per bin
bin_pt <- tapply(pt_ordered, bin_idx, mean)
bin_group <- tapply(groups_ordered, bin_idx, function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
})

cat(sprintf("  Binned into %d pseudotime bins\n", n_bins))

# =============================================================================
# 5. Scale genes (z-score across bins)
# =============================================================================
cat("Scaling gene expression ...\n")
bin_scaled <- t(scale(t(bin_mat)))
bin_scaled[bin_scaled > 2.5] <- 2.5
bin_scaled[bin_scaled < -2.5] <- -2.5
bin_scaled[is.na(bin_scaled)] <- 0

# =============================================================================
# 6. Cluster genes by expression pattern (k-means)
# =============================================================================
cat("Clustering genes ...\n")
set.seed(42)
k <- 6  # number of gene clusters
km <- kmeans(bin_scaled, centers = k, nstart = 25, iter.max = 100)

# Order gene clusters by peak pseudotime position
cluster_peak <- sapply(1:k, function(cl) {
  genes_cl <- names(km$cluster[km$cluster == cl])
  mean_expr <- colMeans(bin_scaled[genes_cl, , drop = FALSE])
  which.max(mean_expr)
})
cluster_order <- order(cluster_peak)

# Reorder genes: within each cluster, sort by peak position
gene_order <- c()
gene_clusters <- c()
for (cl in cluster_order) {
  genes_cl <- names(km$cluster[km$cluster == cl])
  peaks <- apply(bin_scaled[genes_cl, , drop = FALSE], 1, which.max)
  genes_sorted <- genes_cl[order(peaks)]
  gene_order <- c(gene_order, genes_sorted)
  gene_clusters <- c(gene_clusters, rep(cl, length(genes_sorted)))
}

bin_scaled_ordered <- bin_scaled[gene_order, ]
cat(sprintf("  %d genes in %d clusters\n", length(gene_order), k))

# =============================================================================
# 7. Build annotation data frames
# =============================================================================

# Column annotation: pseudotime group (majority per bin)
col_anno <- data.frame(
  Group = factor(bin_group, levels = names(group_colors)),
  row.names = colnames(bin_scaled_ordered)
)

# Row annotation: gene cluster
gene_cluster_labels <- paste0("Module ", match(gene_clusters, cluster_order))
row_anno <- data.frame(
  Module = factor(gene_cluster_labels, levels = paste0("Module ", 1:k)),
  row.names = gene_order
)

# Annotation colors
module_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
names(module_palette) <- paste0("Module ", 1:k)

anno_colors <- list(
  Group = group_colors,
  Module = module_palette
)

# =============================================================================
# 8. Plot heatmap
# =============================================================================
cat("Generating heatmap ...\n")

heatmap_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                     "#F7F7F7",
                                     "#FDDBC7", "#F4A582", "#B2182B"))(100)

# Gaps between gene clusters
cluster_sizes <- table(factor(gene_cluster_labels, levels = paste0("Module ", 1:k)))
row_gaps <- cumsum(cluster_sizes)[-length(cluster_sizes)]

pdf_path <- file.path(out_dir, "fig2J_pseudotime_heatmap.pdf")
cairo_pdf(pdf_path, width = 10, height = 12)
pheatmap(bin_scaled_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = heatmap_colors,
         breaks = seq(-2.5, 2.5, length.out = 101),
         annotation_col = col_anno,
         annotation_row = row_anno,
         annotation_colors = anno_colors,
         gaps_row = row_gaps,
         border_color = NA,
         main = "Gene Expression Dynamics Along Pseudotime",
         fontsize = 10,
         annotation_names_col = TRUE,
         annotation_names_row = TRUE,
         annotation_legend = TRUE)
dev.off()
cat(sprintf("  Saved: %s\n", pdf_path))

# =============================================================================
# 9. Export gene lists per module
# =============================================================================
cat("Exporting gene module lists ...\n")
module_df <- data.frame(
  gene = gene_order,
  module = gene_cluster_labels,
  peak_bin = apply(bin_scaled_ordered, 1, which.max),
  peak_pseudotime = bin_pt[apply(bin_scaled_ordered, 1, which.max)],
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(out_dir, "fig2J_gene_modules.csv"), row.names = FALSE)

# Summary
cat("\nGene module summary:\n")
for (m in paste0("Module ", 1:k)) {
  genes_m <- module_df$gene[module_df$module == m]
  cat(sprintf("  %s: %d genes, median peak pseudotime = %.1f\n",
              m, length(genes_m), median(module_df$peak_pseudotime[module_df$module == m])))
  cat(sprintf("    Top genes: %s\n",
              paste(head(genes_m, 10), collapse = ", ")))
}

cat(sprintf("\nDone! Fig 2J saved to: %s\n", out_dir))
