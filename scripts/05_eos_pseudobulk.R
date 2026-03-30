#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## Pseudobulk analysis: all immune cells (Fig 1B) + eosinophils only
## Aggregates ALL genes per sample, then rlog → PCA/MDS (top 500 variable genes)
## SLURM: --mem=128G --cpus-per-task=4 --time=2:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out_dir <- file.path(BASE, "pseudobulk")
fig_dir <- file.path(BASE, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load annotated object ─────────────────────────────────────────────────────
cat("Loading all_cells_annotated.rds...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
cat(sprintf("  %d genes x %d cells\n", nrow(obj), ncol(obj)))

# ── Helper: Fig 1B aesthetics ─────────────────────────────────────────────────
# Color  = Condition (one distinct color per group, from cond_colors)
# Shape  = Genotype (WT=circle, KO=triangle) — fillable shapes
# Fill   = Infected=filled with condition color, Naive=open (white fill)
add_fig1b_aes <- function(df) {
  df$genotype  <- ifelse(grepl("^WT", df$condition), "WT", "KO")
  df$sex       <- ifelse(grepl("Male", df$condition), "Male", "Female")
  df$infection <- ifelse(grepl("Infected", df$condition), "Infected", "Naive")
  df
}

geno_shapes <- c("WT" = 21, "KO" = 24)

# Fill palette: infected conditions get their color, naive get white (open)
fill_palette <- cond_colors
fill_palette[!grepl("Infected", names(fill_palette))] <- "white"

fig1b_plot <- function(df, x_var, y_var, x_lab, y_lab, title_text) {
  # Ensure condition factor levels match cond_colors order
  present_conds <- intersect(names(cond_colors), unique(df$condition))
  df$condition <- factor(df$condition, levels = names(cond_colors))

  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                 color = condition, fill = condition, shape = genotype)) +
    geom_point(size = 4, stroke = 1.2) +
    scale_shape_manual(values = geno_shapes, name = "Genotype") +
    scale_color_manual(values = cond_colors, labels = cond_labels, name = "Group",
                       breaks = present_conds) +
    scale_fill_manual(values = fill_palette, labels = cond_labels, name = "Group",
                      breaks = present_conds) +
    guides(color = guide_legend(override.aes = list(
             fill = fill_palette[present_conds],
             shape = 21, size = 3)),
           fill = "none") +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_pub() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 14))
}

# ── Helper: pseudobulk pipeline ───────────────────────────────────────────────
run_pseudobulk <- function(sobj, label, prefix) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Pseudobulk: %s (%d cells)\n", label, ncol(sobj)))
  cat(sprintf("========================================\n\n"))

  # Aggregate ALL genes per sample
  cat("Aggregating raw counts per sample (all genes)...\n")
  counts_mat <- GetAssayData(sobj, layer = "counts")
  samples <- unique(sobj$orig.ident)

  pb_list <- lapply(samples, function(s) {
    cells <- which(sobj$orig.ident == s)
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples

  # Remove genes with zero counts across all samples
  pb_counts <- pb_counts[rowSums(pb_counts) > 0, ]
  cat(sprintf("  Pseudobulk matrix: %d genes x %d samples\n",
              nrow(pb_counts), ncol(pb_counts)))

  # Build sample metadata
  sample_meta <- sobj@meta.data %>%
    select(orig.ident, condition) %>%
    distinct()

  col_data <- data.frame(
    sample    = samples,
    condition = sample_meta$condition[match(samples, sample_meta$orig.ident)],
    row.names = samples
  )

  # DESeq2 rlog on ALL genes (recommended over VST for n < 30)
  cat("Running DESeq2 rlog (blind=TRUE) on all genes...\n")
  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                colData   = col_data,
                                design    = ~ 1)
  rlog_obj <- rlog(dds, blind = TRUE)
  rlog_mat <- assay(rlog_obj)
  cat(sprintf("  rlog matrix: %d genes x %d samples\n",
              nrow(rlog_mat), ncol(rlog_mat)))

  # PCA on top 500 most variable genes (standard DESeq2 approach)
  gene_vars <- apply(rlog_mat, 1, var)
  ntop <- min(500, sum(gene_vars > 0))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:ntop]
  cat(sprintf("  Using top %d variable genes for PCA (of %d total)\n",
              ntop, sum(gene_vars > 0)))
  pca_res <- prcomp(t(rlog_mat[top_genes, ]), scale. = FALSE)

  pca_df <- as.data.frame(pca_res$x[, 1:min(5, ncol(pca_res$x))])
  pca_df$sample    <- rownames(pca_df)
  pca_df$condition <- col_data[pca_df$sample, "condition"]
  pca_df <- add_fig1b_aes(pca_df)

  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  # MDS on top 500 variable genes (1 - Spearman correlation)
  cat("Computing MDS (1 - Spearman correlation) on top variable genes...\n")
  cor_mat  <- cor(rlog_mat[top_genes, ], method = "spearman")
  dist_mat <- as.dist(1 - cor_mat)
  mds_res  <- cmdscale(dist_mat, k = 2, eig = TRUE)

  mds_df <- data.frame(
    MDS1      = mds_res$points[, 1],
    MDS2      = mds_res$points[, 2],
    sample    = rownames(mds_res$points),
    condition = col_data[rownames(mds_res$points), "condition"]
  )
  mds_df <- add_fig1b_aes(mds_df)

  # Export PC loadings
  cat("Exporting PC loadings...\n")
  loadings <- pca_res$rotation
  for (pc_idx in 1:min(3, ncol(loadings))) {
    pc_load <- sort(loadings[, pc_idx], decreasing = TRUE)
    load_df <- data.frame(
      gene      = c(names(head(pc_load, 20)), names(tail(pc_load, 20))),
      loading   = c(head(pc_load, 20), tail(pc_load, 20)),
      direction = c(rep("positive", 20), rep("negative", 20))
    )
    write.csv(load_df, file.path(out_dir,
              sprintf("%s_pc%d_loadings.csv", prefix, pc_idx)),
              row.names = FALSE)
  }

  # --- Figures ---------------------------------------------------------------

  # MDS plot (Fig 1B style) — cairo_pdf for Unicode sex symbols
  cat("Plotting MDS...\n")
  cairo_pdf(file.path(fig_dir, sprintf("%s_mds.pdf", prefix)), width = 9, height = 7)
  print(fig1b_plot(mds_df, "MDS1", "MDS2",
                   "Dimension 1", "Dimension 2", label))
  dev.off()

  # PCA plot (Fig 1B style)
  cat("Plotting PCA...\n")
  cairo_pdf(file.path(fig_dir, sprintf("%s_pca.pdf", prefix)), width = 9, height = 7)
  print(fig1b_plot(pca_df, "PC1", "PC2",
                   sprintf("PC1 (%.1f%%)", pct_var[1]),
                   sprintf("PC2 (%.1f%%)", pct_var[2]),
                   label))
  dev.off()

  # Save pseudobulk counts
  cat("Saving pseudobulk counts...\n")
  write.csv(pb_counts, file.path(out_dir,
            sprintf("%s_counts.csv", prefix)))
  write.csv(pca_df, file.path(out_dir,
            sprintf("%s_pca_coordinates.csv", prefix)), row.names = FALSE)
  write.csv(mds_df, file.path(out_dir,
            sprintf("%s_mds_coordinates.csv", prefix)), row.names = FALSE)

  cat(sprintf("Done: %s\n", label))
}

# ── Run pseudobulk: All immune cells (Fig 1B) ────────────────────────────────
immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils", "Mast cells")
immune <- subset(obj, celltype %in% immune_types)
run_pseudobulk(immune, "Immune Cells", "05_immune")
rm(immune); gc()

# ── Run pseudobulk: Eosinophils only ─────────────────────────────────────────
eos <- subset(obj, celltype == "Eosinophils")
run_pseudobulk(eos, "Eosinophils", "05_eos")
rm(eos, obj); gc()

cat("\nDone! 05_eos_pseudobulk.R completed successfully.\n")
