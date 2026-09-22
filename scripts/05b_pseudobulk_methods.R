#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

## Try multiple pseudobulk normalization/dim-reduction combinations
## to find best separation of biological conditions
## SLURM: --mem=128G --cpus-per-task=4 --time=2:00:00 -p batch

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading all_cells_annotated.rds...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

# ── Fig 1B aesthetics ────────────────────────────────────────────────────────
add_aes <- function(df) {
  df$genotype  <- ifelse(grepl("^WT", df$condition), "WT", "KO")
  df$sex       <- ifelse(grepl("Male", df$condition), "Male", "Female")
  df$infection <- ifelse(grepl("Infected", df$condition), "Infected", "Naive")
  df
}

sex_colors  <- c("Male" = "#2166AC", "Female" = "#B2182B")
geno_shapes <- c("WT" = 21, "KO" = 24)

make_plot <- function(df, x_var, y_var, x_lab, y_lab, title_text) {
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                 color = sex, fill = sex, shape = genotype)) +
    geom_point(size = 3.5, stroke = 1.0, aes(alpha = infection)) +
    scale_shape_manual(values = geno_shapes) +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    scale_alpha_manual(values = c("Naive" = 0.4, "Infected" = 1)) +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_pub(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold"))
}

# ── Helper: aggregate pseudobulk counts ──────────────────────────────────────
aggregate_pb <- function(sobj) {
  counts_mat <- GetAssayData(sobj, layer = "counts")
  samples <- unique(sobj$orig.ident)
  pb_list <- lapply(samples, function(s) {
    cells <- which(sobj$orig.ident == s)
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples
  pb_counts <- pb_counts[rowSums(pb_counts) > 0, ]

  sample_meta <- sobj@meta.data %>%
    select(orig.ident, condition) %>%
    distinct()
  col_data <- data.frame(
    sample = samples,
    condition = sample_meta$condition[match(samples, sample_meta$orig.ident)],
    row.names = samples
  )
  list(counts = pb_counts, col_data = col_data)
}

# ── Helper: get top variable genes ───────────────────────────────────────────
top_var_genes <- function(mat, ntop) {
  gene_vars <- apply(mat, 1, var)
  ntop <- min(ntop, sum(gene_vars > 0))
  names(sort(gene_vars, decreasing = TRUE))[1:ntop]
}

# ── Run all methods for a given cell subset ──────────────────────────────────
run_all_methods <- function(sobj, label, prefix) {
  cat(sprintf("\n====== %s (%d cells) ======\n", label, ncol(sobj)))

  pb <- aggregate_pb(sobj)
  pb_counts <- pb$counts
  col_data  <- pb$col_data
  cat(sprintf("  Pseudobulk: %d genes x %d samples\n", nrow(pb_counts), ncol(pb_counts)))

  plots <- list()

  # --- Normalization methods ---

  # 1. DESeq2 VST
  cat("  VST...\n")
  dds <- DESeqDataSetFromMatrix(round(pb_counts), col_data, ~ 1)
  vst_mat <- assay(vst(dds, blind = TRUE))

  # 2. DESeq2 rlog
  cat("  rlog...\n")
  rlog_mat <- assay(rlog(dds, blind = TRUE))

  # 3. edgeR logCPM
  cat("  logCPM...\n")
  dge <- DGEList(counts = round(pb_counts))
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm_mat <- cpm(dge, log = TRUE, prior.count = 1)

  # 4. Simple log2(CPM+1)
  cat("  log2(CPM+1)...\n")
  cpm_raw <- t(t(pb_counts) / colSums(pb_counts) * 1e6)
  log2cpm_mat <- log2(cpm_raw + 1)

  norm_methods <- list(
    "VST"          = vst_mat,
    "rlog"         = rlog_mat,
    "edgeR_logCPM" = logcpm_mat,
    "log2CPM1"     = log2cpm_mat
  )

  # --- Gene selection x Normalization x Dim reduction ---
  ntops <- c(200, 500, 1000, 2000)

  for (norm_name in names(norm_methods)) {
    mat <- norm_methods[[norm_name]]

    for (ntop in ntops) {
      top_g <- top_var_genes(mat, ntop)
      mat_sub <- mat[top_g, ]
      tag <- sprintf("%s_top%d", norm_name, ntop)

      # PCA (no scaling — normalization already handles it)
      pca <- prcomp(t(mat_sub), scale. = FALSE)
      pca_df <- as.data.frame(pca$x[, 1:min(3, ncol(pca$x))])
      pca_df$sample <- rownames(pca_df)
      pca_df$condition <- col_data[pca_df$sample, "condition"]
      pca_df <- add_aes(pca_df)
      pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

      plots[[paste0("PCA_", tag)]] <- make_plot(
        pca_df, "PC1", "PC2",
        sprintf("PC1 (%.0f%%)", pct[1]),
        sprintf("PC2 (%.0f%%)", pct[2]),
        sprintf("PCA %s", tag))

      # PCA with scaling
      pca_s <- prcomp(t(mat_sub), scale. = TRUE)
      pca_s_df <- as.data.frame(pca_s$x[, 1:min(3, ncol(pca_s$x))])
      pca_s_df$sample <- rownames(pca_s_df)
      pca_s_df$condition <- col_data[pca_s_df$sample, "condition"]
      pca_s_df <- add_aes(pca_s_df)
      pct_s <- round(100 * pca_s$sdev^2 / sum(pca_s$sdev^2), 1)

      plots[[paste0("PCAscaled_", tag)]] <- make_plot(
        pca_s_df, "PC1", "PC2",
        sprintf("PC1 (%.0f%%)", pct_s[1]),
        sprintf("PC2 (%.0f%%)", pct_s[2]),
        sprintf("PCA(scaled) %s", tag))

      # MDS Spearman
      cor_sp <- cor(mat_sub, method = "spearman")
      mds_sp <- cmdscale(as.dist(1 - cor_sp), k = 2)
      mds_sp_df <- data.frame(MDS1 = mds_sp[, 1], MDS2 = mds_sp[, 2],
                               sample = rownames(mds_sp))
      mds_sp_df$condition <- col_data[mds_sp_df$sample, "condition"]
      mds_sp_df <- add_aes(mds_sp_df)

      plots[[paste0("MDS_Spearman_", tag)]] <- make_plot(
        mds_sp_df, "MDS1", "MDS2", "MDS1", "MDS2",
        sprintf("MDS(Spearman) %s", tag))

      # MDS Pearson
      cor_pe <- cor(mat_sub, method = "pearson")
      mds_pe <- cmdscale(as.dist(1 - cor_pe), k = 2)
      mds_pe_df <- data.frame(MDS1 = mds_pe[, 1], MDS2 = mds_pe[, 2],
                               sample = rownames(mds_pe))
      mds_pe_df$condition <- col_data[mds_pe_df$sample, "condition"]
      mds_pe_df <- add_aes(mds_pe_df)

      plots[[paste0("MDS_Pearson_", tag)]] <- make_plot(
        mds_pe_df, "MDS1", "MDS2", "MDS1", "MDS2",
        sprintf("MDS(Pearson) %s", tag))
    }
  }

  cat(sprintf("  Generated %d plots\n", length(plots)))

  # --- Arrange into pages: group by dim-reduction method ---
  # 4 norm x 4 ntop = 16 panels per dim-reduction method
  dr_methods <- c("PCA", "PCAscaled", "MDS_Spearman", "MDS_Pearson")

  pdf(file.path(fig_dir, sprintf("%s_method_comparison.pdf", prefix)),
      width = 20, height = 20)
  for (dr in dr_methods) {
    dr_plots <- plots[grep(paste0("^", dr, "_"), names(plots))]
    if (length(dr_plots) > 0) {
      # Order: by norm method then ntop
      p_combined <- wrap_plots(dr_plots, ncol = 4) +
        plot_annotation(
          title = sprintf("%s: %s — Method Comparison", label, dr),
          theme = theme(plot.title = element_text(face = "bold", size = 14))
        )
      print(p_combined)
    }
  }
  dev.off()
  cat(sprintf("  Saved: %s_method_comparison.pdf\n", prefix))
}

# ── Run for immune cells and eosinophils ─────────────────────────────────────
immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils", "Mast cells")
immune <- subset(obj, celltype %in% immune_types)
run_all_methods(immune, "Immune Cells", "05b_immune")
rm(immune); gc()

eos <- subset(obj, celltype == "Eosinophils")
run_all_methods(eos, "Eosinophils", "05b_eos")
rm(eos, obj); gc()

cat("\nDone! 05b_pseudobulk_methods.R\n")
