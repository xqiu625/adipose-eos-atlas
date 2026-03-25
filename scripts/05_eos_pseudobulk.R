#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=eos_pseudobulk
#SBATCH --output=eos_pseudobulk_%j.out

# ==============================================================================
# 05_eos_pseudobulk.R
# Pseudobulk analysis of eosinophils: PCA, MDS, and batch assessment
# ==============================================================================

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

# --- Load eosinophil object ---------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))

# --- Derive GEM reaction mapping from sample names ----------------------------
cat("Deriving GEM reaction assignments...\n")
sample_meta <- eos@meta.data %>%
  select(orig.ident, condition) %>%
  distinct()

sample_meta <- sample_meta %>%
  mutate(
    replicate_num = as.integer(sub(".*_(\\d+)$", "\\1", orig.ident)),
    is_infected   = grepl("Infected", condition),
    genotype      = sub("^(WT|KO)_.*", "\\1", condition),
    gem_reaction  = case_when(
      is_infected                              ~ "GEM2",
      !is_infected & replicate_num == 1        ~ "GEM1",
      !is_infected & genotype == "WT" &
        replicate_num %in% 2:3                 ~ "GEM3",
      !is_infected & genotype == "KO" &
        replicate_num %in% 2:3                 ~ "GEM4",
      TRUE                                     ~ "Unknown"
    )
  )

# Add gem_reaction back to Seurat object
gem_map <- setNames(sample_meta$gem_reaction, sample_meta$orig.ident)
eos$gem_reaction <- gem_map[eos$orig.ident]

cat("  GEM assignments:\n")
print(table(sample_meta$orig.ident, sample_meta$gem_reaction))

# --- Aggregate raw counts per sample (pseudobulk) ----------------------------
cat("Aggregating raw counts per sample...\n")
counts_mat <- GetAssayData(eos, slot = "counts")
samples <- unique(eos$orig.ident)

pb_list <- lapply(samples, function(s) {
  cells <- which(eos$orig.ident == s)
  Matrix::rowSums(counts_mat[, cells, drop = FALSE])
})
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- samples

# Remove genes with zero counts across all samples
pb_counts <- pb_counts[rowSums(pb_counts) > 0, ]
cat(sprintf("  Pseudobulk matrix: %d genes x %d samples\n",
            nrow(pb_counts), ncol(pb_counts)))

# --- DESeq2 VST normalization -------------------------------------------------
cat("Running DESeq2 VST (blind=TRUE)...\n")
col_data <- data.frame(
  sample    = samples,
  condition = sample_meta$condition[match(samples, sample_meta$orig.ident)],
  gem       = sample_meta$gem_reaction[match(samples, sample_meta$orig.ident)],
  row.names = samples
)

dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                              colData   = col_data,
                              design    = ~ 1)
vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)
cat("  VST normalization complete.\n")

# --- PCA on top 500 most variable genes ---------------------------------------
cat("Running PCA on top 500 variable genes...\n")
gene_vars <- apply(vst_mat, 1, var)
top500 <- names(sort(gene_vars, decreasing = TRUE))[1:min(500, length(gene_vars))]
pca_res <- prcomp(t(vst_mat[top500, ]), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:min(5, ncol(pca_res$x))])
pca_df$sample    <- rownames(pca_df)
pca_df$condition <- col_data[pca_df$sample, "condition"]
pca_df$gem       <- col_data[pca_df$sample, "gem"]
pca_df$genotype  <- ifelse(grepl("^WT", pca_df$condition), "WT", "KO")

pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

# Export PC loadings (top genes driving each PC)
cat("Exporting PC loadings...\n")
loadings <- pca_res$rotation
for (pc_idx in 1:min(3, ncol(loadings))) {
  pc_name <- paste0("PC", pc_idx)
  pc_load <- sort(loadings[, pc_idx], decreasing = TRUE)
  top_pos <- head(pc_load, 20)
  top_neg <- tail(pc_load, 20)
  load_df <- data.frame(
    gene    = c(names(top_pos), names(top_neg)),
    loading = c(top_pos, top_neg),
    direction = c(rep("positive", length(top_pos)), rep("negative", length(top_neg)))
  )
  write.csv(load_df, file.path(out_dir, sprintf("pc%d_loadings.csv", pc_idx)),
            row.names = FALSE)
}

# --- MDS with 1 - Spearman correlation distance ------------------------------
cat("Computing MDS (1 - Spearman correlation)...\n")
cor_mat  <- cor(vst_mat[top500, ], method = "spearman")
dist_mat <- as.dist(1 - cor_mat)
mds_res  <- cmdscale(dist_mat, k = 2, eig = TRUE)

mds_df <- data.frame(
  MDS1      = mds_res$points[, 1],
  MDS2      = mds_res$points[, 2],
  sample    = rownames(mds_res$points),
  condition = col_data[rownames(mds_res$points), "condition"]
)

# --- Figures ------------------------------------------------------------------

# PCA plot (color by condition, shape by genotype, label samples)
cat("Plotting PCA...\n")
pdf(file.path(fig_dir, "05_pseudobulk_pca.pdf"), width = 9, height = 7)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = condition,
                     label = sample)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(size = 4) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = cond_colors) +
    scale_shape_manual(values = cond_shapes) +
    labs(x = sprintf("PC1 (%.1f%%)", pct_var[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var[2]),
         title = "Eosinophil Pseudobulk PCA") +
    theme_pub() +
    theme(legend.position = "right")
)
dev.off()

# MDS plot
cat("Plotting MDS...\n")
pdf(file.path(fig_dir, "05_pseudobulk_mds.pdf"), width = 9, height = 7)
print(
  ggplot(mds_df, aes(x = MDS1, y = MDS2, color = condition, label = sample)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(size = 4) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = cond_colors) +
    labs(x = "MDS1", y = "MDS2",
         title = "Eosinophil Pseudobulk MDS (1 - Spearman)") +
    theme_pub()
)
dev.off()

# PCA colored by GEM reaction (batch check)
cat("Plotting PCA by GEM reaction...\n")
pdf(file.path(fig_dir, "05_pseudobulk_pca_gem.pdf"), width = 9, height = 7)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = gem, label = sample)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(size = 4) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    labs(x = sprintf("PC1 (%.1f%%)", pct_var[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var[2]),
         title = "Eosinophil Pseudobulk PCA — Batch (GEM Reaction)") +
    theme_pub()
)
dev.off()

# --- Save outputs -------------------------------------------------------------
cat("Saving pseudobulk counts and coordinates...\n")
write.csv(pb_counts, file.path(out_dir, "pseudobulk_counts.csv"))
write.csv(pca_df,    file.path(out_dir, "pseudobulk_coordinates.csv"), row.names = FALSE)

cat("Done! 05_eos_pseudobulk.R completed successfully.\n")
