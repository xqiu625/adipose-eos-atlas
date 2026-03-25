#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --job-name=qc_supplement
#SBATCH --output=qc_supplement_%j.out

# ==============================================================================
# 17_qc_supplement.R
# Supplemental QC figures showing data quality
# Addresses Reviewer 1 point 3a
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

SCRIPTS <- dirname(sys.frame(1)$ofile)
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures", "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Sample metadata ──────────────────────────────────────────────────────────
sample_info <- data.frame(
  sample_id = 1:16,
  condition = c("WT_Male_HFD_1",       "WT_Female_HFD_1",
                "KO_Male_HFD_1",       "KO_Female_HFD_1",
                "WT_Female_HFDInfected_1", "WT_Female_HFDInfected_2",
                "KO_Female_HFDInfected_1", "KO_Female_HFDInfected_2",
                "WT_Male_HFD_2",       "WT_Male_HFD_3",
                "WT_Female_HFD_2",     "WT_Female_HFD_3",
                "KO_Male_HFD_2",       "KO_Male_HFD_3",
                "KO_Female_HFD_2",     "KO_Female_HFD_3"),
  gem_reaction = c(rep("GEM1", 4), rep("GEM2", 4),
                   rep("GEM3", 4), rep("GEM4", 4)),
  stringsAsFactors = FALSE
)

# ── Load post-QC Seurat object ──────────────────────────────────────────────
cat("Loading seurat_qc_filtered.rds...\n")
sobj <- readRDS(file.path(BASE, "seurat_qc_filtered.rds"))
cat(sprintf("  %d genes x %d cells\n", nrow(sobj), ncol(sobj)))

# ── Load DecontX contamination scores ────────────────────────────────────────
cat("Loading DecontX scores...\n")
decontx_scores <- read.csv(file.path(BASE, "decontx", "decontx_all_scores.csv"),
                           stringsAsFactors = FALSE)
cat(sprintf("  %d entries\n", nrow(decontx_scores)))

# ── Prepare metadata ─────────────────────────────────────────────────────────
meta <- sobj@meta.data
meta$cell_barcode <- rownames(meta)

# Order samples consistently
sample_order <- sample_info$condition
meta$orig.ident <- factor(meta$orig.ident, levels = sample_order)

# ── Figure 1: Violin plot — nFeature_RNA per sample ─────────────────────────
cat("Plotting nFeature_RNA violin...\n")
p_nfeat <- ggplot(meta, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, show.legend = FALSE) +
  labs(x = "Sample", y = "nFeature_RNA",
       title = "Genes Detected per Cell (post-QC)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ── Figure 2: Violin plot — nCount_RNA per sample ───────────────────────────
cat("Plotting nCount_RNA violin...\n")
p_ncount <- ggplot(meta, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, show.legend = FALSE) +
  scale_y_log10() +
  labs(x = "Sample", y = "nCount_RNA (log10)",
       title = "UMI Counts per Cell (post-QC)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ── Figure 3: Violin plot — percent.mt per sample ───────────────────────────
cat("Plotting percent.mt violin...\n")
p_mt <- ggplot(meta, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  labs(x = "Sample", y = "Percent Mitochondrial",
       title = "Mitochondrial Content per Cell (post-QC)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ── Figure 4: Scatter — nCount vs nFeature colored by percent.mt ────────────
cat("Plotting nCount vs nFeature scatter...\n")
p_scatter <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_color_viridis_c(name = "% mt") +
  scale_x_log10() +
  facet_wrap(~ orig.ident, ncol = 4, scales = "free") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA",
       title = "UMI Count vs Gene Count (colored by % mito)") +
  theme_pub() +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.position = "right")

# ── Figure 5: DecontX contamination histogram per sample ─────────────────────
cat("Plotting DecontX contamination histograms...\n")
decontx_scores$condition <- factor(decontx_scores$condition,
                                   levels = sample_order)

p_decontx <- ggplot(decontx_scores, aes(x = contamination)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ condition, ncol = 4, scales = "free_y") +
  labs(x = "DecontX Contamination Score", y = "Count",
       title = "DecontX Contamination Distribution per Sample (threshold = 0.3)") +
  theme_pub() +
  theme(strip.text = element_text(size = 6),
        axis.text = element_text(size = 6))

# ── Figure 6: Cell count bar plot across filtering stages ────────────────────
cat("Computing cell counts at each filtering stage...\n")

# Raw cell counts from DecontX scores (pre-filtering)
raw_counts <- decontx_scores %>%
  group_by(condition) %>%
  summarise(n_raw = n(), .groups = "drop")

# After DecontX filtering (contamination <= 0.3)
after_decontx <- decontx_scores %>%
  filter(contamination <= 0.3) %>%
  group_by(condition) %>%
  summarise(n_after_decontx = n(), .groups = "drop")

# Post-QC counts (from the loaded object, which is post-QC + post-doublet)
post_final <- meta %>%
  group_by(orig.ident) %>%
  summarise(n_final = n(), .groups = "drop") %>%
  rename(condition = orig.ident)
post_final$condition <- as.character(post_final$condition)

# Merge all stages
count_stages <- raw_counts %>%
  left_join(after_decontx, by = "condition") %>%
  left_join(post_final, by = "condition") %>%
  replace_na(list(n_after_decontx = 0, n_final = 0))

# Pivot to long format for plotting
count_long <- count_stages %>%
  pivot_longer(cols = c(n_raw, n_after_decontx, n_final),
               names_to = "stage", values_to = "n_cells") %>%
  mutate(stage = factor(stage,
                        levels = c("n_raw", "n_after_decontx", "n_final"),
                        labels = c("Raw", "After DecontX", "After QC + Doublet")),
         condition = factor(condition, levels = sample_order))

p_barplot <- ggplot(count_long, aes(x = condition, y = n_cells, fill = stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Raw" = "#4393C3",
                                "After DecontX" = "#FDAE61",
                                "After QC + Doublet" = "#D73027")) +
  labs(x = "Sample", y = "Number of Cells", fill = "Stage",
       title = "Cell Counts Across Filtering Stages") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ── Combined multi-panel figure ──────────────────────────────────────────────
cat("Assembling combined supplement figure...\n")

pdf(file.path(fig_dir, "fig_qc_supplement.pdf"), width = 16, height = 24)
print(
  (p_nfeat / p_ncount / p_mt / p_barplot) +
    plot_annotation(
      title = "Supplementary Figure: Quality Control Metrics",
      tag_levels = "A"
    )
)
dev.off()

# Save the scatter and decontx panels separately (they need more space)
pdf(file.path(fig_dir, "fig_qc_scatter_per_sample.pdf"), width = 14, height = 14)
print(p_scatter)
dev.off()

pdf(file.path(fig_dir, "fig_qc_decontx_histograms.pdf"), width = 14, height = 10)
print(p_decontx)
dev.off()

# ── Figure 7: Summary CSV ────────────────────────────────────────────────────
cat("Computing per-sample summary table...\n")

# Per-sample medians from post-QC object
qc_summary <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    median_nFeature = median(nFeature_RNA),
    median_nCount   = median(nCount_RNA),
    median_pct_mt   = median(percent.mt),
    .groups = "drop"
  ) %>%
  rename(condition = orig.ident)
qc_summary$condition <- as.character(qc_summary$condition)

# Median contamination from DecontX scores (pre-filter, all cells)
contam_summary <- decontx_scores %>%
  group_by(condition) %>%
  summarise(median_contamination = median(contamination), .groups = "drop")

# Merge everything
summary_table <- count_stages %>%
  left_join(qc_summary, by = "condition") %>%
  left_join(contam_summary, by = "condition") %>%
  rename(
    n_cells_raw          = n_raw,
    n_cells_after_decontx = n_after_decontx,
    n_cells_after_qc_doublet = n_final
  )

# Add GEM reaction info
summary_table <- summary_table %>%
  left_join(sample_info[, c("condition", "gem_reaction")], by = "condition")

# Reorder columns
summary_table <- summary_table %>%
  select(condition, gem_reaction, n_cells_raw, n_cells_after_decontx,
         n_cells_after_qc_doublet, median_nFeature, median_nCount,
         median_pct_mt, median_contamination)

csv_out <- file.path(fig_dir, "qc_summary_table.csv")
write.csv(summary_table, csv_out, row.names = FALSE)
cat(sprintf("Saved summary table: %s\n", csv_out))

# Print summary for log
cat("\n========================================\n")
cat("QC Summary Table\n")
cat("========================================\n")
print(as.data.frame(summary_table))

cat("\n========================================\n")
cat("Overall totals:\n")
cat(sprintf("  Raw cells:                    %d\n", sum(summary_table$n_cells_raw)))
cat(sprintf("  After DecontX:                %d\n", sum(summary_table$n_cells_after_decontx)))
cat(sprintf("  After QC + doublet removal:   %d\n", sum(summary_table$n_cells_after_qc_doublet)))
cat("========================================\n")

cat("\nDone! 17_qc_supplement.R completed successfully.\n")
