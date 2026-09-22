#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=ambient_evidence
#SBATCH --output=ambient_evidence_%j.out

# ==============================================================================
# 21_ambient_rna_evidence.R
# Evidence that DecontX effectively removed ambient RNA and that residual
# Cd14/macrophage marker expression is not contamination
# (ambient-RNA evidence and QC)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Matrix)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures", "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("21_ambient_rna_evidence.R\n")
cat("DecontX ambient RNA evidence\n")
cat("========================================\n\n")

# ==============================================================================
# Load data
# ==============================================================================
cat("Loading data...\n")

# Contamination scores
scores_path <- file.path(BASE, "decontx", "decontx_all_scores.csv")
cat(sprintf("  Contamination scores: %s\n", scores_path))
contam_scores <- read.csv(scores_path, stringsAsFactors = FALSE)
cat(sprintf("  %d cells with contamination scores\n", nrow(contam_scores)))

# All-cells annotated object (has celltype labels)
cat("  Loading all_cells_annotated.rds...\n")
all_cells <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
cat(sprintf("  %d cells, %d cell types\n", ncol(all_cells),
            length(unique(all_cells$broad_celltype))))

# Eos subclustered object
cat("  Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
# Remove C6 if still present (too few cells, n=14)
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  eos$eos_cluster <- droplevels(eos$eos_cluster)
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}
cat(sprintf("  %d eosinophils, %d clusters\n", ncol(eos),
            length(unique(eos$eos_cluster))))

# ==============================================================================
# Analysis 1: DecontX contamination by cell type
# ==============================================================================
cat("\n--- Analysis 1: Contamination by cell type ---\n")

# Match contamination scores to cell type annotations
# contam_scores has barcode and sample_id columns
# all_cells has barcodes as column names

# The barcodes in contam_scores are raw CellRanger barcodes,
# while all_cells colnames are condition_barcode format.
# We need to match them via the cell metadata.
# The all_cells object should have orig.ident and original barcodes embedded.

# Strategy: use the all_cells metadata which has contamination info if available,
# or match by reconstructing the barcode mapping.
# Since contam_scores has barcode + sample_id, and all_cells colnames are
# condition_barcode, we can try to match.

# Build a lookup from condition+barcode -> celltype
cell_meta <- all_cells@meta.data
cell_meta$cell_id <- rownames(cell_meta)

# Extract the raw barcode from the cell_id (last component after condition prefix)
cell_meta$raw_barcode <- sub("^.*_([ACGT]+-\\d+)$", "\\1", cell_meta$cell_id)

# Build sample_id -> condition mapping from the decontx script's sample_info
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
  stringsAsFactors = FALSE
)

# Add condition to contam_scores and create matching key
contam_scores <- contam_scores %>%
  left_join(sample_info, by = "sample_id") %>%
  rename(orig_condition = condition.x, sample_condition = condition.y)

# Create cell_id to match all_cells colnames: condition_barcode
contam_scores$cell_id <- paste0(contam_scores$sample_condition, "_",
                                contam_scores$barcode)

# Merge with cell type info
contam_with_type <- contam_scores %>%
  inner_join(
    cell_meta %>% select(cell_id, broad_celltype),
    by = "cell_id"
  )

cat(sprintf("  Matched %d / %d cells to cell type annotations\n",
            nrow(contam_with_type), nrow(contam_scores)))

if (nrow(contam_with_type) == 0) {
  cat("  WARNING: No matches found. Trying alternative matching strategy...\n")
  # Alternative: match by raw barcode + orig.ident
  contam_with_type <- contam_scores %>%
    inner_join(
      cell_meta %>% select(cell_id, raw_barcode, broad_celltype,
                           orig.ident = orig.ident),
      by = c("barcode" = "raw_barcode", "sample_condition" = "orig.ident")
    )
  cat(sprintf("  Alternative matching: %d cells\n", nrow(contam_with_type)))
}

# Figure 1: Violin plot of contamination per cell type
cat("  Plotting contamination violin by cell type...\n")

if (nrow(contam_with_type) > 0) {
  # Order cell types by median contamination
  type_order <- contam_with_type %>%
    group_by(broad_celltype) %>%
    summarise(med = median(contamination), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(broad_celltype)

  contam_with_type$broad_celltype <- factor(contam_with_type$broad_celltype,
                                            levels = type_order)

  p1 <- ggplot(contam_with_type, aes(x = broad_celltype, y = contamination,
                                      fill = broad_celltype)) +
    geom_violin(scale = "width", alpha = 0.7, show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.5) +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red", linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.32, label = "Threshold = 0.3",
             hjust = 0, color = "red", size = 3) +
    labs(x = "Cell Type", y = "DecontX Contamination Score",
         title = "Ambient RNA Contamination by Cell Type") +
    theme_pub(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(file.path(fig_dir, "fig_ambient_evidence_contamination_by_type.pdf"),
      width = 8, height = 5)
  print(p1)
  dev.off()
  cat("  Saved: fig_ambient_evidence_contamination_by_type.pdf\n")
} else {
  cat("  Skipping figure (no matched cells)\n")
  p1 <- NULL
}

# ==============================================================================
# Analysis 2: Macrophage marker expression before vs after DecontX
# ==============================================================================
cat("\n--- Analysis 2: Macrophage marker expression before vs after DecontX ---\n")

mac_markers <- c("Cd14", "Adgre1", "Csf1r", "Cd68", "Lyz2")

# Load per-sample DecontX RDS files to get original and corrected counts
# for eosinophil cells only
cat("  Loading per-sample DecontX results for eos cells...\n")

# Get eos cell barcodes and their sample assignments
eos_meta <- eos@meta.data
eos_meta$cell_id <- rownames(eos_meta)
eos_meta$raw_barcode <- sub("^.*_([ACGT]+-\\d+)$", "\\1", eos_meta$cell_id)

# For each sample, load decontx RDS and extract original + corrected counts
# for eosinophil cells
before_list <- list()
after_list  <- list()

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  cond <- sample_info$condition[i]

  rds_path <- file.path(BASE, "decontx", sprintf("decontx_sample_%d.rds", sid))
  if (!file.exists(rds_path)) {
    cat(sprintf("    Sample %d: RDS not found, skipping\n", sid))
    next
  }

  cat(sprintf("    Sample %d (%s)...\n", sid, cond))
  decon <- readRDS(rds_path)

  # The decontx result has: counts (corrected), and we need original counts
  # Original filtered counts were used as input to decontX
  # corrected counts are in decon$counts

  corrected_counts <- decon$counts

  # Identify eos cells from this sample
  eos_in_sample <- eos_meta %>% filter(orig.ident == cond)
  if (nrow(eos_in_sample) == 0) next

  # Match barcodes
  common_bc <- intersect(eos_in_sample$raw_barcode, colnames(corrected_counts))
  if (length(common_bc) == 0) next

  # Get corrected counts for mac markers
  available_markers <- intersect(mac_markers, rownames(corrected_counts))
  if (length(available_markers) == 0) next

  corrected_sub <- as.matrix(corrected_counts[available_markers, common_bc, drop = FALSE])

  # Map to cell_ids for consistent naming
  bc_to_cellid <- setNames(eos_in_sample$cell_id, eos_in_sample$raw_barcode)
  colnames(corrected_sub) <- bc_to_cellid[common_bc]

  after_list[[as.character(sid)]] <- corrected_sub
}

# For "before" DecontX: use the RNA counts in the current eos object
# (which ARE post-DecontX), we need original.
# Alternative approach: re-load raw filtered counts per sample.
# Actually the Seurat pipeline uses DecontX-corrected counts already.
# So we need to load the original CellRanger filtered counts.

# Simpler approach: load per-sample decontx RDS which stores corrected counts,
# and compare with what is in the current Seurat object's raw RNA counts.
# But the Seurat object data slot has normalized (log-transformed) values.

# Best approach: load original CellRanger filtered counts for comparison.
# The paths are defined in the decontx script.
data_root <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data"
gem_map_list <- list(
  "RR0104" = 1:4,
  "RR0508" = 5:8,
  "RR0912" = 9:12,
  "RR1316" = 13:16
)

get_sample_path <- function(sample_id) {
  gem <- NULL
  for (g in names(gem_map_list)) {
    if (sample_id %in% gem_map_list[[g]]) {
      gem <- g
      break
    }
  }
  file.path(data_root, gem, "outs", "per_sample_outs",
            as.character(sample_id), "count", "sample_filtered_feature_bc_matrix")
}

cat("  Loading original (pre-DecontX) counts for eos cells...\n")

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  cond <- sample_info$condition[i]

  eos_in_sample <- eos_meta %>% filter(orig.ident == cond)
  if (nrow(eos_in_sample) == 0) next

  orig_path <- get_sample_path(sid)
  if (!dir.exists(orig_path)) {
    cat(sprintf("    Sample %d: original data not found at %s, skipping\n",
                sid, orig_path))
    next
  }

  cat(sprintf("    Sample %d (%s): loading original counts...\n", sid, cond))
  orig_counts <- tryCatch({
    Read10X(data.dir = orig_path)
  }, error = function(e) {
    cat(sprintf("    Error loading: %s\n", e$message))
    NULL
  })
  if (is.null(orig_counts)) next

  common_bc <- intersect(eos_in_sample$raw_barcode, colnames(orig_counts))
  if (length(common_bc) == 0) next

  available_markers <- intersect(mac_markers, rownames(orig_counts))
  if (length(available_markers) == 0) next

  orig_sub <- as.matrix(orig_counts[available_markers, common_bc, drop = FALSE])

  bc_to_cellid <- setNames(eos_in_sample$cell_id, eos_in_sample$raw_barcode)
  colnames(orig_sub) <- bc_to_cellid[common_bc]

  before_list[[as.character(sid)]] <- orig_sub
}

# Combine before and after matrices
if (length(before_list) > 0 && length(after_list) > 0) {
  before_mat <- do.call(cbind, before_list)
  after_mat  <- do.call(cbind, after_list)

  # Align cells present in both
  common_cells <- intersect(colnames(before_mat), colnames(after_mat))
  cat(sprintf("  Eos cells with both original and corrected counts: %d\n",
              length(common_cells)))

  if (length(common_cells) > 0) {
    before_mat <- before_mat[, common_cells, drop = FALSE]
    after_mat  <- after_mat[, common_cells, drop = FALSE]

    # Build long-format data for plotting
    available_markers <- rownames(before_mat)
    plot_data <- data.frame()

    for (gene in available_markers) {
      before_vals <- before_mat[gene, ]
      after_vals  <- after_mat[gene, ]

      df_gene <- data.frame(
        cell = rep(common_cells, 2),
        gene = gene,
        expression = c(before_vals, after_vals),
        stage = rep(c("Before DecontX", "After DecontX"), each = length(common_cells)),
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, df_gene)
    }

    plot_data$stage <- factor(plot_data$stage,
                              levels = c("Before DecontX", "After DecontX"))

    # Figure 2: Paired violin/boxplot for macrophage markers
    cat("  Plotting macrophage marker comparison...\n")

    p2 <- ggplot(plot_data, aes(x = stage, y = expression, fill = stage)) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.5) +
      facet_wrap(~ gene, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Before DecontX" = "#E84A5F",
                                   "After DecontX" = "#5BA4A4")) +
      labs(x = "", y = "Raw Counts",
           title = "Macrophage Marker Expression in Eosinophils",
           subtitle = "Before vs After DecontX Correction") +
      theme_pub(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none"
      )

    pdf(file.path(fig_dir, "fig_ambient_evidence_mac_markers.pdf"),
        width = 12, height = 5)
    print(p2)
    dev.off()
    cat("  Saved: fig_ambient_evidence_mac_markers.pdf\n")

    # Compute summary statistics for the table
    marker_summary <- plot_data %>%
      group_by(gene, stage) %>%
      summarise(
        mean_expr = mean(expression),
        median_expr = median(expression),
        pct_nonzero = 100 * mean(expression > 0),
        .groups = "drop"
      ) %>%
      arrange(gene, stage)

    cat("  Macrophage marker summary (before vs after DecontX):\n")
    print(as.data.frame(marker_summary))

    # --- Analysis 3: Top genes most reduced by DecontX in eosinophils -----------
    cat("\n--- Analysis 3: Ambient RNA gene profile ---\n")
    cat("  Identifying genes most reduced by DecontX...\n")

    # Use full gene sets from before/after matrices
    # Reload with all genes (not just mac markers) for a representative subset
    # For efficiency, compute reduction for available markers + top variable genes

    # Use the before and after counts for ALL genes from one representative sample
    # to identify the ambient RNA profile. We'll use the full-gene before/after
    # comparison from the per-sample DecontX RDS files.

    # Pick the sample with the most eos cells
    eos_per_sample <- eos_meta %>%
      group_by(orig.ident) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))

    top_sample_cond <- eos_per_sample$orig.ident[1]
    top_sample_id <- sample_info$sample_id[sample_info$condition == top_sample_cond]

    cat(sprintf("  Using sample %d (%s, %d eos cells) for gene profile\n",
                top_sample_id, top_sample_cond, eos_per_sample$n[1]))

    # Load original counts
    orig_path <- get_sample_path(top_sample_id)
    orig_full <- tryCatch(Read10X(data.dir = orig_path), error = function(e) NULL)

    # Load corrected counts
    rds_path <- file.path(BASE, "decontx", sprintf("decontx_sample_%d.rds", top_sample_id))
    decon_full <- tryCatch(readRDS(rds_path), error = function(e) NULL)

    if (!is.null(orig_full) && !is.null(decon_full)) {
      corrected_full <- decon_full$counts

      # Get eos barcodes for this sample
      eos_bc <- eos_meta %>%
        filter(orig.ident == top_sample_cond) %>%
        pull(raw_barcode)
      common_bc <- intersect(eos_bc, intersect(colnames(orig_full),
                                                colnames(corrected_full)))

      if (length(common_bc) > 0) {
        common_genes <- intersect(rownames(orig_full), rownames(corrected_full))
        orig_sub <- orig_full[common_genes, common_bc]
        corr_sub <- corrected_full[common_genes, common_bc]

        # Mean expression per gene: before vs after
        mean_before <- rowMeans(as.matrix(orig_sub))
        mean_after  <- rowMeans(as.matrix(corr_sub))

        # Absolute reduction
        reduction <- mean_before - mean_after

        # Top 20 most reduced genes
        top_reduced <- sort(reduction, decreasing = TRUE)[1:20]
        top_reduced_df <- data.frame(
          gene = names(top_reduced),
          mean_before = mean_before[names(top_reduced)],
          mean_after = mean_after[names(top_reduced)],
          reduction = as.numeric(top_reduced),
          stringsAsFactors = FALSE
        )
        top_reduced_df$gene <- factor(top_reduced_df$gene,
                                      levels = rev(top_reduced_df$gene))

        cat("  Top 20 genes most reduced by DecontX in eosinophils:\n")
        print(top_reduced_df)

        # Figure 3: Bar plot of top reduced genes
        p3 <- ggplot(top_reduced_df, aes(x = gene, y = reduction)) +
          geom_col(fill = "#E84A5F", alpha = 0.8) +
          geom_text(aes(label = sprintf("%.3f", reduction)),
                    hjust = -0.1, size = 3) +
          coord_flip() +
          labs(x = "Gene", y = "Mean Expression Reduction (Before - After DecontX)",
               title = "Top 20 Genes Most Reduced by DecontX in Eosinophils",
               subtitle = sprintf("Sample: %s (%d cells)", top_sample_cond,
                                  length(common_bc))) +
          theme_pub(base_size = 11) +
          theme(axis.text.y = element_text(face = "italic"))

        pdf(file.path(fig_dir, "fig_ambient_evidence_top_reduced_genes.pdf"),
            width = 8, height = 6)
        print(p3)
        dev.off()
        cat("  Saved: fig_ambient_evidence_top_reduced_genes.pdf\n")
      } else {
        cat("  No common barcodes for full gene analysis\n")
        p3 <- NULL
        top_reduced_df <- NULL
      }
    } else {
      cat("  Could not load full gene data for ambient profile\n")
      p3 <- NULL
      top_reduced_df <- NULL
    }

  } else {
    cat("  No overlapping cells between before/after. Skipping comparison.\n")
    p2 <- NULL
    p3 <- NULL
  }
} else {
  cat("  Insufficient data for before/after comparison.\n")
  cat("  Before samples loaded: %d, After samples loaded: %d\n",
      length(before_list), length(after_list))
  p2 <- NULL
  p3 <- NULL
}

# ==============================================================================
# Analysis 4: Per-sample contamination summary
# ==============================================================================
cat("\n--- Analysis 4: Per-sample contamination summary ---\n")

sample_summary <- contam_scores %>%
  group_by(sample_id, sample_condition) %>%
  summarise(
    n_cells = n(),
    mean_contam = mean(contamination),
    median_contam = median(contamination),
    pct_above_threshold = 100 * mean(contamination > 0.3),
    .groups = "drop"
  ) %>%
  arrange(sample_id)

cat("  Per-sample contamination summary:\n")
print(as.data.frame(sample_summary))

# Figure 4a: Bar plot of mean contamination per sample
p4a <- ggplot(sample_summary, aes(x = reorder(sample_condition, sample_id),
                                   y = mean_contam)) +
  geom_col(fill = "#6CA6CD", alpha = 0.8) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  labs(x = "Sample", y = "Mean Contamination Score",
       title = "Mean DecontX Contamination per Sample") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Figure 4b: Bar plot of % cells above threshold per sample
p4b <- ggplot(sample_summary, aes(x = reorder(sample_condition, sample_id),
                                   y = pct_above_threshold)) +
  geom_col(fill = "#E84A5F", alpha = 0.8) +
  labs(x = "Sample", y = "% Cells Above Threshold (0.3)",
       title = "Cells Removed per Sample (Contamination > 0.3)") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Figure 4c: Scatter of contamination vs nFeature_RNA
# Merge contamination with cell QC metrics from all_cells object
if (nrow(contam_with_type) > 0 && "nFeature_RNA" %in% colnames(cell_meta)) {
  scatter_data <- contam_with_type %>%
    inner_join(
      cell_meta %>% select(cell_id, nFeature_RNA),
      by = "cell_id"
    )

  p4c <- ggplot(scatter_data, aes(x = nFeature_RNA, y = contamination)) +
    geom_point(alpha = 0.05, size = 0.3, color = "#2F3C7E") +
    geom_smooth(method = "lm", color = "#E84A5F", se = TRUE) +
    geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
    labs(x = "nFeature_RNA", y = "Contamination Score",
         title = "Contamination vs Gene Complexity") +
    theme_pub(base_size = 11)
} else {
  p4c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "nFeature_RNA data not available") +
    theme_void()
}

# Combined per-sample summary figure
p4_combined <- (p4a | p4b) / p4c +
  plot_annotation(
    title = "Per-Sample Contamination Summary",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf(file.path(fig_dir, "fig_ambient_evidence_per_sample.pdf"), width = 14, height = 10)
print(p4_combined)
dev.off()
cat("  Saved: fig_ambient_evidence_per_sample.pdf\n")

# ==============================================================================
# Figure 5: Combined supplemental figure
# ==============================================================================
cat("\n--- Figure 5: Combined supplemental figure ---\n")

# Assemble available panels
panels <- list()
panel_labels <- c()

if (!is.null(p1)) { panels <- c(panels, list(p1)); panel_labels <- c(panel_labels, "A") }
if (!is.null(p2)) { panels <- c(panels, list(p2)); panel_labels <- c(panel_labels, "B") }
if (!is.null(p3)) { panels <- c(panels, list(p3)); panel_labels <- c(panel_labels, "C") }
panels <- c(panels, list(p4a + p4b))
panel_labels <- c(panel_labels, "D")

if (length(panels) >= 2) {
  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = "Supplemental: DecontX Ambient RNA Removal Evidence",
      tag_levels = list(panel_labels),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.tag = element_text(face = "bold", size = 14)
      )
    )

  height <- 5 * length(panels)
  pdf(file.path(fig_dir, "fig_ambient_evidence_combined.pdf"),
      width = 14, height = min(height, 30))
  print(combined)
  dev.off()
  cat("  Saved: fig_ambient_evidence_combined.pdf\n")
} else {
  cat("  Not enough panels for combined figure\n")
}

# ==============================================================================
# Save summary CSV
# ==============================================================================
cat("\nSaving ambient RNA summary CSV...\n")

summary_out <- list()
summary_out$per_sample <- sample_summary

if (exists("marker_summary") && !is.null(marker_summary)) {
  summary_out$marker_comparison <- as.data.frame(marker_summary)
}
if (exists("top_reduced_df") && !is.null(top_reduced_df)) {
  summary_out$top_reduced_genes <- top_reduced_df
}

# Write per-sample summary as main CSV
write.csv(sample_summary, file.path(BASE, "ambient_rna_summary.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: %s/ambient_rna_summary.csv\n", BASE))

# Write marker comparison if available
if (exists("marker_summary") && !is.null(marker_summary)) {
  write.csv(marker_summary,
            file.path(BASE, "ambient_rna_marker_comparison.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved: %s/ambient_rna_marker_comparison.csv\n", BASE))
}

# --- Summary ------------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(sprintf("  Total cells scored:    %d\n", nrow(contam_scores)))
cat(sprintf("  Matched to cell type:  %d\n", nrow(contam_with_type)))
cat(sprintf("  Eosinophils analyzed:  %d\n", ncol(eos)))
cat(sprintf("  Overall mean contam:   %.4f\n", mean(contam_scores$contamination)))
cat("\nOutput files:\n")
cat(sprintf("  %s/ambient_rna_summary.csv\n", BASE))
cat(sprintf("  %s/fig_ambient_evidence_*.pdf\n", fig_dir))
cat("\nDone! 21_ambient_rna_evidence.R completed successfully.\n")
