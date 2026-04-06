#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 30_cross_dataset_similarity.R
# Cross-dataset cluster similarity: Qiu C0-C5 vs Hu depots vs Gurtner BM/Blood
# Method: Spearman correlation of average gene expression profiles
# Input:  three_way_integrated.rds + hu_adipose_eos.rds (for depot info)
# Output: similarity heatmap PDF + CSV
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "integration", "similarity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load three-way integrated object
# ==============================================================================
cat("Loading three-way integrated object...\n")
tw <- readRDS(file.path(BASE, "integration", "three_way", "three_way_integrated.rds"))
if (inherits(tw[["RNA"]], "Assay5")) tw[["RNA"]] <- JoinLayers(tw[["RNA"]])
cat(sprintf("  %d cells\n", ncol(tw)))

# ==============================================================================
# 2. Recover Hu depot info from original object
# ==============================================================================
cat("Loading original Hu object for depot info...\n")
hu <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds")

# Find depot column — search for column with BAT/eWAT/scWAT values
depot_col <- NULL
for (col in colnames(hu@meta.data)) {
  vals <- as.character(unique(hu@meta.data[[col]]))
  if (any(c("BAT", "eWAT", "scWAT") %in% vals)) {
    depot_col <- col
    break
  }
}

if (!is.null(depot_col)) {
  cat(sprintf("  Found depot column: '%s'\n", depot_col))
  hu_depot <- data.frame(
    hu_barcode = colnames(hu),
    depot = as.character(hu@meta.data[[depot_col]]),
    stringsAsFactors = FALSE
  )
  print(table(hu_depot$depot))
} else {
  cat("  WARNING: No depot column found in Hu object\n")
  hu_depot <- data.frame(
    hu_barcode = colnames(hu),
    depot = "Adipose_Healthy",
    stringsAsFactors = FALSE
  )
}
rm(hu); gc()

# Map Hu depot to three-way barcodes
# Barcodes in the three-way object may have various prefixes from Seurat merge
hu_barcodes_tw <- colnames(tw)[tw$dataset == "Hu_2026_Healthy"]
hu_barcodes_orig <- hu_depot$hu_barcode
cat(sprintf("  Hu cells in three-way: %d\n", length(hu_barcodes_tw)))
cat(sprintf("  Hu cells in original: %d\n", length(hu_barcodes_orig)))
cat("  Example three-way barcodes:", paste(head(hu_barcodes_tw, 3), collapse=", "), "\n")
cat("  Example original barcodes:", paste(head(hu_barcodes_orig, 3), collapse=", "), "\n")

# Try multiple matching strategies
depot_map <- setNames(rep(NA_character_, ncol(tw)), colnames(tw))

# Strategy 1: exact match with "Hu_" prefix
match1 <- match(hu_barcodes_tw, paste0("Hu_", hu_barcodes_orig))
n1 <- sum(!is.na(match1))
cat(sprintf("  Strategy 1 (Hu_ prefix): %d matches\n", n1))
if (n1 > 0) {
  depot_map[hu_barcodes_tw[!is.na(match1)]] <- hu_depot$depot[match1[!is.na(match1)]]
}

# Strategy 2: strip prefix from three-way barcodes, match directly
if (n1 < length(hu_barcodes_tw) * 0.5) {
  stripped <- sub("^Hu_", "", hu_barcodes_tw)
  match2 <- match(stripped, hu_barcodes_orig)
  n2 <- sum(!is.na(match2))
  cat(sprintf("  Strategy 2 (strip Hu_): %d matches\n", n2))
  if (n2 > n1) {
    depot_map[hu_barcodes_tw[!is.na(match2)]] <- hu_depot$depot[match2[!is.na(match2)]]
  }
}

# Strategy 3: strip all prefixes up to last underscore before barcode
if (sum(!is.na(depot_map[hu_barcodes_tw])) < length(hu_barcodes_tw) * 0.5) {
  # Try matching by barcode suffix (e.g., ACGTACGT-1)
  tw_suffix <- sub("^.*_([ACGT]+-\\d+)$", "\\1", hu_barcodes_tw)
  orig_suffix <- sub("^.*_([ACGT]+-\\d+)$", "\\1", hu_barcodes_orig)
  match3 <- match(tw_suffix, orig_suffix)
  n3 <- sum(!is.na(match3))
  cat(sprintf("  Strategy 3 (barcode suffix): %d matches\n", n3))
  if (n3 > sum(!is.na(depot_map[hu_barcodes_tw]))) {
    depot_map[hu_barcodes_tw[!is.na(match3)]] <- hu_depot$depot[match3[!is.na(match3)]]
  }
}

matched_hu <- sum(!is.na(depot_map[hu_barcodes_tw]))
cat(sprintf("  Final Hu depot matches: %d / %d\n", matched_hu, length(hu_barcodes_tw)))
cat("  Depot distribution:\n")
print(table(depot_map[hu_barcodes_tw], useNA = "ifany"))

# ==============================================================================
# 3. Assign comparison groups (direct assignment — avoids case_when factor issues)
# ==============================================================================
cat("Assigning comparison groups...\n")

ds <- as.character(tw$dataset)
org <- as.character(tw$organ)
cl <- as.character(tw$original_cluster)

tw$sim_group <- NA_character_

# Qiu: use original cluster (C0-C5)
qiu_idx <- which(ds == "Qiu_2025_Adipose")
tw$sim_group[qiu_idx] <- cl[qiu_idx]
cat(sprintf("  Qiu assigned: %d\n", length(qiu_idx)))

# Hu: use depot (BAT, eWAT, scWAT)
hu_idx <- which(ds == "Hu_2026_Healthy")
tw$sim_group[hu_idx] <- paste0("Hu.", depot_map[colnames(tw)[hu_idx]])
cat(sprintf("  Hu assigned: %d (example: %s)\n", length(hu_idx),
            paste(head(tw$sim_group[hu_idx], 3), collapse=", ")))

# Gurtner BM and Blood
bm_idx <- which(org == "Bone Marrow")
tw$sim_group[bm_idx] <- "Gurtner.BM"
cat(sprintf("  Gurtner BM assigned: %d\n", length(bm_idx)))

blood_idx <- which(org == "Blood")
tw$sim_group[blood_idx] <- "Gurtner.Blood"
cat(sprintf("  Gurtner Blood assigned: %d\n", length(blood_idx)))

# Keep only groups we want
keep_groups <- c(paste0("C", 0:5), "Hu.eWAT", "Hu.scWAT", "Hu.BAT",
                 "Gurtner.BM", "Gurtner.Blood")
tw$sim_group[!tw$sim_group %in% keep_groups] <- NA

cat("  Group sizes:\n")
print(table(tw$sim_group, useNA = "ifany"))

# Subset to valid groups
tw_sub <- subset(tw, cells = colnames(tw)[!is.na(tw$sim_group)])
tw_sub$sim_group <- factor(tw_sub$sim_group, levels = keep_groups)
Idents(tw_sub) <- "sim_group"
cat(sprintf("  Keeping %d cells in %d groups\n", ncol(tw_sub), length(keep_groups)))

# ==============================================================================
# 4. Compute Spearman correlation of average expression
# ==============================================================================
cat("Computing average expression...\n")
tw_sub <- FindVariableFeatures(tw_sub, nfeatures = 2000, verbose = FALSE)
hvg <- VariableFeatures(tw_sub)

avg <- AverageExpression(tw_sub, features = hvg, assays = "RNA", layer = "data")$RNA
avg <- as.matrix(avg)
cat(sprintf("  Average expression matrix: %d genes x %d groups\n", nrow(avg), ncol(avg)))
cat("  Columns:", paste(colnames(avg), collapse = ", "), "\n")

cor_mat <- cor(avg, method = "spearman")
write.csv(cor_mat, file.path(outdir, "cross_dataset_spearman_correlation.csv"))

# ==============================================================================
# 5. Heatmap
# ==============================================================================
cat("Plotting similarity heatmap...\n")

# Annotation
ann_row <- data.frame(
  Source = case_when(
    grepl("^C[0-5]", rownames(cor_mat)) ~ "Qiu HFD",
    grepl("^Hu", rownames(cor_mat))     ~ "Hu Healthy",
    TRUE                                 ~ "Gurtner"
  ),
  row.names = rownames(cor_mat)
)
ann_colors <- list(
  Source = c("Qiu HFD" = "#E63946", "Hu Healthy" = "#4DAF4A", "Gurtner" = "#1D3557")
)

# Display labels (prettier names)
display_names <- c(
  "C0" = "C0", "C1" = "C1", "C2" = "C2", "C3" = "C3", "C4" = "C4", "C5" = "C5",
  "Hu.eWAT" = "Hu eWAT", "Hu.scWAT" = "Hu scWAT", "Hu.BAT" = "Hu BAT",
  "Gurtner.BM" = "Gurtner BM", "Gurtner.Blood" = "Gurtner Blood"
)
row_labels <- display_names[rownames(cor_mat)]
row_labels[is.na(row_labels)] <- rownames(cor_mat)[is.na(row_labels)]

pdf(file.path(outdir, "cross_dataset_similarity_heatmap.pdf"), width = 10, height = 8)
pheatmap(cor_mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(0.5, 1, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_row = ann_row, annotation_col = ann_row,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 8,
         fontsize_row = 11, fontsize_col = 11,
         labels_row = row_labels, labels_col = row_labels,
         main = "Cross-Dataset Eosinophil Similarity (Spearman)",
         border_color = NA)
dev.off()
cat("  Saved: cross_dataset_similarity_heatmap.pdf\n")

cat("\n========== Cross-dataset similarity analysis complete ==========\n")
cat("Output directory:", outdir, "\n")
