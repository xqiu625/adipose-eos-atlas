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

# Gurtner: use published cluster labels (Gurtner et al. 2022, Nature)
# Mapping: 4=Precursors, 3=Immature, 16=Cycling, 2=Circulating, 1=A-Eos, 5=B-Eos
gurtner_label_map <- c(
  "4"  = "Gurtner.Precursors",
  "3"  = "Gurtner.Immature",
  "16" = "Gurtner.Cycling",
  "2"  = "Gurtner.Circulating",
  "1"  = "Gurtner.A-Eos",
  "5"  = "Gurtner.B-Eos"
)
gurtner_idx <- which(ds == "Gurtner_2022_Blood_BM")
gurtner_cl <- cl[gurtner_idx]
tw$sim_group[gurtner_idx] <- gurtner_label_map[gurtner_cl]
cat(sprintf("  Gurtner assigned: %d\n", sum(!is.na(tw$sim_group[gurtner_idx]))))

# Keep only relevant groups (BM + Blood + intestinal + Hu depots + C0-C5)
keep_groups <- c(paste0("C", 0:5),
                 "Hu.eWAT", "Hu.scWAT", "Hu.BAT",
                 "Gurtner.Precursors", "Gurtner.Immature",
                 "Gurtner.Circulating", "Gurtner.A-Eos", "Gurtner.B-Eos")
tw$sim_group[!tw$sim_group %in% keep_groups] <- NA
cat(sprintf("  Keeping %d groups\n", length(keep_groups)))

cat("  Group sizes:\n")
print(table(tw$sim_group, useNA = "ifany"))

# Subset to valid groups
tw_sub <- subset(tw, cells = colnames(tw)[!is.na(tw$sim_group)])
tw_sub$sim_group <- factor(tw_sub$sim_group, levels = keep_groups)
Idents(tw_sub) <- "sim_group"
cat(sprintf("  Keeping %d cells in %d groups\n", ncol(tw_sub), length(keep_groups)))

# ==============================================================================
# 4. Compute similarity using Harmony embeddings (batch-corrected)
# ==============================================================================
cat("Computing average Harmony embeddings per group...\n")

# Use Harmony embeddings (batch-corrected) instead of raw expression
harmony_emb <- Embeddings(tw_sub, "harmony")
cat(sprintf("  Harmony embeddings: %d cells x %d dims\n", nrow(harmony_emb), ncol(harmony_emb)))

# Average Harmony embedding per group
groups <- as.character(Idents(tw_sub))
unique_groups <- levels(Idents(tw_sub))
unique_groups <- unique_groups[unique_groups %in% unique(groups)]

avg_harmony <- do.call(rbind, lapply(unique_groups, function(g) {
  colMeans(harmony_emb[groups == g, , drop = FALSE])
}))
rownames(avg_harmony) <- unique_groups
cat(sprintf("  Average embedding matrix: %d groups x %d dims\n", nrow(avg_harmony), ncol(avg_harmony)))
cat("  Groups:", paste(rownames(avg_harmony), collapse = ", "), "\n")

# Spearman correlation on Harmony embeddings
cor_mat <- cor(t(avg_harmony), method = "spearman")
write.csv(cor_mat, file.path(outdir, "cross_dataset_spearman_correlation_harmony.csv"))

# ==============================================================================
# 5. Dendrogram (Hu-paper style) + heatmap
# ==============================================================================

# Display labels — auto-generate from group names
make_display <- function(x) gsub("\\.", " ", x)

# Label colors by dataset source
get_label_color <- function(x) {
  ifelse(grepl("^C[0-5]$", x), "#E63946",
  ifelse(grepl("^Hu", x), "#4DAF4A", "#1D3557"))
}

# Hierarchical clustering on correlation distance
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")

# --- Dendrogram (main figure, Hu-paper style) ---
cat("Plotting dendrogram...\n")
library(dendextend)

dend <- as.dendrogram(hc)

# Map labels to display names and colors
orig_labels <- labels(dend)
labels(dend) <- make_display(orig_labels)
labels_colors(dend) <- sapply(orig_labels, get_label_color)
labels_cex(dend) <- 1.1

pdf(file.path(outdir, "cross_dataset_dendrogram.pdf"), width = 12, height = 7)
par(mar = c(8, 4, 4, 2))
plot(dend, main = "Eosinophil scRNA-seq Cross-Dataset Similarity",
     ylab = "1 - Spearman correlation", las = 2)
# Legend
legend("topright",
       legend = c("Qiu et al. 2025 (HFD adipose)",
                  "Hu et al. 2026 (healthy adipose)",
                  "Gurtner et al. 2022 (BM/Blood/GI)"),
       text.col = c("#E63946", "#4DAF4A", "#1D3557"),
       bty = "n", cex = 1.0)
dev.off()
cat("  Saved: cross_dataset_dendrogram.pdf\n")

# --- Also save heatmap as supplement ---
cat("Plotting similarity heatmap (supplement)...\n")
ann_row <- data.frame(
  Source = ifelse(grepl("^C[0-5]$", rownames(cor_mat)), "Qiu HFD",
           ifelse(grepl("^Hu", rownames(cor_mat)), "Hu Healthy", "Gurtner")),
  row.names = rownames(cor_mat)
)
ann_colors <- list(
  Source = c("Qiu HFD" = "#E63946", "Hu Healthy" = "#4DAF4A", "Gurtner" = "#1D3557")
)
row_labels <- make_display(rownames(cor_mat))

n_groups <- ncol(cor_mat)
pdf(file.path(outdir, "cross_dataset_similarity_heatmap.pdf"),
    width = max(10, n_groups * 0.8), height = max(8, n_groups * 0.7))
pheatmap(cor_mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(0.5, 1, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_row = ann_row, annotation_col = ann_row,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = max(5, 10 - n_groups * 0.3),
         fontsize_row = 10, fontsize_col = 10,
         labels_row = row_labels, labels_col = row_labels,
         main = "Cross-Dataset Eosinophil Similarity (Spearman)",
         border_color = NA)
dev.off()
cat("  Saved: cross_dataset_similarity_heatmap.pdf\n")

cat("\n========== Cross-dataset similarity analysis complete ==========\n")
cat("Output directory:", outdir, "\n")
