#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=eos_mac_combined
#SBATCH --output=eos_mac_combined_%j.out

# ==============================================================================
# 09_eos_mac_combined_umap.R
# Merge eosinophils and macrophages for a combined UMAP visualization
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Broad palette for eos vs mac overview
broad_colors <- c("Eosinophil" = "#E84A5F", "Macrophage" = "#2166AC")

# --- Load objects -------------------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
# Remove C6 if still present (too few cells, n=14)
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}

cat("Loading mac_subtyped.rds...\n")
mac <- readRDS(file.path(BASE, "mac_subtyped.rds"))

# --- Add celltype labels ------------------------------------------------------
cat("Adding detailed celltype labels...\n")

# Eosinophils: Eos_C0, Eos_C1, ...
eos$detailed_celltype <- paste0("Eos_", eos$eos_cluster)
eos$broad_celltype_combined <- "Eosinophil"

# Macrophages: Mac_M1, Mac_M2_Surf, Mac_M2_Func (already in mac_subtype)
mac$detailed_celltype <- mac$mac_subtype
mac$broad_celltype_combined <- "Macrophage"

cat(sprintf("  Eosinophils: %d cells, %d subtypes\n",
            ncol(eos), length(unique(eos$detailed_celltype))))
cat(sprintf("  Macrophages: %d cells, %d subtypes\n",
            ncol(mac), length(unique(mac$detailed_celltype))))

# --- Merge objects ------------------------------------------------------------
cat("Merging eosinophils and macrophages...\n")
combined <- merge(eos, mac)
cat(sprintf("  Combined: %d cells\n", ncol(combined)))

# --- Re-process combined object -----------------------------------------------
cat("Normalizing combined object...\n")
combined <- NormalizeData(combined, verbose = FALSE)

cat("Finding variable features...\n")
combined <- FindVariableFeatures(combined, nfeatures = 2000, verbose = FALSE)

cat("Scaling data...\n")
combined <- ScaleData(combined, verbose = FALSE)

cat("Running PCA...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

cat("Running UMAP...\n")
combined <- RunUMAP(combined, dims = 1:20, verbose = FALSE)

# --- Figures ------------------------------------------------------------------

# Filter detailed colors to only those present
present_types <- sort(unique(combined$detailed_celltype))
celltype_colors_use <- celltype_colors[intersect(names(celltype_colors), present_types)]

# UMAP by detailed celltype
cat("Plotting UMAP by detailed celltype...\n")
pdf(file.path(fig_dir, "09_combined_umap_detailed.pdf"), width = 9, height = 7)
print(
  DimPlot(combined, group.by = "detailed_celltype",
          cols = celltype_colors_use, pt.size = 0.3) +
    ggtitle("Eosinophil + Macrophage Combined UMAP (Detailed Subtypes)") +
    theme_umap()
)
dev.off()

# UMAP by broad celltype
cat("Plotting UMAP by broad celltype...\n")
pdf(file.path(fig_dir, "09_combined_umap_broad.pdf"), width = 7, height = 6)
print(
  DimPlot(combined, group.by = "broad_celltype_combined",
          cols = broad_colors, pt.size = 0.3) +
    ggtitle("Eosinophil vs Macrophage") +
    theme_umap()
)
dev.off()

# UMAP by condition
cat("Plotting UMAP by condition...\n")
pdf(file.path(fig_dir, "09_combined_umap_condition.pdf"), width = 8, height = 6)
print(
  DimPlot(combined, group.by = "condition", pt.size = 0.3) +
    ggtitle("Combined UMAP by Condition") +
    theme_umap()
)
dev.off()

# Composition stacked bar (detailed celltype by condition)
cat("Plotting composition stacked bar...\n")
comp_df <- combined@meta.data %>%
  group_by(condition, detailed_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(frac = n / sum(n))

pdf(file.path(fig_dir, "09_combined_composition.pdf"), width = 10, height = 6)
print(
  ggplot(comp_df, aes(x = condition, y = frac, fill = detailed_celltype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = celltype_colors_use) +
    labs(x = "Condition", y = "Fraction", fill = "Celltype",
         title = "Eosinophil + Macrophage Composition by Condition") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

# --- Save output --------------------------------------------------------------
cat("Saving eos_mac_combined.rds...\n")
saveRDS(combined, file.path(BASE, "eos_mac_combined.rds"))

cat("Done! 09_eos_mac_combined_umap.R completed successfully.\n")
