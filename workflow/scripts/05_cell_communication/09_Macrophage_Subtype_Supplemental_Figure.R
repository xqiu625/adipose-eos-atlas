# Macrophage Subtype Supplemental Figure
# Purpose: Generate UMAP and DotPlot showing M1, M2 Functional, M2 Surface markers
# Date: 2024-12-08
#
# Output:
#   1. UMAP of macrophage subtypes (excluding Other Macs)
#   2. DotPlot of marker genes for each subtype

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# =============================================================================
# SETUP
# =============================================================================

DATA_DIR <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed"
OUTPUT_DIR <- "/bigdata/nairlab/shared/XQiu/EOS_project/R20_20251208/macrophage_supplemental"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_object_with_mac_subtypes_251027.rds"))

# Extract only macrophages
macrophages <- subset(seurat_object, subset = cell_type == "Macrophages")

cat("Total macrophages:", ncol(macrophages), "\n")
cat("Subtypes present:\n")
print(table(macrophages$macrophage_type))

# =============================================================================
# FILTER: Remove "Other_Macrophages"
# =============================================================================

cat("\nFiltering out Other_Macrophages...\n")

# Keep only M1, M2_Functional, M2_Surface
mac_filtered <- subset(macrophages,
                       subset = macrophage_type %in% c("M1_Macrophages",
                                                        "M2_Macrophages_Functional",
                                                        "M2_Macrophages_Surface_Only"))

cat("Filtered macrophages:", ncol(mac_filtered), "\n")
cat("Subtypes after filtering:\n")
print(table(mac_filtered$macrophage_type))

# Rename for cleaner labels
mac_filtered$macrophage_subtype <- case_when(
  mac_filtered$macrophage_type == "M1_Macrophages" ~ "M1",
  mac_filtered$macrophage_type == "M2_Macrophages_Functional" ~ "M2 Functional",
  mac_filtered$macrophage_type == "M2_Macrophages_Surface_Only" ~ "M2 Surface",
  TRUE ~ mac_filtered$macrophage_type
)

# Set factor levels for ordering
mac_filtered$macrophage_subtype <- factor(mac_filtered$macrophage_subtype,
                                           levels = c("M1", "M2 Functional", "M2 Surface"))

# Set identity
Idents(mac_filtered) <- "macrophage_subtype"

# =============================================================================
# DEFINE COLORS (as specified)
# =============================================================================

mac_colors <- c(
  "M1" = "#cc426c",
  "M2 Functional" = "#8fc43e",
  "M2 Surface" = "#58429a"
)

# =============================================================================
# DEFINE MARKER GENES
# =============================================================================

# Genes requested (with gene symbol mapping)
# Cd11c = Itgax
marker_genes <- c(
  # M1 markers
  "Itgax",    # CD11c
  "Cd80",
  "Tnf",
  "Il6",
  "Ccl2",
  # M2 Functional markers
  "Il10",
  "Arg1",
  "Chil3",    # Cd301 = Clec10a, but Chil3 is the functional marker
  # M2 Surface markers
  "Mrc1",
  "Lyve1",
  "Cd163"
)

# Full list as requested by user
marker_genes_full <- c(
  "Itgax",    # CD11c
  "Cd80",
  "Tnf",
  "Il6",
  "Ccl2",
  "Il10",
  "Arg1",
  "Clec10a",  # CD301
  "Chil3",
  "Mrc1",
  "Lyve1",
  "Cd163"
)

# Check which genes are present
genes_present <- marker_genes_full[marker_genes_full %in% rownames(mac_filtered)]
genes_missing <- marker_genes_full[!marker_genes_full %in% rownames(mac_filtered)]

cat("\nMarker genes present:", length(genes_present), "\n")
print(genes_present)
if (length(genes_missing) > 0) {
  cat("\nMarker genes missing:", length(genes_missing), "\n")
  print(genes_missing)
}

# =============================================================================
# FIGURE 1: UMAP
# =============================================================================

cat("\nCreating UMAP plot...\n")

umap_plot <- DimPlot(
  mac_filtered,
  reduction = "umap",
  group.by = "macrophage_subtype",
  cols = mac_colors,
  pt.size = 0.5,
  label = FALSE
) +
  ggtitle("Macrophage Subtypes") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Save UMAP
ggsave(
  file.path(OUTPUT_DIR, "macrophage_subtypes_umap.pdf"),
  umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(OUTPUT_DIR, "macrophage_subtypes_umap.png"),
  umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("UMAP saved!\n")

# =============================================================================
# FIGURE 2: DOTPLOT
# =============================================================================

cat("\nCreating DotPlot...\n")

# Order genes as specified by user:
# Cd11c (Itgax), Cd80, Tnf, Il6, Ccl2, Il10, Arg1, Cd301, Chil3, Mrc1, Lyve1, Cd301b, Cd163
# Note: Cd301 = Clec10a, Cd301b = Mgl2

genes_for_dotplot <- c("Itgax", "Cd80", "Tnf", "Il6", "Ccl2",
                        "Il10", "Arg1", "Clec10a", "Chil3",
                        "Mrc1", "Lyve1", "Mgl2", "Cd163")

# Filter to only present genes
genes_for_dotplot <- genes_for_dotplot[genes_for_dotplot %in% rownames(mac_filtered)]

cat("Genes for DotPlot:", length(genes_for_dotplot), "\n")
print(genes_for_dotplot)

# Create DotPlot
dot_plot <- DotPlot(
  mac_filtered,
  features = genes_for_dotplot,
  cols = c("lightgrey", "darkred"),
  dot.scale = 8,
  group.by = "macrophage_subtype"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  ggtitle("Marker Gene Expression")

# Save DotPlot
ggsave(
  file.path(OUTPUT_DIR, "macrophage_markers_dotplot.pdf"),
  dot_plot,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(OUTPUT_DIR, "macrophage_markers_dotplot.png"),
  dot_plot,
  width = 10,
  height = 4,
  dpi = 300
)

cat("DotPlot saved!\n")

# =============================================================================
# FIGURE 3: COMBINED SUPPLEMENTAL FIGURE
# =============================================================================

cat("\nCreating combined supplemental figure...\n")

# Combine UMAP and DotPlot
combined_plot <- umap_plot / dot_plot +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(
    title = "Macrophage Subtype Classification",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# Save combined figure
ggsave(
  file.path(OUTPUT_DIR, "macrophage_supplemental_figure.pdf"),
  combined_plot,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(OUTPUT_DIR, "macrophage_supplemental_figure.png"),
  combined_plot,
  width = 10,
  height = 10,
  dpi = 300
)

cat("Combined supplemental figure saved!\n")

# =============================================================================
# ALTERNATIVE: HORIZONTAL LAYOUT
# =============================================================================

cat("\nCreating horizontal layout version...\n")

# Horizontal layout
combined_horizontal <- (umap_plot | dot_plot) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(
  file.path(OUTPUT_DIR, "macrophage_supplemental_figure_horizontal.pdf"),
  combined_horizontal,
  width = 14,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(OUTPUT_DIR, "macrophage_supplemental_figure_horizontal.png"),
  combined_horizontal,
  width = 14,
  height = 5,
  dpi = 300
)

cat("Horizontal layout saved!\n")

# =============================================================================
# CELL COUNT SUMMARY
# =============================================================================

cat("\n=== CELL COUNT SUMMARY ===\n")
cell_counts <- as.data.frame(table(mac_filtered$macrophage_subtype))
colnames(cell_counts) <- c("Subtype", "Cell_Count")
cell_counts$Percentage <- round(cell_counts$Cell_Count / sum(cell_counts$Cell_Count) * 100, 1)
print(cell_counts)

write.csv(cell_counts,
          file.path(OUTPUT_DIR, "macrophage_subtype_cell_counts.csv"),
          row.names = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Generated files:\n")
cat("  1. macrophage_subtypes_umap.pdf/png\n")
cat("  2. macrophage_markers_dotplot.pdf/png\n")
cat("  3. macrophage_supplemental_figure.pdf/png (vertical layout)\n")
cat("  4. macrophage_supplemental_figure_horizontal.pdf/png (horizontal layout)\n")
cat("  5. macrophage_subtype_cell_counts.csv\n")
cat("\nColors used:\n")
cat("  M1 Mac: #cc426c\n")
cat("  M2 Functional: #8fc43e\n")
cat("  M2 Surface: #58429a\n")

# SCP command
cat("\n=== SCP COMMAND ===\n")
cat(paste0('scp -r xinruqiu@cluster.hpcc.ucr.edu:', OUTPUT_DIR, ' .\n'))
