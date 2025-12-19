library(Seurat)
library(ggplot2)
library(dplyr)


# ======================
# Visualization: Violin plots for macrophages
# ======================

# Read the saved Seurat object with M1/M2 scores
seurat_obj <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/seurat_object_M1_M2_scored.rds")

fig_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Only plot for macrophages
seurat_obj$active.ident <- seurat_obj$cell_type
macrophages <- subset(seurat_obj, idents = "Macrophages")

# M1 violin plot
vln_m1 <- VlnPlot(macrophages, features = "M1_Score1", pt.size = 0) +
  ggtitle("M1 Signature Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(fig_dir, "violin_M1_score_macrophages.pdf"),
       plot = vln_m1, width = 8, height = 6)

# M2 violin plot
vln_m2 <- VlnPlot(macrophages, features = "M2_Score1", pt.size = 0) +
  ggtitle("M2 Signature Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(fig_dir, "violin_M2_score_macrophages.pdf"),
       plot = vln_m2, width = 8, height = 6)

cat("Scoring complete. Output saved.\n")


# ======================
# Visualization: Violin plots for all cell types
# ======================

# Reset identity to cell_type and plot all cell types
Idents(seurat_obj) <- "cell_type"

# M1 violin plot for all cell types
vln_m1_all <- VlnPlot(seurat_obj, features = "M1_Score1", pt.size = 0) +
  ggtitle("M1 Signature Score Across Cell Types") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M1_score_all_cells.pdf"),
       plot = vln_m1_all, width = 12, height = 6)

# M2 violin plot for all cell types  
vln_m2_all <- VlnPlot(seurat_obj, features = "M2_Score1", pt.size = 0) +
  ggtitle("M2 Signature Score Across Cell Types") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M2_score_all_cells.pdf"),
       plot = vln_m2_all, width = 12, height = 6)

# ======================
# umap plot for M1 and M2 scores
# ======================

# M1 score UMAP
umap_m1 <- FeaturePlot(seurat_obj, features = "M1_Score1",
                       reduction = "umap", pt.size = 0.3, order = TRUE) +
  ggtitle("UMAP of M1 Signature Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(fig_dir, "umap_M1_score.pdf"),
       plot = umap_m1, width = 8, height = 6)


# M2 score UMAP
umap_m2 <- FeaturePlot(seurat_obj, features = "M2_Score1",
                       reduction = "umap", pt.size = 0.3, order = TRUE) +
  ggtitle("UMAP of M2 Signature Score") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(fig_dir, "umap_M2_score.pdf"),
       plot = umap_m2, width = 8, height = 6)


# ======================
# Visualization: Violin plots for all conditions and all cells
# ======================
Idents(seurat_obj) <- "condition"
# M1 score violin plot for all conditions
vln_m1_all_conditions <- VlnPlot(seurat_obj, features = "M1_Score1", pt.size = 0) +
  ggtitle("M1 Signature Score Across Conditions All Cells") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M1_score_all_conditions_all_cells.pdf"),
       plot = vln_m1_all_conditions, width = 12, height = 6)

# M2 score violin plot for all conditions
vln_m2_all_conditions <- VlnPlot(seurat_obj, features = "M2_Score1", pt.size = 0) +
  ggtitle("M2 Signature Score Across Conditions All Cells") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M2_score_all_conditions_all_cells.pdf"),
       plot = vln_m2_all_conditions, width = 12, height = 6)

# ======================
# Visualization: Violin plots for all conditions and only macrophages
# ======================
Idents(seurat_obj) <- "condition"
macrophages <- subset(seurat_obj, cell_type == "Macrophages")
# M1 score violin plot for all conditions
vln_m1_all_conditions <- VlnPlot(macrophages, features = "M1_Score1", pt.size = 0) +
  ggtitle("M1 Signature Score Across Conditions Macrophages") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M1_score_all_conditions_macrophages.pdf"),
       plot = vln_m1_all_conditions, width = 12, height = 6)

# M2 score violin plot for all conditions
vln_m2_all_conditions <- VlnPlot(macrophages, features = "M2_Score1", pt.size = 0) +
  ggtitle("M2 Signature Score Across Conditions Macrophages") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(fig_dir, "violin_M2_score_all_conditions_macrophages.pdf"),
       plot = vln_m2_all_conditions, width = 12, height = 6)