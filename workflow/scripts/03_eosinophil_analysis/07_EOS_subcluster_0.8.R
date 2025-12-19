library(Seurat)
library(dplyr)
library(ggplot2)

dpi = 300

# Load seurat object
seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(seurat_file)

# 1. Subset the eosinophils
eos_subset <- subset(seurat_object, cell_type == "Eosinophils")

# 2. Normalize the data
eos_subset <- NormalizeData(eos_subset)

# 3. Find variable features
eos_subset <- FindVariableFeatures(eos_subset, selection.method = "vst", nfeatures = 2000)

# 4. Scale the data
all_genes <- rownames(eos_subset)
eos_subset <- ScaleData(eos_subset, features = all_genes)

# 5. Perform PCA
eos_subset <- RunPCA(eos_subset, features = VariableFeatures(object = eos_subset))

# 6. Determine the 'dimensionality' of the dataset
plot <- ElbowPlot(eos_subset, ndims = 50)
plot_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/ElbowPlot_EOS_0.8.png"
ggsave(plot_file, plot = plot, width = 4, height = 4)

# Based on the elbow plot, choose the number of dimensions for clustering
n_dims <- 20

# 7. Cluster the cells
eos_subset <- FindNeighbors(eos_subset, dims = 1:n_dims)
eos_subset <- FindClusters(eos_subset, resolution = 0.8)
# Resolution 0.8 yields 7 clusters (C0-C6)

# 8. Run UMAP for visualization
eos_subset <- RunUMAP(eos_subset, dims = 1:n_dims)

# Save the processed object
eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
saveRDS(eos_subset, file = eos_subset_file)

# Define cluster colors for eosinophil subtypes
eos_colors <- c(
    "0" = "#6CA6CD",
    "1" = "#E84A5F",
    "2" = "#FF6F1F",
    "3" = "#CBAACB",
    "4" = "#B03A48",
    "5" = "#5BA4A4",
    "6" = "#2F3C7E"
)

# Create UMAP plot by cluster
plot <- DimPlot(eos_subset, reduction = "umap", group.by = "seurat_clusters",
                label = TRUE, repel = TRUE, cols = eos_colors)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/EOS_UMAP_plot_by_cluster_0.8.png"
ggsave(plot_output_file, plot = plot, width = 10, height = 8)

# Create UMAP plot split by condition
umap_plot <- DimPlot(eos_subset,
                     reduction = "umap",
                     split.by = "condition",
                     ncol = 3,
                     label = TRUE,
                     label.size = 6,
                     pt.size = 1.5,
                     cols = eos_colors) +
                     ggtitle("Eosinophil Subclusters by Condition") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_condition_0.8.pdf",
       plot = umap_plot,
       width = 12,
       height = 8,
       device = cairo_pdf)

# Create UMAP plot split by Sex
umap_plot <- DimPlot(eos_subset,
                     reduction = "umap",
                     split.by = "Sex",
                     label = TRUE,
                     label.size = 6,
                     pt.size = 1.5,
                     cols = eos_colors) +
                     ggtitle("Eosinophil Subclusters by Sex") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_Sex_0.8.pdf",
       plot = umap_plot,
       width = 10,
       height = 8,
       device = cairo_pdf)

# Create UMAP plot split by genotype
umap_plot <- DimPlot(eos_subset,
                     reduction = "umap",
                     split.by = "genotype",
                     label = TRUE,
                     label.size = 6,
                     pt.size = 1.5,
                     cols = eos_colors) +
                     ggtitle("Eosinophil Subclusters by Genotype") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_genotype_0.8.pdf",
       plot = umap_plot,
       width = 10,
       height = 8,
       device = cairo_pdf)

# Create UMAP plot split by infection status
umap_plot <- DimPlot(eos_subset,
                     reduction = "umap",
                     split.by = "infection_status",
                     label = TRUE,
                     label.size = 6,
                     pt.size = 1.5,
                     cols = eos_colors) +
                     ggtitle("Eosinophil Subclusters by Infection Status") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_infection_status_0.8.pdf",
       plot = umap_plot,
       width = 10,
       height = 8,
       device = cairo_pdf)

# 9. Find markers for each cluster
eos_markers <- FindAllMarkers(eos_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save all markers
write.csv(eos_markers, "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8_all_markers.csv", row.names = FALSE)

# View top markers for each cluster
top_markers <- eos_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# Save top markers
write.csv(top_markers, "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8_top5_markers.csv", row.names = FALSE)
