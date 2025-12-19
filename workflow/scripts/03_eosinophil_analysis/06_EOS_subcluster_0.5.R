library(Seurat)
library(dplyr)
library(ggplot2)

dpi = 300
# 1. Subset the eosinophils
# Assuming your eosinophils are labeled as "Eosinophils" in the cell type column
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
# Adjust the number of dimensions based on your data

plot <-  ElbowPlot(eos_subset, ndims = 50)
plot_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/ElbowPlot_EOS.png"
ggsave(plot_file, plot = plot, width = 4, height = 4)

# Based on the elbow plot, choose the number of dimensions for clustering
n_dims <- 20  # Adjust this number based on your elbow plot

# 7. Cluster the cells
eos_subset <- FindNeighbors(eos_subset, dims = 1:n_dims)

# eos_subset <- FindClusters(eos_subset, resolution = 0.1)  # Adjust resolution as needed
#> table(eos_subset@meta.data$RNA_snn_res.0.1)
#   0    1    2
#9555  875  268

eos_subset <- FindClusters(eos_subset, resolution = 0.5)  # Adjust resolution as needed
#   0    1    2    3    4
# 4389 2637 2316 1088  268

# 8. Run UMAP for visualization
eos_subset <- RunUMAP(eos_subset, dims = 1:n_dims)

eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.5.rds"
saveRDS(eos_subset, file = eos_subset_file)


# Create the UMAP plot split by cluster (RNA_snn_res.0.5)
plot <- DimPlot(eos_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/EOS_UMAP_plot_by_cluster_0.5.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 16)

# Create the UMAP plot split by condition
umap_plot <- DimPlot(eos_subset, 
                     reduction = "umap", 
                     split.by = "condition",
                     ncol = 3,
                     label = TRUE, 
                     label.size = 6,
                     pt.size = 1.5) +  # Increase this value to make dots bigger
                     ggtitle("Eosinophil") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_condition_0.5.pdf",
       plot = umap_plot,
       width = 8,
       height = 8,
       device = cairo_pdf)



# Create the UMAP plot split by Sex
umap_plot <- DimPlot(eos_subset, 
                     reduction = "umap", 
                     split.by = "Sex",
                     label = TRUE, 
                     label.size = 6,
                     pt.size = 1.5) +  # Increase this value to make dots bigger
                     ggtitle("Eosinophil") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_Sex_0.5.pdf",
       plot = umap_plot,
       width = 8,
       height = 8,
       device = cairo_pdf)


# Create the UMAP plot split by genotype
umap_plot <- DimPlot(eos_subset, 
                     reduction = "umap", 
                     split.by = "genotype",
                     label = TRUE, 
                     label.size = 6,
                     pt.size = 1.5) +  # Increase this value to make dots bigger
                     ggtitle("Eosinophil") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_genotype_0.5.pdf",
       plot = umap_plot,
       width = 8,
       height = 8,
       device = cairo_pdf)

# Create the UMAP plot split by infect
umap_plot <- DimPlot(eos_subset, 
                     reduction = "umap", 
                     split.by = "infection_status",
                     label = TRUE, 
                     label.size = 6,
                     pt.size = 1.5) +  # Increase this value to make dots bigger
                     ggtitle("Eosinophil") +
                     theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                           legend.text = element_text(size = 12),
                           axis.title = element_text(size = 14))
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_infection_status_0.5.pdf",
       plot = umap_plot,
       width = 8,
       height = 8,
       device = cairo_pdf)


# 10. Find markers for each cluster
eos_markers <- FindAllMarkers(eos_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers for each cluster
top_markers <- eos_markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)

print(top_markers)

