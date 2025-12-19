# cd /bigdata/godziklab/shared/Xinru/302025


library(Seurat)
library(ggplot2)

singler_seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_singleR.rds"
seurat_object <- readRDS(file = singler_seurat_file)

### cluster the cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.8)

seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_singleR_clustered.rds"
saveRDS(seurat_object, file = seurat_file)

# Plot UMAP and color by orig.ident
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
umap_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_orig.ident.png"
ggsave(umap_output_file, plot = umap_plot, width = 10, height = 8)

# Plot UMAP and color by seurat_clusters
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
umap_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_seurat_clusters.png"
ggsave(umap_output_file, plot = umap_plot, width = 10, height = 8)

# Plot UMAP and color by SingleR.labels
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE)
umap_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_SingleR.labels.png"
ggsave(umap_output_file, plot = umap_plot, width = 10, height = 8)



####-----------------------------------use UMAP to find the cell types------------------------------------
library(cowplot)


seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(seurat_file)

# Create the UMAP plot split by cluster
plot <- DimPlot(seurat_object, reduction = "umap", split.by = "seurat_clusters", ncol = 5)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_cluster.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 30)

# Color by eos feature genes
plot <- FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = c("Siglecf", "Il5ra", "Ccr3", "Prg2", "Epx"),
            ncol = 3)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_EOS_feature.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 16)


### Count how many EOS under each condition
# Filter the metadata for seurat_clusters == 1
cluster_1_data <- seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == 1, ]

# Calculate the count of each orig.ident
orig_ident_count <- table(cluster_1_data$orig.ident)

# Convert to a data frame for better display
orig_ident_count_df <- as.data.frame(orig_ident_count)
colnames(orig_ident_count_df) <- c("orig.ident", "Count")


# Color by mac feature genes
plot <- FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = c("Ptprc", "Adgre1", "Itgam", "Fcgr1"),
            ncol = 2)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_MAC_feature.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 16)


# Color by microglia feature genes
# While there is overlap, microglia express some unique genes not found in macrophages, such as P2ry12, Tmem119, and Sall1
plot <- FeaturePlot(seurat_object, 
            reduction = "umap", 
            features = c("P2ry12", "Tmem119", "Sall1"),
            ncol = 2)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_microglia_feature.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 16)

# Create the UMAP plot split by condition
plot <- DimPlot(seurat_object, reduction = "umap", split.by = "condition", ncol = 3,label = TRUE, repel = TRUE)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_condition.png"
ggsave(plot_output_file, plot = plot, width = 16, height = 16)

# Create the UMAP plot split by Sex
plot <- DimPlot(seurat_object, reduction = "umap", split.by = "Sex",label = TRUE, repel = TRUE)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_Sex.png"
ggsave(plot_output_file, plot = plot, width = 8, height = 8)

# Create the UMAP plot split by genotype
plot <- DimPlot(seurat_object, reduction = "umap", split.by = "genotype",label = TRUE, repel = TRUE)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_genotype.png"
ggsave(plot_output_file, plot = plot, width = 8, height = 8)

# Create the UMAP plot split by infect
plot <- DimPlot(seurat_object, reduction = "umap", split.by = "infection_status",label = TRUE, repel = TRUE)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_infection_status.png"
ggsave(plot_output_file, plot = plot, width = 8, height = 8)


                          
