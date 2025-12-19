lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(Matrix)


combined <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/eos_GSM7919060_GSE248899_GSM7156898.rds")
# Print dimensions of count matrices
print("Dimensions of counts.1:")
print(dim(combined[["RNA"]]$counts.1))
# [1] 19059 10698
print("Dimensions of counts.2:")
# [1] 11681  1086
print(dim(combined[["RNA"]]$counts.2))
print("Dimensions of counts.3:")
print(dim(combined[["RNA"]]$counts.3))
# [1] 32344  5679
print("Dimensions of counts.4:")
print(dim(combined[["RNA"]]$counts.4))
# [1] 14525 75099

# Function to process each count matrix separately
process_counts <- function(counts, assay_name) {
  seurat_obj <- CreateSeuratObject(counts = counts, project = assay_name, min.cells = 3, min.features = 200)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
}

# Process each count matrix
seurat_1 <- process_counts(combined[["RNA"]]$counts.1, "Eos")
seurat_2 <- process_counts(combined[["RNA"]]$counts.2, "GSM7919060")
seurat_3 <- process_counts(combined[["RNA"]]$counts.3, "GSE248899")
seurat_4 <- process_counts(combined[["RNA"]]$counts.4, "GSM7156898")

# Add dataset-specific labels
seurat_1$dataset <- "Eos"
seurat_2$dataset <- "GSM7919060"
seurat_3$dataset <- "GSE248899"
seurat_4$dataset <- "GSM7156898"

# Merge the Seurat objects
merged_seurat <- merge(seurat_1, y = c(seurat_2, seurat_3, seurat_4), 
                       add.cell.ids = c("Eos", "GSM7919060", "GSE248899", "GSM7156898"), 
                       project = "merged")

# Find variable features in the merged object
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# Scale the data
merged_seurat <- ScaleData(merged_seurat, features = rownames(merged_seurat))

# Run PCA
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(merged_seurat))

# Run Harmony
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "dataset")

# Run UMAP and clustering on Harmony embeddings
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# add eos_sub meta data
library(dplyr)
eos_sub_metadata <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_sub_metadata.csv", row.names = 1)  # Assuming cell names are in the first column
eos_sub_metadata <- eos_sub_metadata %>% select(cell_type, Clusters, infection_status, genotype, condition, Sex)
eos_sub_metadata$cell_id <- paste0("Eos_", rownames(eos_sub_metadata), "_1")

# Add row names as a column in both data frames
metadata <- merged_seurat@meta.data
metadata$cell_id <- rownames(metadata)

# Now perform the left join
metadata <- left_join(metadata, eos_sub_metadata, by = "cell_id")
row.names(metadata) <- metadata$cell_id
metadata$cell_id <- NULL
# If you want to update the Seurat object with this new metadata
merged_seurat@meta.data <- metadata
merged_seurat@meta.data$Clusters <- ifelse(is.na(merged_seurat@meta.data$Clusters), merged_seurat@meta.data$dataset, merged_seurat@meta.data$Clusters)
merged_seurat@meta.data$condition <- ifelse(is.na(merged_seurat@meta.data$condition), merged_seurat@meta.data$dataset, merged_seurat@meta.data$condition)

# Save the integrated object
saveRDS(merged_seurat, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/harmony_integrated_eos_GSM7919060_GSE248899_GSM7156898.rds")

# Print a summary of the object
print(merged_seurat)

# Print the number of cells from each original assay
table(merged_seurat$dataset)
table(merged_seurat$Clusters)
table(merged_seurat$condition)


#> table(merged_seurat$dataset)

#       Eos  GSE248899 GSM7156898 GSM7919060
#     10698       1006      75097       1086


library(Seurat)
library(ggplot2)
library(cowplot)

setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/")

# Custom function to create UMAP plot
create_umap_plot <- function(seurat_obj, group.by = "dataset", title = NULL) {
  # Extract UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, "umap")
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    UMAP_1 = umap_coords[,1],
    UMAP_2 = umap_coords[,2],
    Group = seurat_obj@meta.data[[group.by]]
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_cowplot() +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  
  return(p)
}

# Create UMAP plot
p1 <- create_umap_plot(merged_seurat, group.by = "dataset", title = "UMAP by Dataset")

# Save the plot
ggsave("immAll_UMAP_plot_by_dataset.png", plot = p1, width = 10, height = 8, units = "in", dpi = 300)


# Create and save a plot by clusters
p2 <- create_umap_plot(merged_seurat, group.by = "Clusters", title = "UMAP by Clusters")
ggsave("immAll_UMAP_plot_by_clusters.png", plot = p2, width = 10, height = 8, units = "in", dpi = 300)

p3 <- create_umap_plot(merged_seurat, group.by = "condition", title = "UMAP by condition")
ggsave("immAll_UMAP_plot_by_condition.png", plot = p3, width = 10, height = 8, units = "in", dpi = 300)

