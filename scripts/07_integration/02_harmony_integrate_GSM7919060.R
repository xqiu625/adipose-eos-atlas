# Load required libraries
library(Seurat)
library(harmony)
library(Matrix)


combined <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/eos_GSM7919060.rds")
# Print dimensions of count matrices
print("Dimensions of counts.1:")
print(dim(combined[["RNA"]]$counts.1))
print("Dimensions of counts.2:")
print(dim(combined[["RNA"]]$counts.2))

# Function to process each count matrix separately
process_counts <- function(counts, assay_name) {
  seurat_obj <- CreateSeuratObject(counts = counts, project = assay_name, min.cells = 3, min.features = 200)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
}

# Process each count matrix
seurat_1 <- process_counts(combined[["RNA"]]$counts.1, "Assay1")
seurat_2 <- process_counts(combined[["RNA"]]$counts.2, "Assay2")

# Find common variable features
features_1 <- VariableFeatures(seurat_1)
features_2 <- VariableFeatures(seurat_2)
common_features <- intersect(features_1, features_2)

# Merge the two Seurat objects
merged_seurat <- merge(seurat_1, seurat_2)

# Scale the data using common features
merged_seurat <- ScaleData(merged_seurat, features = common_features)

# Run PCA
merged_seurat <- RunPCA(merged_seurat, features = common_features)

# Run Harmony
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident")

# Run UMAP and clustering on Harmony embeddings
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# add eos_sub meta data
library(dplyr)
eos_sub_metadata <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_sub_metadata.csv", row.names = 1)  # Assuming cell names are in the first column
eos_sub_metadata <- eos_sub_metadata %>% select(cell_type, Clusters, infection_status, genotype, condition, Sex)
eos_sub_metadata$cell_id <- rownames(eos_sub_metadata)

# Add row names as a column in both data frames
metadata <- merged_seurat@meta.data
metadata$cell_id <- rownames(metadata)

# Now perform the left join
metadata <- left_join(metadata, eos_sub_metadata, by = "cell_id")
row.names(metadata) <- metadata$cell_id
metadata$cell_id <- NULL
# If you want to update the Seurat object with this new metadata
merged_seurat@meta.data <- metadata
merged_seurat@meta.data$Clusters <- ifelse(is.na(merged_seurat@meta.data$Clusters), "GSM7919060", merged_seurat@meta.data$Clusters)
merged_seurat@meta.data$condition <- ifelse(is.na(merged_seurat@meta.data$condition), "GSM7919060", merged_seurat@meta.data$condition)

# Save the integrated object
saveRDS(merged_seurat, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/eos_GSM7919060_harmony_integrated.rds")

# Print a summary of the object
print(merged_seurat)

# Print the number of cells from each original assay
table(merged_seurat$orig.ident)
table(merged_seurat$Clusters)

#Assay2     KO     WT
# 1086   6166   4532

# GSM7919060 The original dataset was quite large, containing 22,597 features (genes) across 6,794,880 cells. However, after applying the filtering criteria, there was a dramatic reduction in both cells and features. Only 1,086 cells met the criterion of having at least 200 features detected, which is a mere 0.016% of the original cell count. This suggests that the vast majority of the original cells had very low gene detection rates, possibly indicating a high proportion of empty droplets or cells with poor RNA capture. On the feature side, 12,455 features were detected in at least 3 cells, representing about 55% of the original features. The final count of cells meeting both criteria (at least 200 features per cell and each feature detected in at least 3 cells) remained at 1,086. 

