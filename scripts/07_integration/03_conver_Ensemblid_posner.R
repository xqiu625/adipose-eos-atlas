# Load required libraries
library(Seurat)
library(harmony)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Load the datasets
eos_this_study <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")
eos_gurtner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/GSE182001_eosinophil_subset.rds")
eos_posner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/eos_posner.rds")

data_dir <- "/bigdata/godziklab/shared/Xinru/302025/dataset/Immunity_24May/GSM7919060"
GSM7919060_seurat_object <- Read10X(data.dir = data_dir)
eos_jorssen <- CreateSeuratObject(counts = GSM7919060_seurat_object, min.cells = 3, min.features = 200)

print(eos_this_study)
print(eos_gurtner) # BD Rhapsody Single-Cell Analysis System
print(eos_posner)
print(eos_jorssen)

# Standardizing Gene Names for Mouse Single-Cell Data
# Function to summarize gene name conversion
summarize_conversion <- function(seurat_obj, name) {
  total_genes <- nrow(seurat_obj)
  ensembl_count <- sum(grepl("^ENSMUS", rownames(seurat_obj)))
  symbol_count <- total_genes - ensembl_count
  
  cat(name, "summary:\n")
  cat("Total genes:", total_genes, "\n")
  cat("Genes converted to symbols:", symbol_count, "\n")
  cat("Genes remaining as Ensembl IDs:", ensembl_count, "\n")
  cat("First few gene names:", head(rownames(seurat_obj)), "\n\n")
}

# Check the results for each dataset
summarize_conversion(eos_this_study, "eos_this_study")
summarize_conversion(eos_gurtner, "eos_gurtner")
summarize_conversion(eos_posner, "eos_posner")
summarize_conversion(eos_jorssen, "eos_jorssen")

# > summarize_conversion(eos_posner, "eos_posner")
# eos_posner summary:
# Total genes: 25232 
# Genes converted to symbols: 0 
# Genes remaining as Ensembl IDs: 25232 
# First few gene names: ENSMUSG00000051951 ENSMUSG00000089699 ENSMUSG00000102343 ENSMUSG00000025900 ENSMUSG00000025902 ENSMUSG00000104328 

# Update eos_posner Seurat object with gene symbols


# Load required libraries
library(Seurat)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)

eos_posner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/eos_posner.rds")


# Step 1: Extract the necessary data
counts_matrix <- GetAssayData(eos_posner, slot = "counts")
metadata <- eos_posner@meta.data
feature_metadata <- eos_posner[["RNA"]][[]]

# Step 2: Get gene symbols for Ensembl IDs
ensembl_ids <- rownames(counts_matrix)
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Load required libraries
library(Seurat)
library(dplyr)

# Assuming eos_posner is your current Seurat object

# Step 1: Extract the necessary data
counts_matrix <- GetAssayData(eos_posner, slot = "counts")
metadata <- eos_posner@meta.data
feature_metadata <- eos_posner[["RNA"]][[]]

# Step 2: Create new feature names from feature_name column, handling duplicates
new_feature_names <- make.unique(as.character(feature_metadata$feature_name))

# Step 3: Rename the rows of the counts matrix
rownames(counts_matrix) <- new_feature_names

# Step 4: Create a new Seurat object
eos_posner_new <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)

# Step 5: Update and add feature metadata
feature_metadata_new <- feature_metadata
rownames(feature_metadata_new) <- new_feature_names
eos_posner_new[["RNA"]] <- AddMetaData(eos_posner_new[["RNA"]], metadata = feature_metadata_new)

# Step 6: Copy over dimensional reductions (if any)
for (reduction in names(eos_posner@reductions)) {
  reduction_data <- Embeddings(eos_posner, reduction)
  eos_posner_new[[reduction]] <- CreateDimReducObject(embeddings = reduction_data, key = Key(eos_posner[[reduction]]), assay = DefaultAssay(eos_posner))
}

# Step 7: Recreate the scaled data and variable features (if present in the original object)
if ("scale.data" %in% slotNames(eos_posner[["RNA"]])) {
  eos_posner_new <- ScaleData(eos_posner_new, features = rownames(eos_posner_new))
}

if (length(VariableFeatures(eos_posner)) > 0) {
  var_features <- VariableFeatures(eos_posner)
  var_features_new <- new_feature_names[match(var_features, feature_metadata$ensembl_ids)]
  VariableFeatures(eos_posner_new) <- var_features_new
}

# Step 8: Verify the new object
print("Dimensions of the new RNA assay:")
print(dim(eos_posner_new[["RNA"]]))

print("New feature names (from feature_name):")
print(head(rownames(eos_posner_new[["RNA"]]), n = 10))

print("Updated feature metadata:")
print(head(eos_posner_new[["RNA"]][[]], n = 10))

# If everything looks good, you can replace the original object
eos_posner <- eos_posner_new

# Optional: Remove the temporary object to free up memory
rm(eos_posner_new)
gc()

# Final verification
print("Final row names:")
print(head(rownames(eos_posner[["RNA"]]), n = 10))

saveRDS(eos_posner, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/eos_posner_updated_geneSymbols.rds")
