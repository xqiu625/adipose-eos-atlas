lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# Now load the libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)

# local Eosinophils 
eos <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")

# GSM7919060 (Cell immunity paper Eosinophils)
data_dir <- "/bigdata/godziklab/shared/Xinru/302025/dataset/Immunity_24May/GSM7919060"
GSM7919060_seurat_object <- Read10X(data.dir = data_dir)
GSM7919060_seurat_object <- CreateSeuratObject(counts = GSM7919060_seurat_object, min.cells = 3, min.features = 200)

# Posner et al. (2022) Current Opinion in Immunology.
posner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/95e0ceb5-45fc-489b-ae13-c2a7e4a6dfa5.rds")
head(posner@meta.data)
table(posner$cell_type)
# Extract eosinophils from Posner et al. dataset
eos_posner <- subset(posner, subset = cell_type == "eosinophil")

# Print some information about the extracted eosinophils
print(paste("Number of eosinophils extracted from Posner et al. dataset:", ncol(eos_posner)))

# Add a dataset identifier to the metadata
eos_posner$dataset <- "Posner_et_al"

# Save the eos_posner Seurat object as an RDS file
output_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/"
saveRDS(eos_posner, file = file.path(output_dir, "eos_posner.rds"))


# GSE182001 (Gurtner Nature paper Eosinophils)
# Load required libraries
library(org.Mm.eg.db)
library(AnnotationDbi)

GSE182001 <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_merged_tissues_seurat_object.rds")
# Get the first few gene names
gene_names <- head(row.names(GSE182001))

# Function to convert gene names
convert_gene_names <- function(genes) {
  # Get NCBI gene IDs
  gene_ids <- mapIds(org.Mm.eg.db,
                     keys = genes,
                     column = "ENTREZID",
                     keytype = "ALIAS",
                     multiVals = "first")
  
  # Get NCBI gene symbols
  gene_symbols <- mapIds(org.Mm.eg.db,
                         keys = gene_ids,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  return(gene_symbols)
}


# Convert all gene names
all_genes <- rownames(GSE182001)
all_ncbi_names <- convert_gene_names(all_genes)

# Create unique row names
make_unique_names <- function(original, ncbi) {
  ifelse(is.na(ncbi), original, 
         ifelse(duplicated(ncbi) | duplicated(ncbi, fromLast = TRUE), 
                paste0(ncbi, "_", original), 
                ncbi))
}

new_rownames <- make_unique_names(all_genes, all_ncbi_names)
rownames(GSE182001) <- new_rownames

# Identify eosinophil marker genes
eos_markers <- c("Siglecf", "Il5ra", "Ccr3", "Epx")
found_markers <- eos_markers[eos_markers %in% rownames(GSE182001)]

if (length(found_markers) < length(eos_markers)) {
  missing_markers <- setdiff(eos_markers, found_markers)
  warning(paste("Some marker genes were not found in the dataset:", 
                paste(missing_markers, collapse = ", ")))
}

# Normalize data if not already normalized
# Perform clustering
GSE182001 <- NormalizeData(GSE182001)
GSE182001 <- FindVariableFeatures(GSE182001, selection.method = "vst", nfeatures = 2000)
GSE182001 <- ScaleData(GSE182001)
GSE182001 <- RunPCA(GSE182001)
GSE182001 <- FindNeighbors(GSE182001, dims = 1:30)
GSE182001 <- FindClusters(GSE182001, resolution = 0.8)
GSE182001 <- RunUMAP(GSE182001, dims = 1:30)


# Check which clusters express eosinophil markers
eos_cluster_expression <- AverageExpression(GSE182001, features = found_markers, group.by = "seurat_clusters")
eos_cluster_expression <- eos_cluster_expression$RNA

# Print eosinophil marker expression in each cluster
print("Eosinophil marker expression in each cluster:")
print(eos_cluster_expression)

# Identify potential eosinophil clusters
eos_score <- colMeans(eos_cluster_expression)
potential_eos_clusters <- which(eos_score > quantile(eos_score, 0.75))

# Subset potential eosinophils
eosinophils <- subset(GSE182001, idents = potential_eos_clusters)
eosinophils
eosinophils$dataset <- "Gurtner_et_al"

# Save the eosinophil subset
saveRDS(eosinophils, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_eosinophil_subset.rds")

# Print a message to confirm the save
print(paste("Eosinophil subset saved with", ncol(eosinophils), "cells"))
