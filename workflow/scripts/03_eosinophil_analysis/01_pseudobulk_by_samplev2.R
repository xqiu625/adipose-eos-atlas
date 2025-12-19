library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(tibble)


# Load the seurat object
seurat_object <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/seurat_object_for_cellchat.rds")

out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'
dir.create(out_dir, recursive = TRUE)

# ===============================================
# PSEUDOBULK FUNCTION FOR PCA (averaged log-normalized expression)
# ===============================================

create_pseudobulk_for_pca <- function(seurat_obj, group_by = "orig.ident", assay = "RNA", slot = "data") {
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  meta <- seurat_obj@meta.data
  meta$group <- meta[[group_by]]
  
  pseudobulk_list <- list()
  for (group in unique(meta$group)) {
    cells <- rownames(meta)[meta$group == group]
    if (length(cells) > 1) {
      pseudobulk_list[[group]] <- Matrix::rowMeans(expr[, cells])
    } else {
      pseudobulk_list[[group]] <- expr[, cells]
    }
  }
  
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  return(pseudobulk_matrix)
}

# ===============================================
# CREATE PSEUDOBULK & METADATA
# ===============================================

# Create pseudobulk matrix (averaged log-normalized expression)
cat("Creating pseudobulk matrix for PCA...\n")
pseudobulk_matrix <- create_pseudobulk_for_pca(seurat_object, group_by = "orig.ident")

# Filter genes with zero variance
gene_vars <- apply(pseudobulk_matrix, 1, var)
pseudobulk_filtered <- pseudobulk_matrix[gene_vars > 0, ]
cat("Removed", sum(gene_vars == 0), "zero-variance genes\n")

# Create sample metadata
cat("Creating sample metadata...\n")
sample_metadata <- seurat_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    genotype = first(genotype),
    Sex = first(Sex),
    condition = first(condition),
    infection_status = first(infection_status),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  column_to_rownames("orig.ident")

# Align metadata with pseudobulk columns
sample_metadata <- sample_metadata[colnames(pseudobulk_filtered), ]

# Save pseudobulk matrix (filtered for PCA)
saveRDS(pseudobulk_filtered, file.path(out_dir, "pseudobulk_all_celltypes.rds"))

# Save sample metadata
saveRDS(sample_metadata, file.path(out_dir, "sample_metadata_all.rds"))



# Subset Seurat object to only include Macrophages
mac_seurat <- subset(seurat_object, subset = cell_type == "Macrophages")

# Create pseudobulk matrix (averaged log-normalized expression)
cat("Creating pseudobulk matrix for PCA...\n")
pseudobulk_mac <- create_pseudobulk_for_pca(mac_seurat, group_by = "orig.ident")
gene_vars_mac <- apply(pseudobulk_mac, 1, var)
pseudobulk_mac_filtered <- pseudobulk_mac[gene_vars_mac > 0, ]
# Create sample metadata for macrophages
cat("Creating sample metadata for Macrophages...\n")
sample_metadata_mac <- mac_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    genotype = first(genotype),
    Sex = first(Sex),
    condition = first(condition),
    infection_status = first(infection_status),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  column_to_rownames("orig.ident")

# Align metadata with pseudobulk columns for macrophages
sample_metadata_mac <- sample_metadata_mac[colnames(pseudobulk_mac_filtered), ]
# Save macrophage pseudobulk and metadata
saveRDS(pseudobulk_mac_filtered, file.path(out_dir, "pseudobulk_macrophages.rds"))
saveRDS(sample_metadata_mac, file.path(out_dir, "sample_metadata_macrophages.rds"))

# ===============================================
# Subset Seurat object to only include Eosinophils
eos_seurat <- subset(seurat_object, subset = cell_type == "Eosinophils")

# Create pseudobulk matrix (averaged log-normalized expression)
cat("Creating pseudobulk matrix for PCA...\n")
pseudobulk_eos <- create_pseudobulk_for_pca(eos_seurat, group_by = "orig.ident")

# Filter genes with zero variance for eosinophils
gene_vars_eos <- apply(pseudobulk_eos, 1, var)
pseudobulk_eos_filtered <- pseudobulk_eos[gene_vars_eos > 0, ]
cat("Eosinophils: Removed", sum(gene_vars_eos == 0), "zero-variance genes\n")

# Create sample metadata for eosinophils
cat("Creating sample metadata for Eosinophils...\n")
sample_metadata_eos <- eos_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    genotype = first(genotype),
    Sex = first(Sex),
    condition = first(condition),
    infection_status = first(infection_status),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  column_to_rownames("orig.ident")

# Align metadata with pseudobulk columns for eosinophils
sample_metadata_eos <- sample_metadata_eos[colnames(pseudobulk_eos_filtered), ]

# Save eosinophil pseudobulk and metadata
saveRDS(pseudobulk_eos_filtered, file.path(out_dir, "pseudobulk_eosinophils.rds"))
saveRDS(sample_metadata_eos, file.path(out_dir, "sample_metadata_eosinophils.rds"))

# ===============================================
# Subset Seurat object to only include Eosinophils
naive_seurat <- subset(seurat_object, subset = infection_status == "Not_Infected")
