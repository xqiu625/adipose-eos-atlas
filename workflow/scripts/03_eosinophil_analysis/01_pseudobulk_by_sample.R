# Pseudobulk Analysis for Single-Cell RNA-seq Data
# This script aggregates single-cell data to pseudobulk for proper statistical analysis

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(tibble)

# Load the seurat object
seurat_object <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/seurat_object_for_cellchat.rds")

out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/04_PCA_MDSv2'
dir.create(out_dir, recursive = TRUE)

# Subset Seurat object to only include Macrophages
mac_seurat <- subset(seurat_object, subset = cell_type == "Macrophages")
eos_seurat <- subset(seurat_object, subset = cell_type == "Eosinophils")





# ===============================================
# PSEUDOBULK FUNCTION
# ===============================================

create_pseudobulk <- function(seurat_obj, group_by = "orig.ident") {
  
  # Get raw counts
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # Get metadata
  meta <- seurat_obj@meta.data
  meta$group <- meta[[group_by]]
  
  # Create pseudobulk by summing counts within each group
  pseudobulk_list <- list()
  
  for(group in unique(meta$group)) {
    cells_in_group <- rownames(meta)[meta$group == group]
    if(length(cells_in_group) > 1) {
      pseudobulk_list[[group]] <- Matrix::rowSums(counts[, cells_in_group])
    } else {
      pseudobulk_list[[group]] <- counts[, cells_in_group]
    }
  }
  
  # Combine into matrix
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  
  return(pseudobulk_matrix)
}

# ===============================================
# CREATE PSEUDOBULK FOR ALL CELL TYPES
# ===============================================

print("Creating overall pseudobulk...")
pseudobulk_all <- create_pseudobulk(seurat_object, group_by = "orig.ident")

# Create sample metadata
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

sample_metadata <- sample_metadata[colnames(pseudobulk_all), ]

print(paste("Overall pseudobulk:", nrow(pseudobulk_all), "genes x", ncol(pseudobulk_all), "samples"))

# ===============================================
# CREATE CELL TYPE-SPECIFIC PSEUDOBULK
# ===============================================

# Macrophages
print("Creating macrophage pseudobulk...")
mac_seurat <- subset(seurat_object, subset = cell_type == "Macrophages")
pseudobulk_mac <- create_pseudobulk(mac_seurat, group_by = "orig.ident")

mac_metadata <- mac_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    genotype = first(genotype),
    Sex = first(Sex),
    infection_status = first(infection_status),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  column_to_rownames("orig.ident")

mac_metadata <- mac_metadata[colnames(pseudobulk_mac), ]

print(paste("Macrophage pseudobulk:", nrow(pseudobulk_mac), "genes x", ncol(pseudobulk_mac), "samples"))

# Eosinophils
print("Creating eosinophil pseudobulk...")
eos_seurat <- subset(seurat_object, subset = cell_type == "Eosinophils")
pseudobulk_eos <- create_pseudobulk(eos_seurat, group_by = "orig.ident")

eos_metadata <- eos_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    genotype = first(genotype),
    Sex = first(Sex),
    infection_status = first(infection_status),
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  column_to_rownames("orig.ident")

eos_metadata <- eos_metadata[colnames(pseudobulk_eos), ]

print(paste("Eosinophil pseudobulk:", nrow(pseudobulk_eos), "genes x", ncol(pseudobulk_eos), "samples"))

# ===============================================
# SAVE RESULTS
# ===============================================

# Save all pseudobulk data
saveRDS(pseudobulk_all, file.path(out_dir, "pseudobulk_all_celltypes.rds"))
saveRDS(sample_metadata, file.path(out_dir, "sample_metadata_all.rds"))

saveRDS(pseudobulk_mac, file.path(out_dir, "pseudobulk_macrophages.rds"))
saveRDS(mac_metadata, file.path(out_dir, "sample_metadata_macrophages.rds"))

saveRDS(pseudobulk_eos, file.path(out_dir, "pseudobulk_eosinophils.rds"))
saveRDS(eos_metadata, file.path(out_dir, "sample_metadata_eosinophils.rds"))

print("Pseudobulk analysis complete!")
print(paste("Results saved in:", out_dir))
print("Files created:")
print("- pseudobulk_all_celltypes.rds")
print("- sample_metadata_all.rds")
print("- pseudobulk_macrophages.rds") 
print("- sample_metadata_macrophages.rds")
print("- pseudobulk_eosinophils.rds")
print("- sample_metadata_eosinophils.rds")