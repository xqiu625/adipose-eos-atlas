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

# Helper function to create pseudobulk and metadata
create_pseudobulk_and_metadata <- function(seurat_obj, output_prefix) {
  cat("Creating pseudobulk matrix for", output_prefix, "...\n")
  
  # Create pseudobulk matrix
  pseudobulk_matrix <- create_pseudobulk_for_pca(seurat_obj, group_by = "orig.ident")
  
  # Filter genes with zero variance
  gene_vars <- apply(pseudobulk_matrix, 1, var)
  pseudobulk_filtered <- pseudobulk_matrix[gene_vars > 0, ]
  cat(output_prefix, ": Removed", sum(gene_vars == 0), "zero-variance genes\n")
  
  # Create sample metadata
  sample_metadata <- seurat_obj@meta.data %>%
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
  
  # Save files
  saveRDS(pseudobulk_filtered, file.path(out_dir, paste0("pseudobulk_", output_prefix, ".rds")))
  saveRDS(sample_metadata, file.path(out_dir, paste0("sample_metadata_", output_prefix, ".rds")))
  
  cat(output_prefix, "pseudobulk matrix dimensions:", dim(pseudobulk_filtered), "\n")
  cat(output_prefix, "samples:", ncol(pseudobulk_filtered), "\n\n")
  
  return(list(pseudobulk = pseudobulk_filtered, metadata = sample_metadata))
}

# ===============================================
# 1. ALL CELLS - THREE CONDITIONS
# ===============================================

# 1a. All cells - All samples (16 samples)
all_cells_all <- create_pseudobulk_and_metadata(seurat_object, "all_cells_all_samples")

# 1b. All cells - Not infected samples only
seurat_not_infected <- subset(seurat_object, subset = infection_status == "Not_Infected")
all_cells_not_infected <- create_pseudobulk_and_metadata(seurat_not_infected, "all_cells_not_infected")

# 1c. All cells - Infected samples only
seurat_infected <- subset(seurat_object, subset = infection_status == "Infected")
all_cells_infected <- create_pseudobulk_and_metadata(seurat_infected, "all_cells_infected")

# ===============================================
# 2. EOSINOPHILS - THREE CONDITIONS
# ===============================================

# Subset to eosinophils
eos_seurat <- subset(seurat_object, subset = cell_type == "Eosinophils")

# 2a. Eosinophils - All samples
eos_all <- create_pseudobulk_and_metadata(eos_seurat, "eosinophils_all_samples")

# 2b. Eosinophils - Not infected samples only
eos_seurat_not_infected <- subset(eos_seurat, subset = infection_status == "Not_Infected")
eos_not_infected <- create_pseudobulk_and_metadata(eos_seurat_not_infected, "eosinophils_not_infected")

# 2c. Eosinophils - Infected samples only
eos_seurat_infected <- subset(eos_seurat, subset = infection_status == "Infected")
eos_infected <- create_pseudobulk_and_metadata(eos_seurat_infected, "eosinophils_infected")

# ===============================================
# 3. MACROPHAGES - THREE CONDITIONS
# ===============================================

# Subset to macrophages
mac_seurat <- subset(seurat_object, subset = cell_type == "Macrophages")

# 3a. Macrophages - All samples
mac_all <- create_pseudobulk_and_metadata(mac_seurat, "macrophages_all_samples")

# 3b. Macrophages - Not infected samples only
mac_seurat_not_infected <- subset(mac_seurat, subset = infection_status == "Not_Infected")
mac_not_infected <- create_pseudobulk_and_metadata(mac_seurat_not_infected, "macrophages_not_infected")

# 3c. Macrophages - Infected samples only
mac_seurat_infected <- subset(mac_seurat, subset = infection_status == "Infected")
mac_infected <- create_pseudobulk_and_metadata(mac_seurat_infected, "macrophages_infected")

# ===============================================
# SUMMARY
# ===============================================

cat("\n=== COMPREHENSIVE PSEUDOBULK SUMMARY ===\n")
cat("Output directory:", out_dir, "\n\n")

cat("1. ALL CELLS:\n")
cat("   - All samples (16):", dim(all_cells_all$pseudobulk), "\n")
cat("   - Not infected:", dim(all_cells_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(all_cells_infected$pseudobulk), "\n\n")

cat("2. EOSINOPHILS:\n")
cat("   - All samples:", dim(eos_all$pseudobulk), "\n")
cat("   - Not infected:", dim(eos_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(eos_infected$pseudobulk), "\n\n")

cat("3. MACROPHAGES:\n")
cat("   - All samples:", dim(mac_all$pseudobulk), "\n")
cat("   - Not infected:", dim(mac_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(mac_infected$pseudobulk), "\n\n")

cat("FILES CREATED:\n")
files_created <- c(
  "pseudobulk_all_cells_all_samples.rds",
  "pseudobulk_all_cells_not_infected.rds", 
  "pseudobulk_all_cells_infected.rds",
  "pseudobulk_eosinophils_all_samples.rds",
  "pseudobulk_eosinophils_not_infected.rds",
  "pseudobulk_eosinophils_infected.rds",
  "pseudobulk_macrophages_all_samples.rds",
  "pseudobulk_macrophages_not_infected.rds",
  "pseudobulk_macrophages_infected.rds",
  "sample_metadata_all_cells_all_samples.rds",
  "sample_metadata_all_cells_not_infected.rds",
  "sample_metadata_all_cells_infected.rds",
  "sample_metadata_eosinophils_all_samples.rds",
  "sample_metadata_eosinophils_not_infected.rds",
  "sample_metadata_eosinophils_infected.rds",
  "sample_metadata_macrophages_all_samples.rds",
  "sample_metadata_macrophages_not_infected.rds",
  "sample_metadata_macrophages_infected.rds"
)

for (file in files_created) {
  cat("   -", file, "\n")
}

cat("\nTotal files created:", length(files_created), "\n")