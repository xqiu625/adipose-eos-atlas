# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)

# ======================
# Load data and settings
# ======================

# Input file: Seurat object
input_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed/seurat_object_labeled.rds"
seurat_obj <- readRDS(input_file)

# Load M1 and M2 gene sets
m1_genes <- readLines("/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/resources/M1_gene_list.txt")
m2_genes <- readLines("/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/resources/M2_gene_list.txt")

# Remove genes not present in dataset
m1_genes <- intersect(m1_genes, rownames(seurat_obj))
m2_genes <- intersect(m2_genes, rownames(seurat_obj))

# ======================
# Score cells
# ======================

# Add module scores for M1 and M2
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(M1 = m1_genes),
  name = "M1_Score"
)

seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(M2 = m2_genes),
  name = "M2_Score"
)

# ======================
# Save scored object
# ======================

output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/seurat_object_M1_M2_scored.rds"
saveRDS(seurat_obj, output_file)

