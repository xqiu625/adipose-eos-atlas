# Load required packages
library(Seurat)
library(dplyr)

# ================================
# Define M1 and M2 macrophage signatures
# ================================

# Source: eLife 2023;12:e86001 — Supplemental Table 2
# M1 and M2 markers were curated from experimental perturbation models and literature

# M1 markers (pro-inflammatory macrophage genes)
M1_markers <- c(
  "Il1b", "Tnf", "Nos2", "Cxcl9", "Cxcl10", "Il6",
  "Ccl2", "Ccl5", "Ccl3", "Ccl4", "Stat1", "Irf5",
  "Cd86", "Cd80", "Socs3"
)

# M2 markers (anti-inflammatory / tissue-repair macrophage genes)
M2_markers <- c(
  "Arg1", "Mrc1", "Chi3l3", "Retnla", "Il10",
  "Trem2", "Cd163", "Cd206", "Ccl17", "Ccl22",
  "Pparg", "Mgl2", "Stat6", "Socs1"
)

# Save gene lists to file for traceability
dir.create("/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/resources", showWarnings = FALSE, recursive = TRUE)

writeLines(M1_markers, "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/resources/M1_gene_list.txt")
writeLines(M2_markers, "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/04_M1_M2_Scoring/resources/M2_gene_list.txt")

cat("M1 and M2 gene signatures saved successfully.\n")