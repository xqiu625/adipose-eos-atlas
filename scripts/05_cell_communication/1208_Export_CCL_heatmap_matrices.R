# Export CCL Heatmap Matrices - FULL UNFILTERED DATA
# Purpose: Export the FULL unfiltered matrices for CCL_Macrophage_to_Eosinophil
#          and CCL_Eosinophil_to_Macrophage (no pathway or threshold filtering)
# Date: 2024-12-08

library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)

# =============================================================================
# SETUP DIRECTORIES AND LOAD EXISTING DATA
# =============================================================================

out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/CCL_eos_subclusters_updated_mac'
export_dir <- file.path(out_dir, 'exported_matrices')
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

print("=== LOADING EXISTING CCL DATA ===")

# Load the bidirectional eos-mac CCL data
eos_mac_ccl <- read.csv(file.path(out_dir, "Eos_Mac_CCL_Bidirectional.csv"))

print(paste("Loaded", nrow(eos_mac_ccl), "rows of Eos-Mac CCL data"))
print("Direction summary:")
print(table(eos_mac_ccl$Direction))

# =============================================================================
# 1. CCL_Eosinophil_to_Macrophage - FULL UNFILTERED EXPORT
# =============================================================================

print("\n=== EXPORTING CCL_Eosinophil_to_Macrophage FULL MATRIX ===")

eos_to_mac_ccl <- eos_mac_ccl[eos_mac_ccl$Direction == "Eos_to_Mac", ]

# Export raw data (long format) - ALL data including zero strength
write.csv(eos_to_mac_ccl,
          file.path(export_dir, "CCL_Eosinophil_to_Macrophage_FULL_raw.csv"),
          row.names = FALSE)
print(paste("Raw data rows:", nrow(eos_to_mac_ccl)))

# Create summary (mean across conditions) - NO FILTERING
eos_to_mac_summary <- eos_to_mac_ccl %>%
  group_by(Source, Target, Pathway) %>%
  summarise(
    Mean_Strength = mean(Strength),
    Max_Strength = max(Strength),
    Min_Strength = min(Strength),
    N_Conditions = n(),
    .groups = 'drop'
  )

eos_to_mac_summary$Interaction_Label <- paste(
  eos_to_mac_summary$Target,
  eos_to_mac_summary$Pathway,
  sep = "."
)

# Export full summary data (long format)
write.csv(eos_to_mac_summary,
          file.path(export_dir, "CCL_Eosinophil_to_Macrophage_FULL_summary.csv"),
          row.names = FALSE)

# Create FULL heatmap matrix (no filtering)
eos_to_mac_matrix <- acast(eos_to_mac_summary,
                        Source ~ Interaction_Label,
                        value.var = "Mean_Strength",
                        fill = 0)

# Export FULL matrix (wide format)
write.csv(eos_to_mac_matrix,
          file.path(export_dir, "CCL_Eosinophil_to_Macrophage_FULL_matrix.csv"),
          row.names = TRUE)

print(paste("Full matrix dimensions:", nrow(eos_to_mac_matrix), "rows x", ncol(eos_to_mac_matrix), "cols"))
print("Exported: CCL_Eosinophil_to_Macrophage_FULL_raw.csv")
print("Exported: CCL_Eosinophil_to_Macrophage_FULL_summary.csv")
print("Exported: CCL_Eosinophil_to_Macrophage_FULL_matrix.csv")

# =============================================================================
# 2. CCL_Macrophage_to_Eosinophil - FULL UNFILTERED EXPORT
# =============================================================================

print("\n=== EXPORTING CCL_Macrophage_to_Eosinophil FULL MATRIX ===")

mac_to_eos_ccl <- eos_mac_ccl[eos_mac_ccl$Direction == "Mac_to_Eos", ]

# Export raw data (long format) - ALL data including zero strength
write.csv(mac_to_eos_ccl,
          file.path(export_dir, "CCL_Macrophage_to_Eosinophil_FULL_raw.csv"),
          row.names = FALSE)
print(paste("Raw data rows:", nrow(mac_to_eos_ccl)))

# Create summary (mean across conditions) - NO FILTERING
mac_to_eos_summary <- mac_to_eos_ccl %>%
  group_by(Source, Target, Pathway) %>%
  summarise(
    Mean_Strength = mean(Strength),
    Max_Strength = max(Strength),
    Min_Strength = min(Strength),
    N_Conditions = n(),
    .groups = 'drop'
  )

mac_to_eos_summary$Interaction_Label <- paste(
  mac_to_eos_summary$Target,
  mac_to_eos_summary$Pathway,
  sep = "."
)

# Export full summary data (long format)
write.csv(mac_to_eos_summary,
          file.path(export_dir, "CCL_Macrophage_to_Eosinophil_FULL_summary.csv"),
          row.names = FALSE)

# Create FULL heatmap matrix (no filtering)
mac_to_eos_matrix <- acast(mac_to_eos_summary,
                        Source ~ Interaction_Label,
                        value.var = "Mean_Strength",
                        fill = 0)

# Export FULL matrix (wide format)
write.csv(mac_to_eos_matrix,
          file.path(export_dir, "CCL_Macrophage_to_Eosinophil_FULL_matrix.csv"),
          row.names = TRUE)

print(paste("Full matrix dimensions:", nrow(mac_to_eos_matrix), "rows x", ncol(mac_to_eos_matrix), "cols"))
print("Exported: CCL_Macrophage_to_Eosinophil_FULL_raw.csv")
print("Exported: CCL_Macrophage_to_Eosinophil_FULL_summary.csv")
print("Exported: CCL_Macrophage_to_Eosinophil_FULL_matrix.csv")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n=== EXPORT COMPLETE ===")
print(paste("Output directory:", export_dir))
print("\nExported files (FULL UNFILTERED DATA):")
print("  Eosinophil to Macrophage:")
print("    1. CCL_Eosinophil_to_Macrophage_FULL_raw.csv     - All raw data per condition")
print("    2. CCL_Eosinophil_to_Macrophage_FULL_summary.csv - Summary stats (long format)")
print("    3. CCL_Eosinophil_to_Macrophage_FULL_matrix.csv  - Heatmap matrix (wide format)")
print("  Macrophage to Eosinophil:")
print("    4. CCL_Macrophage_to_Eosinophil_FULL_raw.csv     - All raw data per condition")
print("    5. CCL_Macrophage_to_Eosinophil_FULL_summary.csv - Summary stats (long format)")
print("    6. CCL_Macrophage_to_Eosinophil_FULL_matrix.csv  - Heatmap matrix (wide format)")

print("\n=== DATA COLUMN DESCRIPTIONS ===")
print("Raw files: Source, Target, Pathway, Strength, Condition, Direction")
print("Summary files: Source, Target, Pathway, Mean_Strength, Max_Strength, Min_Strength, N_Conditions, Interaction_Label")
print("Matrix files: Rows = Source cell types, Columns = Target.Pathway interactions")
