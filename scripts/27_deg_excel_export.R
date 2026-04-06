#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 27_deg_excel_export.R
# Bundle cluster markers + all DEG CSVs into a single Excel workbook
# Input:  $BASE/eos_cluster_markers.csv, $BASE/deg/*.csv
# Output: $BASE/deg/eos_deg_results.xlsx
# SLURM:  batch, 16G, 1 CPU, 00:30:00
# ==============================================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org",
                   lib = lib_path)
}
library(openxlsx)

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
deg_dir <- file.path(BASE, "deg")
out_file <- file.path(deg_dir, "eos_deg_results.xlsx")

wb <- createWorkbook()

# --- Sheet 1: Cluster markers (FindAllMarkers, from 04_eos_subcluster.R) ------
markers_file <- file.path(BASE, "eos_cluster_markers.csv")
if (file.exists(markers_file)) {
  markers <- read.csv(markers_file)
  addWorksheet(wb, "Cluster_Markers")
  writeData(wb, "Cluster_Markers", markers)
  cat("  Added Cluster_Markers:", nrow(markers), "genes\n")
} else {
  cat("  WARNING: eos_cluster_markers.csv not found\n")
}

# --- Remaining sheets: DEG CSVs from /deg/ ------------------------------------
csv_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
# Sort for consistent ordering
csv_files <- sort(csv_files)

for (f in csv_files) {
  fname <- tools::file_path_sans_ext(basename(f))
  # Truncate sheet name to 31 chars (Excel limit)
  sheet_name <- substr(gsub("^deg_", "", fname), 1, 31)
  df <- read.csv(f)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  cat(sprintf("  Added %s: %d rows\n", sheet_name, nrow(df)))
}

saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\nSaved:", out_file, "\n")
cat("Sheets:", paste(names(wb), collapse = ", "), "\n")
