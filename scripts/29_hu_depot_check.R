#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 29_hu_depot_check.R
# Quick inspection: does the Hu dataset have depot/tissue metadata?
# Output: prints metadata columns and unique values to stdout
# SLURM:  batch, 32G, 1 CPU, 00:30:00
# ==============================================================================

library(Seurat)

hu_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/integration/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/hu_eos/hu_adipose_eos.rds"
)

# Also check for full Hu object (not just eos subset)
hu_full_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_full.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/hu_eos/hu_full.rds"
)

# Try eos subset first
hu <- NULL
for (path in c(hu_candidates, hu_full_candidates)) {
  if (file.exists(path)) {
    cat("Found:", path, "\n")
    hu <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(hu)) break
  }
}

# Also check local data directory
local_hu <- list.files("/bigdata/nairlab/shared/XQiu/EOS_project",
                       pattern = "hu|Hu", recursive = TRUE, full.names = TRUE)
local_hu <- local_hu[grepl("\\.rds$", local_hu)]
cat("\nAll Hu-related RDS files found:\n")
for (f in local_hu) cat("  ", f, "\n")

if (is.null(hu)) {
  cat("\nERROR: Could not load any Hu object\n")
  quit(status = 1)
}

cat(sprintf("\nLoaded: %d cells x %d genes\n", ncol(hu), nrow(hu)))
cat("\n=== Metadata columns ===\n")
cat(paste(colnames(hu@meta.data), collapse = "\n"), "\n")

cat("\n=== Unique values per metadata column ===\n")
for (col in colnames(hu@meta.data)) {
  vals <- unique(hu@meta.data[[col]])
  if (length(vals) <= 30) {
    cat(sprintf("\n%s (%d unique):\n", col, length(vals)))
    cat("  ", paste(head(sort(as.character(vals)), 20), collapse = ", "), "\n")
  } else {
    cat(sprintf("\n%s (%d unique): [too many to show]\n", col, length(vals)))
  }
}

# Look specifically for depot/tissue-related columns
depot_keywords <- c("depot", "tissue", "organ", "site", "location", "fat", "bat", "iwat", "ewat", "gwat")
matches <- grep(paste(depot_keywords, collapse = "|"), colnames(hu@meta.data), ignore.case = TRUE, value = TRUE)
if (length(matches) > 0) {
  cat("\n=== DEPOT-RELATED COLUMNS FOUND ===\n")
  for (col in matches) {
    cat(sprintf("\n%s:\n", col))
    print(table(hu@meta.data[[col]]))
  }
} else {
  cat("\n=== NO DEPOT-RELATED COLUMNS FOUND ===\n")
  cat("Hu data may not have depot metadata, or it may use a different column name.\n")
}
