#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

# ==============================================================================
# 29_hu_depot_check.R
# Quick inspection: does the Hu dataset have depot/tissue metadata?
# SLURM:  batch, 64G, 1 CPU, 00:30:00
# ==============================================================================

library(Seurat)

path <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds"
cat("Loading:", path, "\n")
hu <- readRDS(path)
if (inherits(hu[["RNA"]], "Assay5")) hu[["RNA"]] <- JoinLayers(hu[["RNA"]])

cat(sprintf("Loaded: %d cells x %d genes\n", ncol(hu), nrow(hu)))

cat("\n=== Metadata columns ===\n")
for (col in colnames(hu@meta.data)) {
  vals <- unique(hu@meta.data[[col]])
  if (length(vals) <= 30) {
    cat(sprintf("  %s (%d): %s\n", col, length(vals),
                paste(head(sort(as.character(vals)), 15), collapse = ", ")))
  } else {
    cat(sprintf("  %s (%d unique)\n", col, length(vals)))
  }
}

# Search for depot-related columns
depot_kw <- c("depot", "tissue", "organ", "site", "fat", "bat", "iwat", "ewat", "gwat", "location")
matches <- grep(paste(depot_kw, collapse = "|"), colnames(hu@meta.data), ignore.case = TRUE, value = TRUE)
if (length(matches) > 0) {
  cat("\n=== DEPOT COLUMNS FOUND ===\n")
  for (col in matches) {
    cat(sprintf("\n%s:\n", col))
    print(table(hu@meta.data[[col]]))
  }
} else {
  cat("\n=== NO DEPOT COLUMNS FOUND ===\n")
}

cat("\nDone.\n")
