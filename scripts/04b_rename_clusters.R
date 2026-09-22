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

## 04b_rename_clusters.R -- RETIRED, do not run.
##
## This script used to apply a second, independently hardcoded Seurat-label ->
## manuscript-label permutation to eos_subclustered.rds and overwrite it.
## 04_eos_subcluster.R's "Step 5c: Rename clusters to match manuscript
## numbering" now performs this renaming itself, derived from this run's own
## barcode-matched correspondence to the old manuscript object (not a fixed
## historical map), and already re-runs FindAllMarkers, regenerates the
## UMAP/composition/marker figures, and saves eos_subclustered.rds and
## eos_cluster_markers.csv with the manuscript labels applied.
##
## Running this script after 04_eos_subcluster.R would re-load the
## already-correctly-labeled object and apply 04b's OWN different fixed
## rename_map on top of it -- silently permuting an already-renamed C2-C5 a
## second time. There is no valid use of this script anymore; it stops
## immediately rather than doing that.

stop(
  "04b_rename_clusters.R is retired: 04_eos_subcluster.R now performs ",
  "manuscript-label renaming itself (Step 5c), derived from this run's own ",
  "cluster correspondence. Do not run this script -- it would apply a second, ",
  "different, hardcoded permutation on top of already-renamed clusters and ",
  "corrupt eos_subclustered.rds. See scripts/04_eos_subcluster.R."
)
