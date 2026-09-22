#!/usr/bin/env Rscript
# scSHC (igrabski/sc-SHC@2671f6ee) failed its first install attempt because
# its dependency `scry` (Bioconductor) was not available. Install scry first,
# then retry scSHC into the same pinned library.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
lib_dir <- file.path(rerun_root, "r_env", "renv_library")
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cran.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib_dir)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib_dir)

cat("Installing scry (Bioconductor, scSHC dependency)...\n")
ok_scry <- tryCatch({
  BiocManager::install("scry", lib = lib_dir, update = FALSE, ask = FALSE)
  TRUE
}, error = function(e) { message("scry FAILED: ", conditionMessage(e)); FALSE })
cat("scry install:", if (ok_scry) "OK" else "FAILED", "\n")

cat("Retrying scSHC@2671f6ee5a55572ddd7216de822044f9497d6161...\n")
ok_scshc <- tryCatch({
  remotes::install_github("igrabski/sc-SHC", ref = "2671f6ee5a55572ddd7216de822044f9497d6161", lib = lib_dir, upgrade = "never")
  TRUE
}, error = function(e) { message("scSHC FAILED: ", conditionMessage(e)); FALSE })
cat("scSHC install:", if (ok_scshc) "OK" else "FAILED", "\n")

v <- tryCatch(as.character(packageVersion("scSHC", lib.loc = lib_dir)), error = function(e) "NOT INSTALLED")
cat("scSHC installed version:", v, "\n")

if (!ok_scry || !ok_scshc) quit(status = 1) else quit(status = 0)
