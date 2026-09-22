#!/usr/bin/env Rscript
# Pinned, isolated R environment for the audited HPCC rerun.
#
# Built per investigator directive (2026-09-22): this must be a genuine
# reproduction, not a version-drift sensitivity analysis. sccomp 2.0.0 vs the
# documented 1.7, and a completely missing scSHC, could change scientific
# conclusions (compositional testing, clustering-resolution justification).
#
# Runs under R 4.5.0 (module load R/4.5.0 after `module purge` -- a bare
# `module load R` on this cluster silently resolves to the newer default
# R/4.5.2 instead). Entirely independent from the shared library at
# /bigdata/godziklab/shared/Xinru/R, which this script never writes to.
# Library lives under $RERUN_ROOT/r_env/renv_library.
#
# Version pins were resolved by walking each package's GitHub commit history
# for the commit where DESCRIPTION's Version field matches the README's
# documented published-run pin, since none of CellChat/monocle3/sccomp have a
# git tag at the exact documented version:
#   - CellChat:  sqjin/CellChat@e4f68625 -- DESCRIPTION Version: 1.6.1,
#                unchanged since 2022-11-30 (repo's final commit, 2023-11-06)
#   - monocle3:  cole-trapnell-lab/monocle3@1543aac6 -- DESCRIPTION
#                Version: 1.2.9, the commit immediately before the bump to
#                1.2.10 three minutes later
#   - sccomp:    stemangiola/sccomp@ec575598 -- DESCRIPTION Version: 1.7.16.
#                README only documents "sccomp 1.7" with no patch number;
#                this is the last 1.7.x commit before the repo moved to
#                1.9.x -- a BEST-EFFORT pin, not a verified exact match.
#   - scSHC:     igrabski/sc-SHC@2671f6ee -- this repo's only-ever commit
#                (unchanged since 2024-03-01), so there is no ambiguity.
#   - Seurat 5.1.0, harmony 1.2.1: exact versions available on the CRAN
#     archive via remotes::install_version().
# All other packages (DESeq2, clusterProfiler, celda, clustree, etc.) have no
# exact version documented anywhere in the repo; installed unpinned against
# whatever CRAN/Bioconductor serves for R 4.5.0 at build time.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
r_env_dir <- file.path(rerun_root, "r_env")
lib_dir <- file.path(r_env_dir, "renv_library")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(lib_dir, .libPaths()))
Sys.setenv(R_LIBS_USER = lib_dir)
options(repos = c(CRAN = "https://cran.r-project.org"))

log_file <- file.path(r_env_dir, "install_log.txt")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

install_status <- list()
safe_install <- function(name, expr) {
  log_msg("Installing:", name)
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    log_msg("FAILED:", name, "--", conditionMessage(e))
    FALSE
  })
  install_status[[name]] <<- ok
  ok
}

# ---- Bootstrap: remotes + BiocManager + renv (for lockfile capture) --------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib_dir)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib_dir)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", lib = lib_dir)

library(remotes)
log_msg("Bioconductor version for this R:", as.character(BiocManager::version()))

# ---- Exact-pinned packages --------------------------------------------------
safe_install("Seurat 5.1.0", remotes::install_version("Seurat", version = "5.1.0", lib = lib_dir, upgrade = "never"))
safe_install("harmony 1.2.1", remotes::install_version("harmony", version = "1.2.1", lib = lib_dir, upgrade = "never"))
safe_install("CellChat@e4f68625 (1.6.1)", remotes::install_github("sqjin/CellChat", ref = "e4f68625b074247d619c2e488d33970cc531e17c", lib = lib_dir, upgrade = "never"))
safe_install("monocle3@1543aac6 (1.2.9)", remotes::install_github("cole-trapnell-lab/monocle3", ref = "1543aac6ee45cc99354bc51fbbbcfda9fe83fe61", lib = lib_dir, upgrade = "never"))
safe_install("sccomp@ec575598 (1.7.16, best-effort for documented '1.7')", remotes::install_github("stemangiola/sccomp", ref = "ec575598df7f633623bf577e1d8d7fc0adcdb34d", lib = lib_dir, upgrade = "never"))
safe_install("scSHC@2671f6ee (HEAD, only-ever version)", remotes::install_github("igrabski/sc-SHC", ref = "2671f6ee5a55572ddd7216de822044f9497d6161", lib = lib_dir, upgrade = "never"))

# ---- Bioconductor packages (no exact pin documented; use R 4.5.0's default Bioc release) ----
bioc_pkgs <- c("DESeq2", "clusterProfiler", "celda", "ComplexHeatmap", "edgeR",
               "enrichplot", "org.Mm.eg.db", "Biobase", "SingleCellExperiment")
for (p in bioc_pkgs) {
  safe_install(paste0(p, " (Bioconductor, unpinned)"), BiocManager::install(p, lib = lib_dir, update = FALSE, ask = FALSE))
}

# ---- Remaining CRAN packages (no exact pin documented) ---------------------
cran_pkgs <- c("SeuratWrappers", "clustree", "circlize", "clue", "cluster",
               "cowplot", "dendextend", "devtools", "dplyr", "ggplot2", "ggrepel",
               "gridExtra", "igraph", "mclust", "openxlsx", "patchwork", "pheatmap",
               "reticulate", "scales", "stringr", "tibble", "tidyr", "RColorBrewer",
               "RANN", "Matrix")
for (p in cran_pkgs) {
  if (p == "SeuratWrappers") {
    safe_install(p, remotes::install_github("satijalab/seurat-wrappers", lib = lib_dir, upgrade = "never"))
  } else {
    safe_install(paste0(p, " (CRAN, unpinned)"), install.packages(p, lib = lib_dir))
  }
}

# ---- Provenance capture -----------------------------------------------------
prov_file <- file.path(r_env_dir, "PROVENANCE.md")
sink(prov_file)
cat("# R environment provenance\n\n")
cat("Built:", format(Sys.time()), "\n\n")
cat("## sessionInfo()\n\n```\n")
print(sessionInfo())
cat("```\n\n")
cat("## Bioconductor version\n\n", as.character(BiocManager::version()), "\n\n")
cat("## Pinned package install sources\n\n")
cat("- Seurat 5.1.0: CRAN archive (install_version)\n")
cat("- harmony 1.2.1: CRAN archive (install_version)\n")
cat("- CellChat: GitHub sqjin/CellChat @ e4f68625b074247d619c2e488d33970cc531e17c (DESCRIPTION Version: 1.6.1, unchanged since 2022-11-30)\n")
cat("- monocle3: GitHub cole-trapnell-lab/monocle3 @ 1543aac6ee45cc99354bc51fbbbcfda9fe83fe61 (DESCRIPTION Version: 1.2.9)\n")
cat("- sccomp: GitHub stemangiola/sccomp @ ec575598df7f633623bf577e1d8d7fc0adcdb34d (DESCRIPTION Version: 1.7.16 -- README only documents 'sccomp 1.7' with no patch number; this is the last 1.7.x commit before the repo moved to 1.9.x, a best-effort pin, NOT a verified exact match)\n")
cat("- scSHC: GitHub igrabski/sc-SHC @ 2671f6ee5a55572ddd7216de822044f9497d6161 (repo's only-ever commit as of query time; no version ambiguity)\n")
cat("- All other packages: unpinned, installed at whatever CRAN/Bioconductor served for R 4.5.0 / this Bioc release at build time (no exact version documented anywhere in the repo to reproduce against)\n\n")
cat("## Installed package versions (actual)\n\n")
all_pkgs <- c("Seurat", "SeuratObject", "SeuratWrappers", "CellChat", "monocle3", "sccomp",
              "harmony", "scSHC", "DESeq2", "clusterProfiler", "celda", "ComplexHeatmap",
              "edgeR", "enrichplot", "org.Mm.eg.db", "Biobase", "SingleCellExperiment",
              "clustree", "circlize", "clue", "cluster", "cowplot", "dendextend", "devtools",
              "dplyr", "ggplot2", "ggrepel", "gridExtra", "igraph", "mclust", "openxlsx",
              "patchwork", "pheatmap", "reticulate", "scales", "stringr", "tibble", "tidyr",
              "RColorBrewer", "RANN", "Matrix")
for (p in all_pkgs) {
  v <- tryCatch(as.character(packageVersion(p, lib.loc = lib_dir)), error = function(e) "NOT INSTALLED")
  cat("-", p, ":", v, "\n")
}
cat("\n## Install success/failure per safe_install() call\n\n")
for (nm in names(install_status)) {
  cat("-", nm, ":", if (install_status[[nm]]) "OK" else "FAILED", "\n")
}
sink()

# ---- renv lockfile (best-effort; renv wasn't used to drive install, so this
# snapshots the resulting library rather than governing it) -----------------
tryCatch({
  renv::init(project = r_env_dir, bare = TRUE, restart = FALSE)
  .libPaths(c(lib_dir, .libPaths()))
  renv::snapshot(project = r_env_dir, library = lib_dir, prompt = FALSE, force = TRUE)
}, error = function(e) log_msg("renv snapshot FAILED:", conditionMessage(e)))

log_msg("Setup script complete. See", prov_file)
