#!/usr/bin/env Rscript
# Second fix pass for the pinned R environment install (see
# setup_pinned_r_env.R and its RERUN_LOG.md checkpoint). The first pass
# completed (exit 0:0) but several packages silently failed underneath a
# tryCatch that didn't catch remotes/BiocManager's warning-based failures.
# Root causes, diagnosed from the raw SLURM stdout log:
#
#   - monocle3: missing CRAN dependency Matrix.utils (archived from CRAN).
#   - sccomp: its own dependency `instantiate` was skipped as "already
#     installed" at a stale version (0.2.2.9001 < the >=0.2.3 sccomp needs)
#     because remotes::install_github() skips reinstall when the remote
#     SHA hasn't changed since a prior attempt; needs force = TRUE.
#   - scSHC -> scry -> beachmat: beachmat's compiled .so (found via the
#     admin-installed system Bioconductor library, which several packages
#     were silently satisfied from instead of lib_dir -- see below) requires
#     GLIBCXX_3.4.29, newer than what gcc/8.5.0 (the module-default compiler)
#     provides at runtme.
#   - celda -> scran -> bluster: bluster's C++ fails to compile against the
#     `assorthead` package's bundled tatami headers -- a class-template-
#     argument-deduction pattern gcc 8.5.0's C++17 support doesn't handle.
#   - SeuratWrappers -> Banksy -> aricode: aricode's C++ needs the
#     <execution> header (C++17 parallel STL), not available with gcc 8.5.0's
#     libstdc++.
#   - clusterProfiler, celda, org.Mm.eg.db, Biobase were silently SKIPPED by
#     BiocManager::install() ("package(s) not installed when version(s) same
#     as or greater than current") because a satisfactory version already
#     existed somewhere else in .libPaths() (the admin-installed system Bioc
#     library and/or Xinru's personal ~/R library) -- meaning the "isolated"
#     library was actually leaking on packages that already existed
#     elsewhere. Not true isolation.
#
# Fix: (1) restrict .libPaths() to ONLY lib_dir (base/recommended packages
# still resolve via .Library, R's home library, which is separate from
# .libPaths()) so nothing is silently satisfied by the system or personal
# library; (2) compile with gcc/12.2.0 (module-available, newer libstdc++
# covers the GLIBCXX/CTAD/<execution> gaps) instead of the module-default
# gcc/8.5.0; (3) force = TRUE / update = TRUE throughout so nothing is
# skipped as "already installed" against a library this script isn't
# actually managing.
#
# pak/curl/devtools/usethis are explicitly NOT fixed here: curl's C source
# needs libcurl symbols (CURL_HTTP_VERSION_3, CURLINFO_EFFECTIVE_METHOD) not
# present in this node's system libcurl-devel headers, a separate system
# library gap unrelated to the compiler. devtools is not in the
# investigator's required 12-package list and is not actually used
# functionally anywhere in scripts/ (only named in a help-text string in
# 04_eos_subcluster.R:370), so this gap is documented as non-blocking rather
# than chased further.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
r_env_dir <- file.path(rerun_root, "r_env")
lib_dir <- file.path(r_env_dir, "renv_library")

# Restrict to ONLY the pinned library -- do not let installs be silently
# satisfied by the admin system Bioc library or Xinru's personal ~/R library.
.libPaths(lib_dir)
Sys.setenv(R_LIBS_USER = lib_dir)
options(repos = c(CRAN = "https://cran.r-project.org"))

log_file <- file.path(r_env_dir, "install_log_pass2.txt")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

status <- list()
safe_install <- function(name, expr) {
  log_msg("Installing:", name)
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    log_msg("FAILED (error):", name, "--", conditionMessage(e))
    FALSE
  })
  status[[name]] <<- ok
  ok
}
# Verify against the library itself, not the tryCatch return value alone --
# remotes/BiocManager report some failures as warnings, not errors.
verify_installed <- function(pkg, expected_version = NULL) {
  v <- tryCatch(as.character(packageVersion(pkg, lib.loc = lib_dir)), error = function(e) NA_character_)
  ok <- !is.na(v) && (is.null(expected_version) || v == expected_version)
  log_msg("VERIFY", pkg, ":", if (is.na(v)) "NOT INSTALLED" else v,
          if (!is.null(expected_version)) paste0(" (expected ", expected_version, ")") else "",
          "--", if (ok) "OK" else "FAIL")
  ok
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib_dir)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib_dir)

# ---- monocle3: install missing Matrix.utils from CRAN archive -------------
safe_install("Matrix.utils (CRAN archive, monocle3 dep)",
             remotes::install_version("Matrix.utils", version = "0.9.8", lib = lib_dir, upgrade = "never", force = TRUE))
safe_install("monocle3@1543aac6 retry", remotes::install_github(
  "cole-trapnell-lab/monocle3", ref = "1543aac6ee45cc99354bc51fbbbcfda9fe83fe61",
  lib = lib_dir, upgrade = "never", force = TRUE))
verify_installed("monocle3", "1.2.9")

# ---- sccomp: force instantiate >=0.2.3, then retry sccomp -----------------
safe_install("instantiate (force, sccomp dep)",
             remotes::install_github("stemangiola/instantiate", lib = lib_dir, upgrade = "never", force = TRUE))
safe_install("sccomp@ec575598 retry", remotes::install_github(
  "stemangiola/sccomp", ref = "ec575598df7f633623bf577e1d8d7fc0adcdb34d",
  lib = lib_dir, upgrade = "never", force = TRUE))
verify_installed("sccomp", "1.7.16")

# ---- scSHC -> scry -> beachmat: force fresh compile with gcc/12.2.0 -------
# (CC/CXX are set via Makevars written by the submitting shell script, not
# here, since R's compiler selection is read at package-build subprocess
# launch time from ~/.R/Makevars.)
safe_install("beachmat (force fresh compile)",
             BiocManager::install("beachmat", lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
safe_install("scry (force fresh compile)",
             BiocManager::install("scry", lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
safe_install("scSHC@2671f6ee retry", remotes::install_github(
  "igrabski/sc-SHC", ref = "2671f6ee5a55572ddd7216de822044f9497d6161",
  lib = lib_dir, upgrade = "never", force = TRUE))
verify_installed("scSHC")

# ---- celda -> scran -> bluster: force fresh compile with gcc/12.2.0 -------
safe_install("bluster (force fresh compile)",
             BiocManager::install("bluster", lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
safe_install("scran (force fresh compile)",
             BiocManager::install("scran", lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
safe_install("celda retry", BiocManager::install("celda", lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
verify_installed("celda")

# ---- clusterProfiler, org.Mm.eg.db, Biobase: force fresh install into the
# now-isolated lib_dir (previously silently satisfied elsewhere) -----------
for (p in c("clusterProfiler", "org.Mm.eg.db", "Biobase")) {
  safe_install(paste0(p, " (force, isolated)"),
               BiocManager::install(p, lib = lib_dir, update = FALSE, ask = FALSE, force = TRUE))
  verify_installed(p)
}

# ---- SeuratWrappers -> Banksy -> aricode: force fresh compile with gcc/12 --
safe_install("aricode (force fresh compile)",
             remotes::install_version("aricode", version = "1.1.0", lib = lib_dir, upgrade = "never", force = TRUE))
safe_install("leidenAlg (force)",
             remotes::install_version("leidenAlg", lib = lib_dir, upgrade = "never", force = TRUE))
safe_install("Banksy (force)", remotes::install_github(
  "prabhakarlab/Banksy", lib = lib_dir, upgrade = "never", force = TRUE))
safe_install("SeuratWrappers retry", remotes::install_github(
  "satijalab/seurat-wrappers", lib = lib_dir, upgrade = "never", force = TRUE))
verify_installed("SeuratWrappers")

# ---- Final verification sweep for everything this pass touched ------------
cat("\n=== Pass 2 verification summary ===\n")
final_check <- c("Matrix.utils", "monocle3", "instantiate", "sccomp", "beachmat", "scry",
                  "scSHC", "bluster", "scran", "celda", "clusterProfiler", "org.Mm.eg.db",
                  "Biobase", "aricode", "leidenAlg", "Banksy", "SeuratWrappers")
all_ok <- TRUE
for (p in final_check) {
  ok <- verify_installed(p)
  all_ok <- all_ok && ok
}
cat(if (all_ok) "ALL PASS-2 TARGETS INSTALLED\n" else "SOME PASS-2 TARGETS STILL MISSING -- see log\n")

log_msg("Fix pass 2 complete.")
if (!all_ok) quit(status = 1) else quit(status = 0)
