#!/usr/bin/env Rscript
# End-to-end test that a real pipeline script's exact lib_path/.libPaths()
# pattern, run under RERUN_R_LIB, resolves packages from the pinned renv
# library and NOT from the shared library at
# /bigdata/godziklab/shared/Xinru/R -- the real failure mode this guards
# against is the shared library silently winning because a script re-prepends
# it after some other environment-level fix already set .libPaths().

# This mirrors scripts/04_eos_subcluster.R lines 1-3 exactly.
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

cat("RERUN_R_LIB env var:", Sys.getenv("RERUN_R_LIB"), "\n")
cat("Resolved lib_path:", lib_path, "\n")
cat(".libPaths():\n")
print(.libPaths())

shared_lib <- "/bigdata/godziklab/shared/Xinru/R"
renv_lib <- Sys.getenv("RERUN_R_LIB")

stopifnot("lib_path must equal RERUN_R_LIB, not fall back to shared lib" = lib_path == renv_lib)
stopifnot("renv lib must be first in .libPaths()" = .libPaths()[1] == renv_lib)

check_pkg_source <- function(pkg) {
  path <- tryCatch(find.package(pkg), error = function(e) NA_character_)
  from_shared <- !is.na(path) && startsWith(path, shared_lib)
  from_renv <- !is.na(path) && startsWith(path, renv_lib)
  cat(sprintf("  %-20s path=%s  from_renv=%s  from_shared=%s\n",
              pkg, path, from_renv, from_shared))
  list(pkg = pkg, path = path, from_renv = from_renv, from_shared = from_shared)
}

cat("\n=== Package source resolution ===\n")
pkgs <- c("Seurat", "harmony", "CellChat", "monocle3", "sccomp", "scSHC", "DESeq2")
results <- lapply(pkgs, check_pkg_source)

any_from_shared <- any(vapply(results, function(r) isTRUE(r$from_shared), logical(1)))
all_from_renv <- all(vapply(results, function(r) isTRUE(r$from_renv), logical(1)))

cat("\n=== Verdict ===\n")
if (any_from_shared) {
  cat("FAIL: at least one package silently resolved from the shared library.\n")
  quit(status = 1)
} else if (!all_from_renv) {
  cat("FAIL: at least one package did not resolve from the renv library (missing?).\n")
  quit(status = 1)
} else {
  cat("PASS: all checked packages resolve from the pinned renv library; shared library is not shadowing anything.\n")
  quit(status = 0)
}
