#!/usr/bin/env Rscript
# Final provenance/status capture for the pinned R environment: sessionInfo(),
# renv::status(), confirming the lockfile is in sync with the library. Run
# only after setup_pinned_r_env.R, fix_pass_2.R, the sccomp/scSHC retries,
# and the cmdstan install have all completed.
#
# Two bugs fixed here after an initial run produced a 4MB/140k-line
# FINAL_STATUS.md:
#  1. `remotes` was never actually installed INTO the isolated library --
#     every script that needed it found a copy via the broader .libPaths()
#     search before that search was restricted, so `requireNamespace()`
#     always returned TRUE and the explicit install was skipped. Force it in
#     here.
#  2. `print(renv::status(...))` under `sink()` printed the raw status object
#     (full DESCRIPTION-record lists for every package) instead of renv's
#     normal human-readable console summary. Fixed by capturing with
#     `capture.output()` and writing just the resulting text lines.
#  3. `renv::snapshot()` defaults to only recording packages renv's static
#     source-code scanner detects as "used" (a `library(pkg)` call somewhere
#     in the project). None of our install scripts have such calls -- they
#     install by string name through wrapper functions -- so the first
#     snapshot recorded almost nothing and every subsequent status() call
#     showed the whole library as "not installed". Fixed with
#     `type = "all"`, which snapshots everything actually present in the
#     library regardless of static usage detection.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
r_env_dir <- file.path(rerun_root, "r_env")
lib_dir <- file.path(r_env_dir, "renv_library")
.libPaths(lib_dir)
options(repos = c(CRAN = "https://cran.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE, lib.loc = lib_dir)) {
  install.packages("remotes", lib = lib_dir)
}
if (!requireNamespace("renv", quietly = TRUE, lib.loc = lib_dir)) {
  install.packages("renv", lib = lib_dir)
}
remotes_ok <- requireNamespace("remotes", quietly = TRUE, lib.loc = lib_dir)

out_file <- file.path(r_env_dir, "FINAL_STATUS.md")
sink(out_file)
cat("# Final R environment status\n\n")
cat("Captured:", format(Sys.time()), "\n\n")

cat("## R.version.string\n\n", R.version.string, "\n\n")

cat("## remotes installed into pinned library\n\n", remotes_ok, "\n\n")

cat("## sessionInfo()\n\n```\n")
print(sessionInfo())
cat("```\n\n")

cat("## renv::status() (post type=\"all\" snapshot)\n\n```\n")
snap_ok <- tryCatch({
  renv::snapshot(project = r_env_dir, library = lib_dir, type = "all", prompt = FALSE, force = TRUE)
  TRUE
}, error = function(e) { cat("renv::snapshot() FAILED:", conditionMessage(e), "\n"); FALSE })

status_text <- tryCatch(
  capture.output(renv::status(project = r_env_dir)),
  error = function(e) paste("renv::status() FAILED:", conditionMessage(e))
)
cat(paste(status_text, collapse = "\n"), "\n")
cat("```\n\n")

cat("## renv.lock present\n\n")
lockfile <- file.path(r_env_dir, "renv.lock")
cat(lockfile, ":", file.exists(lockfile), "\n\n")

sink()
cat("Wrote", out_file, "\n")
cat("snapshot ok:", snap_ok, "\n")
