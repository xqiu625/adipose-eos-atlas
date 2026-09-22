#!/usr/bin/env Rscript
# Final provenance/status capture for the pinned R environment: sessionInfo(),
# renv::status(), confirming the lockfile is in sync with the library. Run
# only after setup_pinned_r_env.R (and fix_scshc_install.R) have completed.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
r_env_dir <- file.path(rerun_root, "r_env")
lib_dir <- file.path(r_env_dir, "renv_library")
.libPaths(c(lib_dir, .libPaths()))

out_file <- file.path(r_env_dir, "FINAL_STATUS.md")
sink(out_file)
cat("# Final R environment status\n\n")
cat("Captured:", format(Sys.time()), "\n\n")

cat("## R.version.string\n\n", R.version.string, "\n\n")

cat("## sessionInfo()\n\n```\n")
print(sessionInfo())
cat("```\n\n")

cat("## renv::status()\n\n```\n")
tryCatch({
  renv::snapshot(project = r_env_dir, library = lib_dir, prompt = FALSE, force = TRUE)
  print(renv::status(project = r_env_dir))
}, error = function(e) cat("renv::status() FAILED:", conditionMessage(e), "\n"))
cat("```\n\n")

cat("## renv.lock present\n\n")
lockfile <- file.path(r_env_dir, "renv.lock")
cat(lockfile, ":", file.exists(lockfile), "\n\n")

sink()
cat("Wrote", out_file, "\n")
