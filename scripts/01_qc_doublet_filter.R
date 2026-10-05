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

## Cell QC and merge after pre-DecontX doublet removal
## SLURM: --mem=128G --cpus-per-task=8 --time=4:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ── Shared theme ─────────────────────────────────────────────────────────────
SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- Sys.getenv(
  "EOS_ATLAS_BASE",
  unset = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
)
CONTAM_THRESHOLD <- 0.2

# Sample-to-condition mapping
sample_info <- data.frame(
  sample_id = 1:16,
  condition = c("WT_Male_HFD_1",       "WT_Female_HFD_1",
                "KO_Male_HFD_1",       "KO_Female_HFD_1",
                "WT_Female_HFDInfected_1", "WT_Female_HFDInfected_2",
                "KO_Female_HFDInfected_1", "KO_Female_HFDInfected_2",
                "WT_Male_HFD_2",       "WT_Male_HFD_3",
                "WT_Female_HFD_2",     "WT_Female_HFD_3",
                "KO_Male_HFD_2",       "KO_Male_HFD_3",
                "KO_Female_HFD_2",     "KO_Female_HFD_3"),
  stringsAsFactors = FALSE
)

# GEM reaction mapping
gem_lookup <- c(rep("GEM1", 4), rep("GEM2", 4), rep("GEM3", 4), rep("GEM4", 4))

# ── Helper: parse metadata from condition string ──────────────────────────────
parse_metadata <- function(cond) {
  # e.g. "WT_Male_HFD_1" -> genotype=WT, sex=Male, infection=HFD, replicate=1
  # e.g. "KO_Female_HFDInfected_2" -> genotype=KO, sex=Female, infection=HFDInfected
  parts <- strsplit(cond, "_")[[1]]
  data.frame(
    genotype  = parts[1],
    sex       = parts[2],
    infection = parts[3],
    stringsAsFactors = FALSE
  )
}

# ── Step 1: Load DecontX results and filter by contamination ──────────────────
cat("========================================\n")
cat("Step 1: Load DecontX results & filter contamination\n")
cat("========================================\n\n")
cat("Analysis base:", BASE, "\n\n")

decontx_dir <- file.path(BASE, "decontx")
sample_objects <- list()

total_before <- 0
total_after  <- 0

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  cond <- sample_info$condition[i]

  rds_file <- file.path(decontx_dir, sprintf("decontx_sample_%d.rds", sid))
  cat(sprintf("[Sample %d] %s\n", sid, cond))

  tryCatch({
    decon <- readRDS(rds_file)
    if (!isTRUE(decon$doublet_removed_before_decontx)) {
      stop("Input was not generated after pre-DecontX doublet removal: ", rds_file)
    }
    counts <- decon$counts
    contam <- decon$contamination
    scrublet_score <- decon$scrublet_score
    if (is.null(scrublet_score) || length(scrublet_score) != ncol(counts)) {
      stop("Missing or misaligned Scrublet scores in ", rds_file)
    }
    if (is.null(names(scrublet_score)) ||
        !identical(names(scrublet_score), colnames(counts))) {
      stop("Scrublet score names do not match count-matrix columns in ", rds_file)
    }

    n_before <- ncol(counts)
    keep <- contam <= CONTAM_THRESHOLD
    n_after <- sum(keep)

    cat(sprintf("  Cells before contam filter: %d\n", n_before))
    cat(sprintf("  Cells after contam filter:  %d (removed %d, %.1f%%)\n",
                n_after, n_before - n_after,
                100 * (n_before - n_after) / n_before))

    counts_clean <- counts[, keep, drop = FALSE]

    # Rename barcodes to avoid collisions across samples
    colnames(counts_clean) <- paste0(cond, "_", colnames(counts_clean))

    sample_objects[[as.character(sid)]] <- list(
      counts    = counts_clean,
      scrublet_score = scrublet_score[keep],
      condition = cond,
      sample_id = sid
    )

    total_before <- total_before + n_before
    total_after  <- total_after + n_after

  }, error = function(e) {
    stop(sprintf("ERROR processing sample %d (%s): %s", sid, cond, e$message),
         call. = FALSE)
  })
}

if (length(sample_objects) != nrow(sample_info)) {
  stop("QC input did not complete for all 16 samples; refusing to continue.")
}

cat(sprintf("\nTotal cells before contam filter: %d\n", total_before))
cat(sprintf("Total cells after contam filter:  %d\n\n", total_after))

# ── Step 2: Merge into one Seurat object ──────────────────────────────────────
cat("========================================\n")
cat("Step 2: Merge into Seurat object\n")
cat("========================================\n\n")

# Combine counts from all samples
# cbind() does not align rows by gene name. Refuse silent gene permutation.
reference_genes <- rownames(sample_objects[[1]]$counts)
if (is.null(reference_genes) || anyDuplicated(reference_genes)) {
  stop("Missing or duplicated reference gene identifiers before sample merge")
}
for (sid in names(sample_objects)) {
  if (!identical(rownames(sample_objects[[sid]]$counts), reference_genes)) {
    stop("Gene rows differ in identity/order for sample ", sid,
         "; explicitly align features upstream before merging")
  }
}
all_counts <- do.call(cbind, lapply(sample_objects, function(x) x$counts))
cat(sprintf("Combined matrix: %d genes x %d cells\n", nrow(all_counts), ncol(all_counts)))

# Create Seurat object
sobj <- CreateSeuratObject(counts = all_counts, project = "adipose_atlas")

# Assign metadata directly via @meta.data to avoid Seurat v5 name-matching issues
# Barcodes are: condition_originalbarcode (e.g. WT_Male_HFD_1_ACGTACGT-1)
# Extract condition by removing barcode suffix after the last underscore
ident_vec <- sub("_[^_]+$", "", colnames(sobj))
sobj@meta.data$orig.ident <- ident_vec
sobj$scrublet_score <- unlist(lapply(sample_objects, function(x) x$scrublet_score),
                                use.names = FALSE)
sobj$scrublet_doublet <- FALSE

cat(sprintf("Seurat object: %d genes x %d cells\n", nrow(sobj), ncol(sobj)))

# ── Step 3: Add metadata ──────────────────────────────────────────────────────
cat("\nAdding metadata columns...\n")

# Condition = orig.ident without trailing replicate number
sobj@meta.data$condition <- sub("_\\d+$", "", ident_vec)

# Parse genotype, sex, infection from condition string
meta_parsed <- do.call(rbind, lapply(ident_vec, parse_metadata))
sobj@meta.data$genotype  <- meta_parsed$genotype
sobj@meta.data$sex       <- meta_parsed$sex
sobj@meta.data$infection <- meta_parsed$infection

# GEM reaction
sample_to_gem <- setNames(gem_lookup, sample_info$condition)
sobj@meta.data$gem_reaction <- unname(sample_to_gem[ident_vec])

cat("Metadata columns added: condition, genotype, sex, infection, gem_reaction\n")
cat(sprintf("Samples: %s\n", paste(unique(sobj$orig.ident), collapse = ", ")))

# ── Step 4: QC filtering ──────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Step 3: QC filtering\n")
cat("========================================\n\n")

n_before_qc <- ncol(sobj)
cat(sprintf("Cells before QC: %d\n", n_before_qc))

# Compute percent mitochondrial
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")

# Apply QC filters
keep_qc <- sobj$nFeature_RNA >= 200 &
           sobj$nFeature_RNA <= 5000 &
           sobj$percent.mt <= 10

cat(sprintf("  nFeature_RNA >= 200:  %d pass\n", sum(sobj$nFeature_RNA >= 200)))
cat(sprintf("  nFeature_RNA <= 5000: %d pass\n", sum(sobj$nFeature_RNA <= 5000)))
cat(sprintf("  percent.mt <= 10:     %d pass\n", sum(sobj$percent.mt <= 10)))
cat(sprintf("  All QC filters:       %d pass\n", sum(keep_qc)))

sobj <- subset(sobj, cells = colnames(sobj)[keep_qc])
n_after_qc <- ncol(sobj)
cat(sprintf("Cells after QC: %d (removed %d, %.1f%%)\n\n",
            n_after_qc, n_before_qc - n_after_qc,
            100 * (n_before_qc - n_after_qc) / n_before_qc))

# ── Step 5: Final summary and save ────────────────────────────────────────────
cat("\n========================================\n")
cat("Final Summary\n")
cat("========================================\n")
cat(sprintf("Final Seurat object: %d genes x %d cells\n", nrow(sobj), ncol(sobj)))
cat("\nCells per sample:\n")
print(table(sobj$orig.ident))
cat("\nCells per condition:\n")
print(table(sobj$condition))

out_file <- file.path(BASE, "seurat_qc_filtered.rds")
cat(sprintf("\nSaving: %s\n", out_file))
saveRDS(sobj, out_file)

cat("Done.\n")
