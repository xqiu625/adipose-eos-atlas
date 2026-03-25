#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## QC filtering + doublet removal on DecontX-cleaned cells
## SLURM: --mem=128G --cpus-per-task=8 --time=4:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(reticulate)
})

# ── Shared theme ─────────────────────────────────────────────────────────────
SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
CONTAM_THRESHOLD <- 0.3

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

# ── Load Scrublet via reticulate ───────────────────────────────────────────────
cat("Importing Scrublet via reticulate...\n")
use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env", required = TRUE)
scrublet <- tryCatch({
  import("scrublet")
}, error = function(e) {
  stop("Failed to import scrublet: ", e$message,
       "\nMake sure scrublet is installed in the active Python environment.")
})
scipy_sparse <- import("scipy.sparse")
np <- import("numpy")

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
    counts <- decon$counts
    contam <- decon$contamination

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
      condition = cond,
      sample_id = sid
    )

    total_before <- total_before + n_before
    total_after  <- total_after + n_after

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

cat(sprintf("\nTotal cells before contam filter: %d\n", total_before))
cat(sprintf("Total cells after contam filter:  %d\n\n", total_after))

# ── Step 2: Merge into one Seurat object ──────────────────────────────────────
cat("========================================\n")
cat("Step 2: Merge into Seurat object\n")
cat("========================================\n\n")

# Combine counts from all samples
all_counts <- do.call(cbind, lapply(sample_objects, function(x) x$counts))
cat(sprintf("Combined matrix: %d genes x %d cells\n", nrow(all_counts), ncol(all_counts)))

# Create Seurat object
sobj <- CreateSeuratObject(counts = all_counts, project = "adipose_atlas")

# Assign metadata directly via @meta.data to avoid Seurat v5 name-matching issues
# Barcodes are: condition_originalbarcode (e.g. WT_Male_HFD_1_ACGTACGT-1)
# Extract condition by removing barcode suffix after the last underscore
ident_vec <- sub("_[^_]+$", "", colnames(sobj))
sobj@meta.data$orig.ident <- ident_vec

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

# ── Step 5: Doublet detection with Scrublet ───────────────────────────────────
cat("========================================\n")
cat("Step 4: Doublet detection (Scrublet)\n")
cat("========================================\n\n")

samples <- unique(sobj$orig.ident)
doublet_scores <- rep(NA_real_, ncol(sobj))
doublet_calls  <- rep(NA, ncol(sobj))
names(doublet_scores) <- colnames(sobj)
names(doublet_calls)  <- colnames(sobj)

for (s in samples) {
  cat(sprintf("[%s] ", s))

  cells_idx <- which(sobj$orig.ident == s)
  n_cells <- length(cells_idx)
  cat(sprintf("%d cells -> ", n_cells))

  tryCatch({
    # Extract counts for this sample
    counts_sub <- GetAssayData(sobj, layer = "counts")[, cells_idx]

    # Convert to scipy sparse matrix (CSC format expected by Scrublet)
    counts_t <- t(counts_sub)  # Scrublet expects cells x genes
    counts_scipy <- scipy_sparse$csc_matrix(
      tuple(
        np$array(counts_t@x),
        np$array(as.integer(counts_t@i)),
        np$array(as.integer(counts_t@p))
      ),
      shape = tuple(as.integer(nrow(counts_t)), as.integer(ncol(counts_t)))
    )

    # Run Scrublet
    scrub <- scrublet$Scrublet(counts_scipy, expected_doublet_rate = 0.06)
    result <- scrub$scrub_doublets()

    scores <- result[[1]]
    predicted <- result[[2]]

    doublet_scores[cells_idx] <- as.numeric(scores)
    doublet_calls[cells_idx]  <- as.logical(predicted)

    n_doublets <- sum(as.logical(predicted))
    cat(sprintf("%d doublets (%.1f%%)\n", n_doublets,
                100 * n_doublets / n_cells))

  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
  })
}

# Add to metadata
sobj$scrublet_score   <- doublet_scores
sobj$scrublet_doublet <- doublet_calls

# Remove predicted doublets
n_before_doublet <- ncol(sobj)
sobj <- subset(sobj, cells = colnames(sobj)[!sobj$scrublet_doublet])
n_after_doublet <- ncol(sobj)

cat(sprintf("\nCells before doublet removal: %d\n", n_before_doublet))
cat(sprintf("Cells after doublet removal:  %d (removed %d, %.1f%%)\n",
            n_after_doublet, n_before_doublet - n_after_doublet,
            100 * (n_before_doublet - n_after_doublet) / n_before_doublet))

# ── Step 6: Final summary and save ────────────────────────────────────────────
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
