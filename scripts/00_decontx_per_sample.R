#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## DecontX ambient RNA removal per sample, after pre-DecontX doublet removal
## SLURM: --mem=128G --cpus-per-task=8 --time=4:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(celda)
  library(Matrix)
})

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- Sys.getenv(
  "EOS_ATLAS_BASE",
  unset = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
)
CONTAM_THRESHOLD <- 0.2
DOUBLET_DIR <- file.path(BASE, "doublets")

# Sample-to-condition mapping (16 samples)
sample_info <- data.frame(
  sample_id   = 1:16,
  condition    = c("WT_Male_HFD_1",       "WT_Female_HFD_1",
                   "KO_Male_HFD_1",       "KO_Female_HFD_1",
                   "WT_Female_HFDInfected_1", "WT_Female_HFDInfected_2",
                   "KO_Female_HFDInfected_1", "KO_Female_HFDInfected_2",
                   "WT_Male_HFD_2",       "WT_Male_HFD_3",
                   "WT_Female_HFD_2",     "WT_Female_HFD_3",
                   "KO_Male_HFD_2",       "KO_Male_HFD_3",
                   "KO_Female_HFD_2",     "KO_Female_HFD_3"),
  stringsAsFactors = FALSE
)

# CellRanger data paths by GEM reaction
data_root <- Sys.getenv(
  "EOS_DATA_ROOT",
  unset = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data"
)
gem_map <- list(
  "RR0104" = 1:4,
  "RR0508" = 5:8,
  "RR0912" = 9:12,
  "RR1316" = 13:16
)

# Build per-sample paths
get_sample_paths <- function(sample_id) {
  gem <- NULL
  for (g in names(gem_map)) {
    if (sample_id %in% gem_map[[g]]) {
      gem <- g
      break
    }
  }
  base_path <- file.path(data_root, gem, "outs", "per_sample_outs",
                         as.character(sample_id), "count")
  list(
    filtered = file.path(base_path, "sample_filtered_feature_bc_matrix"),
    raw      = file.path(base_path, "sample_raw_feature_bc_matrix")
  )
}

# ── Output directory ───────────────────────────────────────────────────────────
out_dir <- file.path(BASE, "decontx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("DecontX ambient RNA removal per sample\n")
cat("========================================\n")
cat("Output directory:", out_dir, "\n")
cat("Cell Ranger data root:", data_root, "\n")
cat("Contamination threshold:", CONTAM_THRESHOLD, "\n\n")

# ── Run DecontX per sample ─────────────────────────────────────────────────────
all_scores <- list()

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  cond <- sample_info$condition[i]

  cat(sprintf("[Sample %d/%d] %s (ID: %d)\n", i, nrow(sample_info), cond, sid))

  paths <- get_sample_paths(sid)

  tryCatch({
    # Load filtered counts
    cat("  Loading filtered counts:", paths$filtered, "\n")
    filtered_counts <- Read10X(data.dir = paths$filtered)

    # Load raw counts (background)
    cat("  Loading raw counts:", paths$raw, "\n")
    raw_counts <- Read10X(data.dir = paths$raw)

    cat(sprintf("  Filtered: %d genes x %d cells\n",
                nrow(filtered_counts), ncol(filtered_counts)))
    cat(sprintf("  Raw: %d genes x %d cells\n",
                nrow(raw_counts), ncol(raw_counts)))

    # Apply Scrublet calls produced from the original, uncorrected filtered matrix.
    doublet_file <- file.path(DOUBLET_DIR, sprintf("scrublet_sample_%d.csv", sid))
    if (!file.exists(doublet_file)) {
      stop("Missing pre-DecontX doublet calls: ", doublet_file,
           ". Run 00a_doublet_detection_per_sample.py first.")
    }
    doublets <- read.csv(doublet_file, stringsAsFactors = FALSE)
    required_cols <- c("barcode", "scrublet_score", "scrublet_doublet")
    if (!all(required_cols %in% colnames(doublets))) {
      stop("Doublet file lacks required columns: ", doublet_file)
    }
    if (anyDuplicated(doublets$barcode)) stop("Duplicate barcodes in ", doublet_file)
    if (!setequal(doublets$barcode, colnames(filtered_counts))) {
      stop("Doublet-call barcodes do not exactly match the filtered matrix for sample ", sid)
    }
    doublets <- doublets[match(colnames(filtered_counts), doublets$barcode), ]
    call_text <- tolower(as.character(doublets$scrublet_doublet))
    if (!all(call_text %in% c("true", "false"))) {
      stop("Invalid Scrublet calls in ", doublet_file)
    }
    calls <- call_text == "true"
    if (anyNA(calls)) stop("NA doublet calls in ", doublet_file)
    keep_singlet <- !calls
    cat(sprintf("  Pre-DecontX Scrublet: removing %d/%d predicted doublets (%.2f%%)\n",
                sum(calls), length(calls), 100 * mean(calls)))
    filtered_counts <- filtered_counts[, keep_singlet, drop = FALSE]
    doublets <- doublets[keep_singlet, , drop = FALSE]

    # Subset raw to same genes as filtered (DecontX requires matching rows)
    shared_genes <- intersect(rownames(filtered_counts), rownames(raw_counts))
    filtered_counts <- filtered_counts[shared_genes, ]
    raw_counts <- raw_counts[shared_genes, ]
    cat(sprintf("  Shared genes: %d\n", length(shared_genes)))

    # Run DecontX
    cat("  Running DecontX...\n")
    decon_result <- decontX(x          = filtered_counts,
                            background = raw_counts,
                            verbose    = FALSE)

    # Extract results (handle both SCE and list return types)
    if (is(decon_result, "SingleCellExperiment")) {
      decon_counts <- decontXcounts(decon_result)
      contamination <- colData(decon_result)$decontX_contamination
    } else {
      # decontX returns a list when given matrix input
      decon_counts <- decon_result$decontXcounts
      contamination <- decon_result$contamination
    }

    cat(sprintf("  Mean contamination: %.4f\n", mean(contamination)))
    cat(sprintf("  Cells above threshold (%.1f): %d / %d (%.1f%%)\n",
                CONTAM_THRESHOLD,
                sum(contamination > CONTAM_THRESHOLD),
                length(contamination),
                100 * sum(contamination > CONTAM_THRESHOLD) / length(contamination)))

    # Save per-sample results
    result <- list(
      counts        = decon_counts,
      contamination = contamination,
      scrublet_score = setNames(doublets$scrublet_score, doublets$barcode),
      doublet_removed_before_decontx = TRUE,
      sample_id     = sid,
      condition     = cond
    )
    out_file <- file.path(out_dir, sprintf("decontx_sample_%d.rds", sid))
    saveRDS(result, out_file)
    cat("  Saved:", out_file, "\n")

    # Collect contamination scores
    score_df <- data.frame(
      barcode       = colnames(filtered_counts),
      sample_id     = sid,
      condition     = cond,
      contamination = contamination,
      scrublet_score = doublets$scrublet_score,
      stringsAsFactors = FALSE
    )
    all_scores[[i]] <- score_df

  }, error = function(e) {
    stop(sprintf("ERROR processing sample %d (%s): %s", sid, cond, e$message),
         call. = FALSE)
  })

  cat("\n")
}

# ── Combine and save all contamination scores ─────────────────────────────────
cat("Combining contamination scores...\n")
if (length(all_scores) != nrow(sample_info) || any(vapply(all_scores, is.null, logical(1)))) {
  stop("DecontX did not complete for all 16 samples; refusing to write partial output.")
}
scores_combined <- do.call(rbind, all_scores)
scores_file <- file.path(out_dir, "decontx_all_scores.csv")
write.csv(scores_combined, scores_file, row.names = FALSE)
cat("Saved:", scores_file, "\n\n")

# ── Summary ────────────────────────────────────────────────────────────────────
cat("========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(sprintf("Total cells processed: %d\n", nrow(scores_combined)))
cat(sprintf("Overall mean contamination: %.4f\n", mean(scores_combined$contamination)))
cat("\nPer-sample summary:\n")
for (sid in sort(unique(scores_combined$sample_id))) {
  sub <- scores_combined[scores_combined$sample_id == sid, ]
  cat(sprintf("  Sample %2d (%s): %d cells, mean contam = %.4f\n",
              sid, sub$condition[1], nrow(sub), mean(sub$contamination)))
}

cat("\nDone.\n")
