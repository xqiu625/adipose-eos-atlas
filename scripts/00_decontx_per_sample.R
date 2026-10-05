#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
# Prefer the isolated SingleCellExperiment/Bioconductor library.
sce_fix_lib <- Sys.getenv("SCE_FIX_LIB", "")
if (nzchar(Sys.getenv("RERUN_R_LIB")) && !nzchar(sce_fix_lib)) {
  stop("SCE_FIX_LIB is required alongside RERUN_R_LIB: without it, SingleCellExperiment/celda may resolve a mismatched Bioconductor-release dependency chain from this cluster's R/4.5.0 site library and object construction fails (RERUN_LOG.md Checkpoints 24 and 26). Refusing to proceed.")
}
.libPaths(c(sce_fix_lib, lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .expected_prefix <- normalizePath(c(sce_fix_lib, Sys.getenv("RERUN_R_LIB")), mustWork = FALSE)
  .resolved_prefix <- normalizePath(.libPaths()[1:2], mustWork = FALSE)
  if (!identical(.resolved_prefix, .expected_prefix)) {
    stop(sprintf(
      "R library audit guard: .libPaths()[1:2] resolved to (%s) but expected (%s) -- the shared/legacy or site library would take precedence over the pinned isolated libraries. Refusing to proceed.",
      paste(.resolved_prefix, collapse = ", "), paste(.expected_prefix, collapse = ", ")))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("SCE_FIX_LIB:", sce_fix_lib, "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

## DecontX ambient RNA removal per sample, after pre-DecontX doublet removal
## SLURM: --mem=128G --cpus-per-task=8 --time=4:00:00 -p batch

# Optional single-sample execution and explicit input manifest.
cli <- commandArgs(trailingOnly = TRUE)
sample_id_arg <- NULL
manifest_path <- NULL
if (length(cli) %% 2L != 0L) stop("Expected --sample-id N and/or --input-manifest PATH")
if (length(cli)) for (j in seq.int(1L, length(cli), by = 2L)) {
  if (cli[j] == "--sample-id") {
    if (!grepl("^([1-9]|1[0-6])$", cli[j + 1L])) stop("sample-id must be 1–16")
    sample_id_arg <- as.integer(cli[j + 1L])
  } else if (cli[j] == "--input-manifest") {
    manifest_path <- cli[j + 1L]
  } else stop("Unknown argument: ", cli[j])
}

# ── Libraries ──────────────────────────────────────────────────────────────────
# Align filtered and background matrices by gene ID.
align_decontx_features <- function(filtered, raw, features) {
  ids <- rownames(filtered)
  raw_ids <- rownames(raw)
  if (is.null(ids) || is.null(raw_ids) || anyNA(ids) || anyNA(raw_ids) ||
      anyDuplicated(ids) || anyDuplicated(raw_ids) ||
      any(!nzchar(ids)) || any(!nzchar(raw_ids)) || ncol(features) < 2L ||
      anyDuplicated(features[[1]]) || anyNA(features[[2]])) {
    stop("Missing/duplicate feature IDs or invalid feature metadata")
  }
  selected <- match(ids, raw_ids)
  annotation <- match(ids, features[[1]])
  if (anyNA(selected) || anyNA(annotation)) stop("Filtered gene ID absent from raw or metadata")
  symbols <- make.unique(as.character(features[[2]][annotation]))
  if (any(!nzchar(symbols))) stop("Empty gene symbol")
  raw <- raw[selected, , drop = FALSE]
  stopifnot(identical(rownames(raw), ids))
  rownames(filtered) <- symbols
  rownames(raw) <- symbols
  list(filtered = filtered, raw = raw, feature_ids = setNames(ids, symbols))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(celda)
  library(Matrix)
})
cat("Analysis packages:", paste(c("Seurat", "celda", "Matrix"),
    vapply(c("Seurat", "celda", "Matrix"), function(p) as.character(packageVersion(p)),
           character(1)), collapse = "; "), "\n")

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- Sys.getenv(
  "EOS_ATLAS_BASE",
  unset = "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
)
CONTAM_THRESHOLD <- 0.2
DOUBLET_DIR <- file.path(BASE, "doublets")

# historical: use raw counts; DecontX removes barcodes shared with singlets.
# exclude_all_original_called_cells: also exclude called doublets from background.
BACKGROUND_POLICY_ALLOWED <- c("historical", "exclude_all_original_called_cells")
BACKGROUND_POLICY <- Sys.getenv("DECONTX_BACKGROUND_POLICY", unset = "historical")
if (!BACKGROUND_POLICY %in% BACKGROUND_POLICY_ALLOWED) {
  stop(sprintf(
    "DECONTX_BACKGROUND_POLICY='%s' is not one of the allowed policies (%s) -- refusing to guess",
    BACKGROUND_POLICY, paste(BACKGROUND_POLICY_ALLOWED, collapse = ", ")
  ))
}
cat("Ambient-background policy:", BACKGROUND_POLICY,
    if (identical(BACKGROUND_POLICY, "historical")) "(official production default, matches methods.md/note.md)"
    else "(candidate doublet-aware branch, NOT yet validated/approved for production)", "\n")

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
input_manifest <- NULL
if (!is.null(manifest_path)) {
  input_manifest <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "pool", "filtered_dir", "raw_dir")
  if (!all(required %in% names(input_manifest)) || nrow(input_manifest) != 16L ||
      anyDuplicated(input_manifest$sample_id) ||
      !setequal(input_manifest$sample_id, 1:16)) stop("Invalid 16-sample input manifest")
  for (pool in names(gem_map)) {
    if (any(input_manifest$pool[match(gem_map[[pool]], input_manifest$sample_id)] != pool)) {
      stop("Incorrect pool mapping: ", pool)
    }
  }
}

# Build per-sample paths
get_sample_paths <- function(sample_id) {
  if (!is.null(input_manifest)) {
    row <- input_manifest[input_manifest$sample_id == sample_id, ]
    return(list(filtered = row$filtered_dir, raw = row$raw_dir))
  }
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
  if (!is.null(sample_id_arg) && sid != sample_id_arg) next
  cond <- sample_info$condition[i]

  cat(sprintf("[Sample %d/%d] %s (ID: %d)\n", i, nrow(sample_info), cond, sid))

  paths <- get_sample_paths(sid)
  out_file <- file.path(out_dir, sprintf("decontx_sample_%d.rds", sid))
  if (!is.null(sample_id_arg) && file.exists(out_file)) stop("Refusing to overwrite: ", out_file)

  tryCatch({
    # Load filtered counts
    cat("  Loading filtered counts:", paths$filtered, "\n")
    filtered_counts <- Read10X(data.dir = paths$filtered, gene.column = 1)

    # Load raw counts (background)
    cat("  Loading raw counts:", paths$raw, "\n")
    raw_counts <- Read10X(data.dir = paths$raw, gene.column = 1)

    cat(sprintf("  Filtered: %d genes x %d cells\n",
                nrow(filtered_counts), ncol(filtered_counts)))
    cat(sprintf("  Raw: %d genes x %d cells\n",
                nrow(raw_counts), ncol(raw_counts)))

    # Raw and filtered have different gene universes. Independent make.unique()
    # on symbols paired one different Ensembl ID in all 16 historical samples.
    # Align IDs first, then restore filtered symbols for downstream compatibility.
    feature_table <- read.delim(gzfile(file.path(paths$filtered, "features.tsv.gz")),
                               header = FALSE, stringsAsFactors = FALSE)
    aligned <- align_decontx_features(filtered_counts, raw_counts, feature_table)
    filtered_counts <- aligned$filtered
    raw_counts <- aligned$raw
    feature_ids <- aligned$feature_ids
    rm(aligned)
    cat("  Feature alignment: exact gene IDs; restored filtered gene symbols\n")

    # Record called barcodes before selecting the background.
    original_called_barcodes <- colnames(filtered_counts)
    raw_barcodes <- colnames(raw_counts)
    if (is.null(original_called_barcodes) || is.null(raw_barcodes) ||
        anyNA(original_called_barcodes) || anyNA(raw_barcodes) ||
        any(!nzchar(original_called_barcodes)) || any(!nzchar(raw_barcodes)) ||
        anyDuplicated(original_called_barcodes) || anyDuplicated(raw_barcodes) ||
        !all(original_called_barcodes %in% raw_barcodes)) {
      stop("Cannot identify background: called/raw barcodes missing, duplicated, or inconsistent")
    }
    if (identical(BACKGROUND_POLICY, "exclude_all_original_called_cells")) {
      # Exclude all called singlets and doublets from the background.
      raw_counts <- raw_counts[, !raw_barcodes %in% original_called_barcodes, drop = FALSE]
      if (ncol(raw_counts) == 0L || sum(raw_counts) <= 0) {
        stop("No informative non-called background remains; review input, do not silently change model")
      }
    }
    # DecontX removes singlet overlap; called doublets remain in this background.
    background_barcodes <- colnames(raw_counts)
    cat(sprintf("  Ambient background (%s): %d barcodes passed to decontX() (of %d original called cells, %d raw droplets total)\n",
                BACKGROUND_POLICY, length(background_barcodes),
                length(original_called_barcodes), length(raw_barcodes)))

    # Apply Scrublet calls produced from the original, uncorrected filtered matrix.
    doublet_file <- file.path(DOUBLET_DIR, sprintf("scrublet_sample_%d.csv", sid))
    if (!file.exists(doublet_file)) {
      stop("Missing pre-DecontX doublet calls: ", doublet_file,
           ". Run 00_a_doublet_detection_per_sample.py first.")
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

    stopifnot(identical(rownames(filtered_counts), rownames(raw_counts)))
    cat(sprintf("  ID-aligned shared genes: %d\n", nrow(filtered_counts)))

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
    if (!identical(dim(decon_counts), dim(filtered_counts)) ||
        !identical(colnames(decon_counts), colnames(filtered_counts)) ||
        !identical(rownames(decon_counts), rownames(filtered_counts))) {
      stop("DecontX changed matrix dimensions or cell/feature alignment")
    }
    if (length(contamination) != ncol(filtered_counts) ||
        any(!is.finite(contamination)) || any(contamination < 0 | contamination > 1)) {
      stop("Invalid DecontX contamination vector")
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
      feature_ids   = feature_ids,
      feature_alignment = "ensembl_id_then_filtered_symbols",
      background_policy = BACKGROUND_POLICY,
      background_barcodes = background_barcodes,
      original_called_barcodes = original_called_barcodes,
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
      background_policy = BACKGROUND_POLICY,
      contamination = contamination,
      scrublet_score = doublets$scrublet_score,
      stringsAsFactors = FALSE
    )
    all_scores[[i]] <- score_df
    if (!is.null(sample_id_arg)) {
      write.csv(score_df, file.path(out_dir, sprintf("decontx_scores_sample_%d.csv", sid)),
                row.names = FALSE)
    }

  }, error = function(e) {
    stop(sprintf("ERROR processing sample %d (%s): %s", sid, cond, e$message),
         call. = FALSE)
  })

  cat("\n")
}

if (!is.null(sample_id_arg)) {
  cat("Single sample complete; shared summary is withheld until all 16 tasks validate.\n")
  quit(save = "no", status = 0L)
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
