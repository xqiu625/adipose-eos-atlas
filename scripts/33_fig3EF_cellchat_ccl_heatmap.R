#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 33_fig3EF_cellchat_ccl_heatmap.R
# Fig 3E: CCL signaling Eos -> Mac (hierarchically clustered heatmap, warm)
# Fig 3F: CCL signaling Mac -> Eos (hierarchically clustered heatmap, cool)
#
# Rows    = CCL ligands (summed across receptors)
# Columns = sender x receiver pairs (hierarchically clustered)
# Values  = communication probability averaged across 6 conditions
#
# NOTE: Bypasses CellChat/igraph loading (GLIBCXX_3.4.29 not on HPCC).
#       Extracts @net$prob and @LR$LRsig via raw deserialization.
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig3/fig3E_ccl_eos_to_mac.pdf
#         revision/figures/fig3/fig3F_ccl_mac_to_eos.pdf
# SLURM:  batch, 32G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision", "figures", "fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Male_HFD",
                     "KO_Female_HFD", "KO_Female_HFDInfected",
                     "WT_Female_HFDInfected")

# =============================================================================
# 1. Load CellChat objects via raw deserialization (bypass igraph)
# =============================================================================
cat("Loading CellChat objects (raw deserialization, no igraph)...\n")

extract_ccl_from_rds <- function(rds_path) {
  # Raw deserialize to avoid S4 class dispatch -> CellChat -> igraph
  con <- gzfile(rds_path, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  # CellChat stores:
  #   @net$prob   — 3D array [source, target, interaction]
  #   @net$pval   — 3D array [source, target, interaction]
  #   @LR$LRsig   — data.frame with pathway_name, interaction_name, ligand, receptor, etc.
  prob_arr <- cc_raw@net$prob
  lrsig    <- cc_raw@LR$LRsig

  cat(sprintf("    prob array: %s, LRsig: %d rows\n",
              paste(dim(prob_arr), collapse = " x "), nrow(lrsig)))
  cat(sprintf("    pathway_name column values: %s\n",
              paste(unique(lrsig$pathway_name), collapse = ", ")))

  # Filter to CCL pathway
  ccl_idx <- which(lrsig$pathway_name == "CCL")
  if (length(ccl_idx) == 0) return(NULL)

  ccl_lrsig <- lrsig[ccl_idx, ]
  cat(sprintf("    CCL interactions: %d\n", nrow(ccl_lrsig)))

  # The 3rd dimension of prob_arr should match LRsig interaction names
  interaction_names <- dimnames(prob_arr)[[3]]
  source_names <- dimnames(prob_arr)[[1]]
  target_names <- dimnames(prob_arr)[[2]]

  # Build data frame of CCL communications
  rows <- list()
  for (i in seq_len(nrow(ccl_lrsig))) {
    int_name <- ccl_lrsig$interaction_name[i]
    ligand   <- ccl_lrsig$ligand[i]
    receptor <- ccl_lrsig$receptor[i]

    # Find matching dimension in prob array
    int_idx <- which(interaction_names == int_name)
    if (length(int_idx) == 0) {
      # Try interaction_name_2 if available
      if ("interaction_name_2" %in% colnames(ccl_lrsig)) {
        int_name2 <- ccl_lrsig$interaction_name_2[i]
        int_idx <- which(interaction_names == int_name2)
      }
    }
    if (length(int_idx) == 0) next

    # Extract prob slice for this interaction
    prob_slice <- prob_arr[, , int_idx]

    # Find non-zero entries
    nz <- which(prob_slice > 0, arr.ind = TRUE)
    if (nrow(nz) == 0) next

    for (j in seq_len(nrow(nz))) {
      rows <- c(rows, list(data.frame(
        source   = source_names[nz[j, 1]],
        target   = target_names[nz[j, 2]],
        ligand   = ligand,
        receptor = receptor,
        prob     = prob_slice[nz[j, 1], nz[j, 2]],
        interaction_name = int_name,
        stringsAsFactors = FALSE
      )))
    }
  }

  rm(cc_raw); gc(verbose = FALSE)

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

all_comm <- list()
for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (!file.exists(f)) { cat(sprintf("  SKIP %s (not found)\n", cond)); next }
  cat(sprintf("  Loading %s ...\n", cond))
  ccl_df <- tryCatch(extract_ccl_from_rds(f),
                      error = function(e) { cat(sprintf("  ERROR %s: %s\n", cond, e$message)); NULL })
  if (!is.null(ccl_df) && nrow(ccl_df) > 0) {
    ccl_df$condition <- cond
    all_comm <- c(all_comm, list(ccl_df))
    cat(sprintf("    %s: %d CCL pairs\n", cond, nrow(ccl_df)))
  }
}

comm_df <- do.call(rbind, all_comm)
cat(sprintf("\nTotal CCL interactions across conditions: %d\n", nrow(comm_df)))

# Standardize ligand names (Title case: Ccl6, Ccl3, etc.)
comm_df$ligand_clean <- paste0(toupper(substr(comm_df$ligand, 1, 1)),
                               tolower(substr(comm_df$ligand, 2,
                                              nchar(comm_df$ligand))))

cat("Unique ligands:", paste(sort(unique(comm_df$ligand_clean)), collapse = ", "), "\n")
cat("Unique sources:", paste(sort(unique(comm_df$source)), collapse = ", "), "\n")
cat("Unique targets:", paste(sort(unique(comm_df$target)), collapse = ", "), "\n")

# --- Cell type definitions ----------------------------------------------------
eos_types <- paste0("Eos_C", 0:5)
mac_types <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")

eos_short <- setNames(paste0("C", 0:5), eos_types)
mac_short <- c("Mac_M1" = "M1", "Mac_M2_Func" = "M2 Funct.", "Mac_M2_Surf" = "M2 Surf.")

# Annotation colors (from 00_theme.R)
eos_anno_col <- eos_colors[paste0("C", 0:5)]
mac_anno_col <- c("M1" = mac_colors[["Mac_M1"]],
                  "M2 Funct." = mac_colors[["Mac_M2_Func"]],
                  "M2 Surf."  = mac_colors[["Mac_M2_Surf"]])

# =============================================================================
# Helper: build directional heatmap matrix
# =============================================================================
build_directional_matrix <- function(comm_df, sender_types, receiver_types) {
  df <- comm_df %>%
    filter(source %in% sender_types & target %in% receiver_types)

  if (nrow(df) == 0) return(NULL)

  # Sum across receptors within each (ligand, source, target, condition)
  df_agg <- df %>%
    group_by(ligand_clean, source, target, condition) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop")

  # Average across conditions
  df_avg <- df_agg %>%
    group_by(ligand_clean, source, target) %>%
    summarise(prob = mean(prob, na.rm = TRUE), .groups = "drop")

  # Create column IDs
  df_avg$col_id <- paste0(df_avg$source, "|", df_avg$target)

  # All possible sender x receiver combinations
  all_pairs <- expand.grid(source = sender_types, target = receiver_types,
                           stringsAsFactors = FALSE)
  all_pairs$col_id <- paste0(all_pairs$source, "|", all_pairs$target)

  # All ligands with signal
  all_ligands <- sort(unique(df_avg$ligand_clean))

  # Build matrix
  mat <- matrix(0, nrow = length(all_ligands), ncol = nrow(all_pairs),
                dimnames = list(all_ligands, all_pairs$col_id))
  for (i in seq_len(nrow(df_avg))) {
    lig <- df_avg$ligand_clean[i]
    cid <- df_avg$col_id[i]
    if (lig %in% rownames(mat) && cid %in% colnames(mat)) {
      mat[lig, cid] <- df_avg$prob[i]
    }
  }

  # Drop all-zero rows
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  list(mat = mat, pairs = all_pairs)
}

# =============================================================================
# Fig 3E: Eos -> Mac (warm palette)
# =============================================================================
cat("\n--- Fig 3E: Eos Subclusters to Mac Subtypes ---\n")

res_e <- build_directional_matrix(comm_df, eos_types, mac_types)
mat_e <- res_e$mat
pairs_e <- res_e$pairs
cat(sprintf("  Ligands: %d, Columns: %d\n", nrow(mat_e), ncol(mat_e)))
cat("  Ligands:", paste(rownames(mat_e), collapse = ", "), "\n")

# Column annotations
col_anno_e <- data.frame(
  `Eos\nSender`  = eos_short[pairs_e$source],
  `Mac\nReceiver` = mac_short[pairs_e$target],
  row.names = pairs_e$col_id,
  check.names = FALSE
)

anno_colors_e <- list(
  `Eos\nSender`  = eos_anno_col,
  `Mac\nReceiver` = mac_anno_col
)

# Warm color palette
warm_pal <- colorRampPalette(c("#FFFFF0", "#FFFFCC", "#FED976",
                               "#FD8D3C", "#E31A1C", "#800026"))(100)

pdf_e <- file.path(out_dir, "fig3E_ccl_eos_to_mac.pdf")
cairo_pdf(pdf_e, width = 10, height = 5)
pheatmap(mat_e,
         color = warm_pal,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = col_anno_e,
         annotation_colors = anno_colors_e,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 11,
         fontsize = 10,
         border_color = "grey80",
         cellwidth = 18, cellheight = 20,
         main = "Eos Subclusters to Mac Subtypes")
dev.off()
cat(sprintf("  Saved: %s\n", pdf_e))

# =============================================================================
# Fig 3F: Mac -> Eos (cool palette)
# =============================================================================
cat("\n--- Fig 3F: Mac Subtypes to Eos Subclusters ---\n")

res_f <- build_directional_matrix(comm_df, mac_types, eos_types)
mat_f <- res_f$mat
pairs_f <- res_f$pairs
cat(sprintf("  Ligands: %d, Columns: %d\n", nrow(mat_f), ncol(mat_f)))
cat("  Ligands:", paste(rownames(mat_f), collapse = ", "), "\n")

# Column annotations
col_anno_f <- data.frame(
  `Mac\nSender`   = mac_short[pairs_f$source],
  `Eos\nReceiver` = eos_short[pairs_f$target],
  row.names = pairs_f$col_id,
  check.names = FALSE
)

anno_colors_f <- list(
  `Mac\nSender`   = mac_anno_col,
  `Eos\nReceiver` = eos_anno_col
)

# Cool color palette
cool_pal <- colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1",
                               "#4292C6", "#2171B5", "#084594"))(100)

pdf_f <- file.path(out_dir, "fig3F_ccl_mac_to_eos.pdf")
cairo_pdf(pdf_f, width = 10, height = 5)
pheatmap(mat_f,
         color = cool_pal,
         cluster_rows = TRUE, cluster_cols = TRUE,
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = col_anno_f,
         annotation_colors = anno_colors_f,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 11,
         fontsize = 10,
         border_color = "grey80",
         cellwidth = 18, cellheight = 20,
         main = "Mac Subtypes to Eos Subclusters")
dev.off()
cat(sprintf("  Saved: %s\n", pdf_f))

# =============================================================================
# Print summary
# =============================================================================
cat("\n========== Summary ==========\n")
cat("Fig 3E (Eos -> Mac):\n")
cat(sprintf("  %d ligands x %d sender-receiver pairs\n", nrow(mat_e), ncol(mat_e)))
cat(sprintf("  Max prob: %.4f\n", max(mat_e)))
cat("Fig 3F (Mac -> Eos):\n")
cat(sprintf("  %d ligands x %d sender-receiver pairs\n", nrow(mat_f), ncol(mat_f)))
cat(sprintf("  Max prob: %.4f\n", max(mat_f)))
cat(sprintf("\nOutput: %s\n", out_dir))
cat("33 complete.\n")
