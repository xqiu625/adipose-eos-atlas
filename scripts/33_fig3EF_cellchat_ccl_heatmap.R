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
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig3/fig3E_ccl_eos_to_mac.pdf
#         revision/figures/fig3/fig3F_ccl_mac_to_eos.pdf
# SLURM:  batch, 32G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
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

# --- Load CellChat objects ----------------------------------------------------
cat("Loading CellChat objects...\n")
cellchat_objects <- list()
for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (file.exists(f)) {
    cellchat_objects[[cond]] <- readRDS(f)
    cat(sprintf("  Loaded %s\n", cond))
  }
}
cat(sprintf("Loaded %d conditions\n", length(cellchat_objects)))

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

# --- Extract CCL communications from all conditions ---------------------------
cat("\nExtracting CCL communications...\n")
all_comm <- list()
for (cond in names(cellchat_objects)) {
  cc <- cellchat_objects[[cond]]
  ccl_comm <- tryCatch(
    subsetCommunication(cc, signaling = "CCL"),
    error = function(e) { cat(sprintf("  %s: %s\n", cond, e$message)); NULL }
  )
  if (!is.null(ccl_comm) && nrow(ccl_comm) > 0) {
    ccl_comm$condition <- cond
    all_comm <- c(all_comm, list(ccl_comm))
    cat(sprintf("  %s: %d CCL interactions\n", cond, nrow(ccl_comm)))
  }
}
comm_df <- do.call(rbind, all_comm)
cat(sprintf("Total CCL interactions across conditions: %d\n", nrow(comm_df)))

# Standardize ligand names (Title case: Ccl6, Ccl3, etc.)
comm_df$ligand_clean <- paste0(toupper(substr(comm_df$ligand, 1, 1)),
                               tolower(substr(comm_df$ligand, 2,
                                              nchar(comm_df$ligand))))

cat("Unique ligands:", paste(sort(unique(comm_df$ligand_clean)), collapse = ", "), "\n")
cat("Unique sources:", paste(sort(unique(comm_df$source)), collapse = ", "), "\n")
cat("Unique targets:", paste(sort(unique(comm_df$target)), collapse = ", "), "\n")

# =============================================================================
# Helper: build directional heatmap matrix
# =============================================================================
build_directional_matrix <- function(comm_df, sender_types, receiver_types) {
  # Filter to desired direction
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
