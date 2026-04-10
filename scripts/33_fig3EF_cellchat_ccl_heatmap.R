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
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
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
  con <- gzfile(rds_path, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  prob_arr <- cc_raw@net$prob
  lrsig    <- cc_raw@LR$LRsig

  cat(sprintf("    prob: %s, LRsig: %d rows\n",
              paste(dim(prob_arr), collapse="x"), nrow(lrsig)))

  ccl_idx <- which(lrsig$pathway_name == "CCL")
  if (length(ccl_idx) == 0) return(NULL)

  ccl_lrsig <- lrsig[ccl_idx, ]
  interaction_names <- dimnames(prob_arr)[[3]]
  source_names <- dimnames(prob_arr)[[1]]
  target_names <- dimnames(prob_arr)[[2]]

  rows <- list()
  for (i in seq_len(nrow(ccl_lrsig))) {
    int_name <- ccl_lrsig$interaction_name[i]
    ligand   <- ccl_lrsig$ligand[i]
    receptor <- ccl_lrsig$receptor[i]

    int_idx <- which(interaction_names == int_name)
    if (length(int_idx) == 0 && "interaction_name_2" %in% colnames(ccl_lrsig)) {
      int_idx <- which(interaction_names == ccl_lrsig$interaction_name_2[i])
    }
    if (length(int_idx) == 0) next

    prob_slice <- prob_arr[, , int_idx]
    nz <- which(prob_slice > 0, arr.ind = TRUE)
    if (nrow(nz) == 0) next

    for (j in seq_len(nrow(nz))) {
      rows <- c(rows, list(data.frame(
        source = source_names[nz[j, 1]],
        target = target_names[nz[j, 2]],
        ligand = ligand, receptor = receptor,
        prob = prob_slice[nz[j, 1], nz[j, 2]],
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
  if (!file.exists(f)) { cat(sprintf("  SKIP %s\n", cond)); next }
  cat(sprintf("  %s ...\n", cond))
  ccl_df <- tryCatch(extract_ccl_from_rds(f),
                      error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
  if (!is.null(ccl_df) && nrow(ccl_df) > 0) {
    ccl_df$condition <- cond
    all_comm <- c(all_comm, list(ccl_df))
    cat(sprintf("    %d CCL pairs\n", nrow(ccl_df)))
  }
}
comm_df <- do.call(rbind, all_comm)
cat(sprintf("\nTotal CCL interactions: %d\n", nrow(comm_df)))

# Standardize ligand names (Title case)
comm_df$ligand_clean <- paste0(toupper(substr(comm_df$ligand, 1, 1)),
                               tolower(substr(comm_df$ligand, 2,
                                              nchar(comm_df$ligand))))

cat("Ligands:", paste(sort(unique(comm_df$ligand_clean)), collapse=", "), "\n")

# --- Cell type definitions ----------------------------------------------------
eos_types <- paste0("Eos_C", 0:5)
mac_types <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")

eos_short <- setNames(paste0("C", 0:5), eos_types)
mac_short <- c("Mac_M1" = "M1", "Mac_M2_Func" = "M2\nFunct.", "Mac_M2_Surf" = "M2\nSurf.")

# Annotation colors
eos_anno_col <- eos_colors[paste0("C", 0:5)]
mac_anno_col <- c("M1" = unname(mac_colors["Mac_M1"]),
                  "M2\nFunct." = unname(mac_colors["Mac_M2_Func"]),
                  "M2\nSurf."  = unname(mac_colors["Mac_M2_Surf"]))

# =============================================================================
# Helper: build directional matrix
# =============================================================================
build_directional_matrix <- function(comm_df, sender_types, receiver_types) {
  df <- comm_df %>%
    filter(source %in% sender_types & target %in% receiver_types)
  if (nrow(df) == 0) return(NULL)

  df_agg <- df %>%
    group_by(ligand_clean, source, target, condition) %>%
    summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop")

  df_avg <- df_agg %>%
    group_by(ligand_clean, source, target) %>%
    summarise(prob = mean(prob, na.rm = TRUE), .groups = "drop")

  df_avg$col_id <- paste0(df_avg$source, "|", df_avg$target)

  all_pairs <- expand.grid(source = sender_types, target = receiver_types,
                           stringsAsFactors = FALSE)
  all_pairs$col_id <- paste0(all_pairs$source, "|", all_pairs$target)

  all_ligands <- sort(unique(df_avg$ligand_clean))

  mat <- matrix(0, nrow = length(all_ligands), ncol = nrow(all_pairs),
                dimnames = list(all_ligands, all_pairs$col_id))
  for (i in seq_len(nrow(df_avg))) {
    lig <- df_avg$ligand_clean[i]
    cid <- df_avg$col_id[i]
    if (lig %in% rownames(mat) && cid %in% colnames(mat))
      mat[lig, cid] <- df_avg$prob[i]
  }
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  list(mat = mat, pairs = all_pairs)
}

# =============================================================================
# Fig 3E: Eos -> Mac (warm palette)
# =============================================================================
cat("\n--- Fig 3E: Eos -> Mac ---\n")

res_e <- build_directional_matrix(comm_df, eos_types, mac_types)
mat_e <- res_e$mat
pairs_e <- res_e$pairs
cat(sprintf("  %d ligands x %d columns\n", nrow(mat_e), ncol(mat_e)))

# Column annotation vectors
sender_e   <- eos_short[pairs_e$source]
receiver_e <- mac_short[pairs_e$target]
names(sender_e) <- pairs_e$col_id
names(receiver_e) <- pairs_e$col_id

# Bottom annotation: colored blocks with sender/receiver identity
bot_e <- HeatmapAnnotation(
  `Eos\nSender`   = sender_e,
  `Mac\nReceiver`  = receiver_e,
  col = list(`Eos\nSender`  = eos_anno_col,
             `Mac\nReceiver` = mac_anno_col),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(5, "mm"),
  show_legend = TRUE,
  which = "column"
)

# Warm color function
col_e <- colorRamp2(
  seq(0, max(mat_e), length.out = 7),
  c("#FFFFF0", "#FFFFCC", "#FED976", "#FEB24C", "#FD8D3C", "#E31A1C", "#800026")
)

ht_e <- Heatmap(mat_e,
  name = " ",
  col = col_e,
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  bottom_annotation = bot_e,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontface = "italic", fontsize = 11),
  column_title = "Eos Subclusters to Mac Subtypes",
  column_title_gp = gpar(fontface = "bold", fontsize = 13),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  heatmap_legend_param = list(
    at = pretty(c(0, max(mat_e)), n = 4),
    border = NA
  ),
  width = unit(ncol(mat_e) * 6, "mm"),
  height = unit(nrow(mat_e) * 8, "mm")
)

pdf_e <- file.path(out_dir, "fig3E_ccl_eos_to_mac.pdf")
cairo_pdf(pdf_e, width = 9, height = 4)
draw(ht_e, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()
cat(sprintf("  Saved: %s\n", pdf_e))

# =============================================================================
# Fig 3F: Mac -> Eos (cool palette)
# =============================================================================
cat("\n--- Fig 3F: Mac -> Eos ---\n")

res_f <- build_directional_matrix(comm_df, mac_types, eos_types)
mat_f <- res_f$mat
pairs_f <- res_f$pairs
cat(sprintf("  %d ligands x %d columns\n", nrow(mat_f), ncol(mat_f)))

sender_f   <- mac_short[pairs_f$source]
receiver_f <- eos_short[pairs_f$target]
names(sender_f) <- pairs_f$col_id
names(receiver_f) <- pairs_f$col_id

bot_f <- HeatmapAnnotation(
  `Mac\nSender`   = sender_f,
  `Eos\nReceiver`  = receiver_f,
  col = list(`Mac\nSender`  = mac_anno_col,
             `Eos\nReceiver` = eos_anno_col),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(5, "mm"),
  show_legend = TRUE,
  which = "column"
)

col_f <- colorRamp2(
  seq(0, max(mat_f), length.out = 7),
  c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#4292C6", "#2171B5", "#084594")
)

ht_f <- Heatmap(mat_f,
  name = " ",
  col = col_f,
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  bottom_annotation = bot_f,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontface = "italic", fontsize = 11),
  column_title = "Mac Subtypes to Eos Subclusters",
  column_title_gp = gpar(fontface = "bold", fontsize = 13),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  heatmap_legend_param = list(
    at = pretty(c(0, max(mat_f)), n = 4),
    border = NA
  ),
  width = unit(ncol(mat_f) * 6, "mm"),
  height = unit(nrow(mat_f) * 8, "mm")
)

pdf_f <- file.path(out_dir, "fig3F_ccl_mac_to_eos.pdf")
cairo_pdf(pdf_f, width = 9, height = 5.5)
draw(ht_f, padding = unit(c(2, 15, 2, 2), "mm"))
dev.off()
cat(sprintf("  Saved: %s\n", pdf_f))

# =============================================================================
# Summary
# =============================================================================
cat("\n========== Summary ==========\n")
cat(sprintf("Fig 3E: %d ligands x %d pairs (max=%.4f)\n",
            nrow(mat_e), ncol(mat_e), max(mat_e)))
cat(sprintf("  Ligands: %s\n", paste(rownames(mat_e), collapse=", ")))
cat(sprintf("Fig 3F: %d ligands x %d pairs (max=%.4f)\n",
            nrow(mat_f), ncol(mat_f), max(mat_f)))
cat(sprintf("  Ligands: %s\n", paste(rownames(mat_f), collapse=", ")))
cat(sprintf("\nOutput: %s\n", out_dir))
cat("33 complete.\n")
