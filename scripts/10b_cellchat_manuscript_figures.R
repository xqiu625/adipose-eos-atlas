#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 10b_cellchat_manuscript_figures.R
# Fig 3A: CCL chord diagrams for ALL conditions
# Fig 3C: CCL signaling roles (sender/receiver/mediator/influencer) per subcluster
# Fig 3B: Difference heatmaps (outgoing/incoming)
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(pheatmap)
  library(circlize)
  library(RColorBrewer)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
fig_dir <- file.path(BASE, "cellchat", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Female_HFD", "KO_Male_HFD",
                     "WT_Female_HFDInfected", "KO_Female_HFDInfected")
cond_display <- c("WT_Female_HFD" = "WT \u2640", "WT_Male_HFD" = "WT \u2642",
                  "KO_Female_HFD" = "KO \u2640", "KO_Male_HFD" = "KO \u2642",
                  "WT_Female_HFDInfected" = "WT \u2640 Inf",
                  "KO_Female_HFDInfected" = "KO \u2640 Inf")

# ==============================================================================
# 1. Load CellChat objects
# ==============================================================================
cat("Loading CellChat objects...\n")
cellchat_objects <- list()
for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (file.exists(f)) {
    cat(sprintf("  Loading %s\n", cond))
    cc <- readRDS(f)
    cc <- netAnalysis_computeCentrality(cc)
    cellchat_objects[[cond]] <- cc
  }
}
cat(sprintf("Loaded %d conditions\n", length(cellchat_objects)))

# Cell type grouping
chord_group <- function(x) {
  dplyr::case_when(
    grepl("^Eos_", x) ~ "Eosinophils",
    x == "Mac_M1"     ~ "M1_Mac",
    x == "Mac_M2_Func"~ "M2_Functional",
    x == "Mac_M2_Surf"~ "M2_Surface",
    TRUE              ~ NA_character_
  )
}
group_colors_chord <- c("Eosinophils" = "#777B7E", "M1_Mac" = "#cb416b",
                        "M2_Functional" = "#8FC50E", "M2_Surface" = "#5633AA")

# ==============================================================================
# Fig 3A: CCL Chord Diagrams — ALL 6 conditions
# ==============================================================================
cat("\nGenerating Fig 3A: CCL Chord Diagrams (all conditions)...\n")

pdf(file.path(fig_dir, "fig3A_ccl_chord_all_conditions.pdf"), width = 18, height = 12)
par(mfrow = c(2, 3), mar = c(1, 1, 3, 1))

for (cond in condition_order) {
  cc <- cellchat_objects[[cond]]
  if (is.null(cc)) next
  if (!"CCL" %in% cc@netP$pathways) {
    plot.new(); title(main = paste("No CCL -", cond_display[cond]))
    next
  }

  tryCatch({
    prob_mat <- cc@netP$prob[,, "CCL"]
    groups <- chord_group(rownames(prob_mat))
    keep <- !is.na(groups)
    prob_sub <- prob_mat[keep, keep]
    group_labels <- groups[keep]

    unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
    agg_mat <- matrix(0, nrow = 4, ncol = 4, dimnames = list(unique_groups, unique_groups))
    for (i in seq_along(group_labels)) {
      for (j in seq_along(group_labels)) {
        agg_mat[group_labels[i], group_labels[j]] <-
          agg_mat[group_labels[i], group_labels[j]] + prob_sub[i, j]
      }
    }

    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, 3), 10))
    chordDiagram(agg_mat,
                 grid.col = group_colors_chord,
                 transparency = 0.25,
                 directional = 1,
                 direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow",
                 annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.15))
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.9)
    }, bg.border = NA)
    title(main = paste("CCL -", cond_display[cond]), cex.main = 1.4)
    circos.clear()
  }, error = function(e) {
    plot.new(); title(main = paste("Error -", cond))
    cat(sprintf("  Warning: %s chord failed: %s\n", cond, e$message))
  })
}
dev.off()
cat("  Saved fig3A_ccl_chord_all_conditions.pdf\n")

# Also save individual chord diagrams
for (cond in condition_order) {
  cc <- cellchat_objects[[cond]]
  if (is.null(cc) || !"CCL" %in% cc@netP$pathways) next

  tryCatch({
    prob_mat <- cc@netP$prob[,, "CCL"]
    groups <- chord_group(rownames(prob_mat))
    keep <- !is.na(groups)
    prob_sub <- prob_mat[keep, keep]
    group_labels <- groups[keep]

    unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
    agg_mat <- matrix(0, nrow = 4, ncol = 4, dimnames = list(unique_groups, unique_groups))
    for (i in seq_along(group_labels)) {
      for (j in seq_along(group_labels)) {
        agg_mat[group_labels[i], group_labels[j]] <-
          agg_mat[group_labels[i], group_labels[j]] + prob_sub[i, j]
      }
    }

    pdf(file.path(fig_dir, paste0("fig3A_ccl_chord_", cond, ".pdf")), width = 9, height = 9)
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, 3), 10))
    chordDiagram(agg_mat, grid.col = group_colors_chord, transparency = 0.25,
                 directional = 1, direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow", annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.15))
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.1)
    }, bg.border = NA)
    title(main = paste("CCL Signaling -", cond_display[cond]), cex.main = 1.4)
    circos.clear()
    dev.off()
  }, error = function(e) cat(sprintf("  Warning: %s individual chord failed\n", cond)))
}
cat("  Saved individual chord diagrams\n")

# ==============================================================================
# Fig 3C: CCL Signaling Roles — per subcluster (sender/receiver/mediator/influencer)
# ==============================================================================
cat("\nGenerating Fig 3C: CCL Signaling Roles (per subcluster)...\n")

# Aggregate centrality across conditions
all_roles <- data.frame()
for (cond in condition_order) {
  cc <- cellchat_objects[[cond]]
  if (is.null(cc)) next
  if (!"CCL" %in% cc@netP$pathways) next

  centr <- cc@netP$centr
  if (is.null(centr) || !"CCL" %in% names(centr)) {
    cat(sprintf("  WARNING: No centrality for CCL in %s\n", cond))
    next
  }
  ccl_centr <- centr$CCL

  for (role in c("outdeg", "indeg", "flowbet", "info")) {
    if (!is.null(ccl_centr[[role]])) {
      role_name <- switch(role,
        "outdeg" = "Sender", "indeg" = "Receiver",
        "flowbet" = "Mediator", "info" = "Influencer")
      vals <- ccl_centr[[role]]
      # Handle both named vector and matrix formats
      if (is.matrix(vals)) {
        vals <- rowSums(vals)
      }
      if (length(vals) > 0 && !is.null(names(vals))) {
        df <- data.frame(
          cell_type = names(vals), value = as.numeric(vals),
          role = role_name, condition = cond,
          stringsAsFactors = FALSE
        )
        all_roles <- rbind(all_roles, df)
      }
    }
  }
}

if (nrow(all_roles) > 0) {
  cat(sprintf("  Collected %d role entries\n", nrow(all_roles)))

  # Average across conditions
  avg_roles <- all_roles %>%
    group_by(cell_type, role) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  # Panel 1: All cell types aggregated
  role_mat <- tidyr::pivot_wider(avg_roles, names_from = role, values_from = mean_value,
                                  values_fill = 0)
  role_mat_df <- as.data.frame(role_mat)
  rownames(role_mat_df) <- role_mat_df$cell_type
  role_mat_df$cell_type <- NULL
  # Reorder columns
  role_cols <- c("Sender", "Receiver", "Mediator", "Influencer")
  role_cols <- role_cols[role_cols %in% colnames(role_mat_df)]
  role_mat_df <- role_mat_df[, role_cols, drop = FALSE]

  pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles.pdf"), width = 8, height = 6)
  pheatmap(as.matrix(role_mat_df),
           color = colorRampPalette(c("white", "#2D6A4F", "#6A0DAD"))(100),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.3f",
           fontsize_number = 8, fontsize_row = 10, fontsize_col = 11,
           main = "CCL Signaling Roles (All Subtypes)",
           border_color = NA)
  dev.off()
  cat("  Saved fig3C_ccl_signaling_roles.pdf\n")

  # Panel 2: C3 focus — compare C3 vs aggregated Mac subtypes
  c3_roles <- avg_roles %>% filter(cell_type == "Eos_C3")
  mac_roles <- avg_roles %>%
    filter(grepl("^Mac_", cell_type)) %>%
    group_by(cell_type, role) %>%
    summarise(mean_value = mean(mean_value), .groups = "drop")

  focus_roles <- rbind(
    c3_roles %>% mutate(cell_type = "Eos_C3"),
    mac_roles
  )
  focus_mat <- tidyr::pivot_wider(focus_roles, names_from = role, values_from = mean_value,
                                   values_fill = 0)
  focus_df <- as.data.frame(focus_mat)
  rownames(focus_df) <- focus_df$cell_type
  focus_df$cell_type <- NULL
  focus_df <- focus_df[, role_cols[role_cols %in% colnames(focus_df)], drop = FALSE]

  pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles_C3_focus.pdf"), width = 8, height = 5)
  pheatmap(as.matrix(focus_df),
           color = colorRampPalette(c("white", "#2D6A4F", "#6A0DAD"))(100),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.3f",
           fontsize_number = 9, fontsize_row = 11, fontsize_col = 12,
           main = "CCL Signaling Roles (Eos_C3 Focus)",
           border_color = NA)
  dev.off()
  cat("  Saved fig3C_ccl_signaling_roles_C3_focus.pdf\n")

  # Save CSV
  write.csv(avg_roles, file.path(fig_dir, "ccl_signaling_roles.csv"), row.names = FALSE)
} else {
  cat("  WARNING: No centrality data found for CCL pathway\n")
}

# ==============================================================================
# Fig 3B: Difference heatmaps (outgoing/incoming, top 5 pathways)
# ==============================================================================
cat("\nGenerating Fig 3B: Difference heatmaps...\n")

selected_pathways <- c("APP", "CCL", "FN1", "CD45", "THBS")
baseline_cond <- "WT_Female_HFD"
cc_base <- cellchat_objects[[baseline_cond]]

if (!is.null(cc_base)) {
  # Get cell type order
  cell_types <- levels(cc_base@idents)
  broad_groups <- unique(chord_group(cell_types))
  broad_groups <- broad_groups[!is.na(broad_groups)]

  # Compute pathway-level signaling strength per condition
  out_list <- list()
  in_list <- list()

  for (cond in condition_order) {
    cc <- cellchat_objects[[cond]]
    if (is.null(cc)) next

    out_vals <- matrix(0, nrow = length(selected_pathways), ncol = length(broad_groups),
                       dimnames = list(selected_pathways, broad_groups))
    in_vals <- out_vals

    for (pw in selected_pathways) {
      if (pw %in% cc@netP$pathways) {
        prob_mat <- cc@netP$prob[,, pw]
        groups <- chord_group(rownames(prob_mat))
        for (g in broad_groups) {
          g_idx <- which(groups == g)
          if (length(g_idx) > 0) {
            out_vals[pw, g] <- sum(prob_mat[g_idx, , drop = FALSE])
            in_vals[pw, g] <- sum(prob_mat[, g_idx, drop = FALSE])
          }
        }
      }
    }
    out_list[[cond]] <- out_vals
    in_list[[cond]] <- in_vals
  }

  # Compute differences from baseline
  base_out <- out_list[[baseline_cond]]
  base_in <- in_list[[baseline_cond]]

  # Separate baseline (sequential colorscale) and differences (diverging)
  diff_conds <- setdiff(condition_order, baseline_cond)
  baseline_label <- cond_display[baseline_cond]
  diff_labels <- sapply(diff_conds, function(x) cond_display[x])

  for (direction in c("Outgoing", "Incoming")) {
    sig_list <- if (direction == "Outgoing") out_list else in_list
    sig_base <- if (direction == "Outgoing") base_out else base_in

    pdf(file.path(fig_dir, paste0("fig3B_", tolower(direction), "_heatmap.pdf")),
        width = length(condition_order) * 3.5, height = 5)

    # Build: baseline (raw) | diff1 | diff2 | ...
    mat_list <- list()
    mat_list[[1]] <- sig_list[[baseline_cond]]
    colnames(mat_list[[1]]) <- paste(baseline_label, "\n", broad_groups)
    for (i in seq_along(diff_conds)) {
      d <- sig_list[[diff_conds[i]]] - sig_base
      colnames(d) <- paste(diff_labels[i], "\n", broad_groups)
      mat_list[[i + 1]] <- d
    }
    combined <- do.call(cbind, mat_list)

    # Two-tone color: baseline uses sequential, diffs use diverging
    # Use single diverging scale but note baseline is raw
    max_val <- max(abs(combined), na.rm = TRUE)
    pheatmap(combined,
             color = colorRampPalette(c("#2166AC", "#4393C3", "#D1E5F0",
                                        "white",
                                        "#FDDBC7", "#D6604D", "#B2182B"))(100),
             breaks = seq(-max_val, max_val, length.out = 101),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize_number = 7, fontsize_row = 10, fontsize_col = 8,
             main = paste(direction, "Signal (WT Female = raw, others = difference)"),
             border_color = "grey80",
             gaps_col = length(broad_groups))  # Gap after baseline

    dev.off()
  }
  cat("  Saved fig3B_outgoing_heatmap.pdf and fig3B_incoming_heatmap.pdf\n")
}

cat("\n10b CellChat figures complete.\n")
cat("Output:", fig_dir, "\n")
