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
# Fig 3A: CCL Chord Diagram — ALL conditions aggregated into ONE diagram
# ==============================================================================
cat("\nGenerating Fig 3A: CCL Chord Diagram (all conditions aggregated)...\n")

# Sum probability matrices across all conditions
unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
agg_total <- matrix(0, nrow = 4, ncol = 4, dimnames = list(unique_groups, unique_groups))

for (cond in condition_order) {
  cc <- cellchat_objects[[cond]]
  if (is.null(cc) || !"CCL" %in% cc@netP$pathways) next

  prob_mat <- cc@netP$prob[,, "CCL"]
  groups <- chord_group(rownames(prob_mat))
  keep <- !is.na(groups)
  prob_sub <- prob_mat[keep, keep]
  group_labels <- groups[keep]

  for (i in seq_along(group_labels)) {
    for (j in seq_along(group_labels)) {
      agg_total[group_labels[i], group_labels[j]] <-
        agg_total[group_labels[i], group_labels[j]] + prob_sub[i, j]
    }
  }
}
# Average across conditions
agg_total <- agg_total / length(cellchat_objects)

pdf(file.path(fig_dir, "fig3A_ccl_chord_diagram.pdf"), width = 9, height = 9)
circos.clear()
circos.par(start.degree = 90, gap.after = c(rep(2, 3), 10))
chordDiagram(agg_total,
             grid.col = group_colors_chord,
             transparency = 0.25,
             directional = 1,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.15))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.1)
}, bg.border = NA)
title(main = "CCL Signaling (All Conditions)", cex.main = 1.4)
circos.clear()
dev.off()
cat("  Saved fig3A_ccl_chord_diagram.pdf\n")

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

# Use CellChat's native netAnalysis_signalingRole_network for the correct format
# Panel 1: Aggregated groups (Eos, M1, M2_Func, M2_Surf)
# Use baseline condition (WT_Female_HFD) as representative
cc_base <- cellchat_objects[[baseline_cond]]

pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles.pdf"), width = 8, height = 5)
tryCatch({
  netAnalysis_signalingRole_network(cc_base, signaling = "CCL",
                                    width = 8, height = 4, font.size = 10)
}, error = function(e) {
  cat(sprintf("  Native plot failed: %s\n  Falling back to manual...\n", e$message))
  # Manual fallback using centrality data
  if (!is.null(cc_base@netP$centr$CCL)) {
    centr <- cc_base@netP$centr$CCL
    role_names <- c("outdeg" = "Sender", "indeg" = "Receiver",
                    "flowbet" = "Mediator", "info" = "Influencer")
    ct <- levels(cc_base@idents)
    mat <- matrix(0, nrow = 4, ncol = length(ct), dimnames = list(
      c("Sender", "Receiver", "Mediator", "Influencer"), ct))
    for (r in names(role_names)) {
      if (!is.null(centr[[r]])) {
        v <- centr[[r]]
        if (is.matrix(v)) v <- rowSums(v)
        mat[role_names[r], names(v)] <- v
      }
    }
    # Aggregate eos subclusters
    eos_ct <- ct[grepl("^Eos_", ct)]
    mac_ct <- ct[grepl("^Mac_", ct)]
    agg_mat <- cbind(
      Eosinophils = rowSums(mat[, eos_ct, drop = FALSE]),
      M1_Mac = if ("Mac_M1" %in% ct) mat[, "Mac_M1"] else rep(0, 4),
      M2_Functional = if ("Mac_M2_Func" %in% ct) mat[, "Mac_M2_Func"] else rep(0, 4),
      M2_Surface = if ("Mac_M2_Surf" %in% ct) mat[, "Mac_M2_Surf"] else rep(0, 4)
    )
    # Color bar annotation
    ann_col <- data.frame(
      CellType = c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface"),
      row.names = colnames(agg_mat))
    ann_colors_ct <- list(CellType = c(
      "Eosinophils" = "#777B7E", "M1_Mac" = "#cb416b",
      "M2_Functional" = "#8FC50E", "M2_Surface" = "#5633AA"))
    pheatmap(agg_mat,
             color = colorRampPalette(c("#F7FCF5", "#74C476", "#006D2C", "#00441B"))(100),
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_col = ann_col, annotation_colors = ann_colors_ct,
             fontsize_row = 12, fontsize_col = 11,
             main = "CCL signaling pathway network",
             border_color = NA)
  }
})
dev.off()
cat("  Saved fig3C_ccl_signaling_roles.pdf\n")

# Panel 2: C3 focus (Eos_C3 vs Mac subtypes)
pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles_C3_focus.pdf"), width = 8, height = 5)
tryCatch({
  # Build matrix with C3 separated
  if (!is.null(cc_base@netP$centr$CCL)) {
    centr <- cc_base@netP$centr$CCL
    role_names <- c("outdeg" = "Sender", "indeg" = "Receiver",
                    "flowbet" = "Mediator", "info" = "Influencer")
    ct <- levels(cc_base@idents)
    mat <- matrix(0, nrow = 4, ncol = length(ct), dimnames = list(
      c("Sender", "Receiver", "Mediator", "Influencer"), ct))
    for (r in names(role_names)) {
      if (!is.null(centr[[r]])) {
        v <- centr[[r]]
        if (is.matrix(v)) v <- rowSums(v)
        mat[role_names[r], names(v)] <- v
      }
    }
    focus_mat <- cbind(
      Eos_C3 = if ("Eos_C3" %in% ct) mat[, "Eos_C3"] else rep(0, 4),
      M1 = if ("Mac_M1" %in% ct) mat[, "Mac_M1"] else rep(0, 4),
      "M2 Funct." = if ("Mac_M2_Func" %in% ct) mat[, "Mac_M2_Func"] else rep(0, 4),
      "M2 Surf." = if ("Mac_M2_Surf" %in% ct) mat[, "Mac_M2_Surf"] else rep(0, 4)
    )
    ann_col <- data.frame(
      CellType = c("Eos_C3", "M1", "M2_Functional", "M2_Surface"),
      row.names = colnames(focus_mat))
    ann_colors_ct <- list(CellType = c(
      "Eos_C3" = "#CBAACB", "M1" = "#cb416b",
      "M2_Functional" = "#8FC50E", "M2_Surface" = "#5633AA"))
    pheatmap(focus_mat,
             color = colorRampPalette(c("#F7FCF5", "#74C476", "#006D2C", "#00441B"))(100),
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_col = ann_col, annotation_colors = ann_colors_ct,
             fontsize_row = 12, fontsize_col = 11,
             main = "CCL signaling pathway network (Eos_C3 Focus)",
             border_color = NA)
  }
}, error = function(e) cat(sprintf("  C3 focus failed: %s\n", e$message)))
dev.off()
cat("  Saved fig3C_ccl_signaling_roles_C3_focus.pdf\n")

# Save CSV of all roles
write.csv(all_roles, file.path(fig_dir, "ccl_signaling_roles.csv"), row.names = FALSE)

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

  # Use CellChat's signaling strength (normalized, not raw prob sums)
  # Extract from netP$centr — outdeg_all/indeg_all for each pathway
  out_list2 <- list()
  in_list2 <- list()

  for (cond in condition_order) {
    cc <- cellchat_objects[[cond]]
    if (is.null(cc)) next

    cell_types_cc <- levels(cc@idents)
    groups_cc <- chord_group(cell_types_cc)

    out_vals <- matrix(0, nrow = length(selected_pathways), ncol = length(broad_groups),
                       dimnames = list(selected_pathways, broad_groups))
    in_vals <- out_vals

    for (pw in selected_pathways) {
      if (pw %in% names(cc@netP$centr)) {
        centr_pw <- cc@netP$centr[[pw]]
        # outdeg = sender strength per cell type
        if (!is.null(centr_pw$outdeg)) {
          od <- centr_pw$outdeg
          if (is.matrix(od)) od <- rowSums(od)
          for (g in broad_groups) {
            g_ct <- cell_types_cc[groups_cc == g]
            g_ct <- g_ct[g_ct %in% names(od)]
            if (length(g_ct) > 0) out_vals[pw, g] <- sum(od[g_ct])
          }
        }
        # indeg = receiver strength per cell type
        if (!is.null(centr_pw$indeg)) {
          id <- centr_pw$indeg
          if (is.matrix(id)) id <- rowSums(id)
          for (g in broad_groups) {
            g_ct <- cell_types_cc[groups_cc == g]
            g_ct <- g_ct[g_ct %in% names(id)]
            if (length(g_ct) > 0) in_vals[pw, g] <- sum(id[g_ct])
          }
        }
      }
    }
    out_list2[[cond]] <- out_vals
    in_list2[[cond]] <- in_vals
  }

  base_out <- out_list2[[baseline_cond]]
  base_in <- in_list2[[baseline_cond]]

  # Plot: 6 small panels in a row (baseline + 5 differences)
  for (direction in c("Outgoing", "Incoming")) {
    sig_list <- if (direction == "Outgoing") out_list2 else in_list2
    sig_base <- if (direction == "Outgoing") base_out else base_in

    pdf(file.path(fig_dir, paste0("fig3B_", tolower(direction), "_heatmap.pdf")),
        width = length(condition_order) * 3.5 + 2, height = 5)

    plot_list <- list()

    # Baseline panel (sequential colorscale)
    base_mat <- sig_list[[baseline_cond]]
    colnames(base_mat) <- c("Eos", "M1", "M2\nFunct.", "M2\nSurf.")
    max_base <- max(base_mat, na.rm = TRUE)
    p_base <- pheatmap(base_mat,
             color = colorRampPalette(c("#F7FCF5", "#74C476", "#00441B", "#FFFF00"))(100),
             breaks = seq(0, max(1.3, max_base * 1.1), length.out = 101),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize_number = 8, fontsize_row = 10, fontsize_col = 9,
             main = paste0("WT Female HFD\n(Baseline)"),
             border_color = "grey70", silent = TRUE)

    # Difference panels (diverging)
    diff_conds <- setdiff(condition_order, baseline_cond)
    diff_plots <- list()
    all_diffs <- c()
    for (cond in diff_conds) {
      d <- sig_list[[cond]] - sig_base
      all_diffs <- c(all_diffs, as.vector(d))
    }
    max_diff <- max(abs(all_diffs), na.rm = TRUE)
    if (max_diff == 0) max_diff <- 0.5

    # Draw all panels together
    library(gridExtra)
    grob_list <- list(p_base[[4]])

    for (cond in diff_conds) {
      d <- sig_list[[cond]] - sig_base
      colnames(d) <- c("Eos", "M1", "M2\nFunct.", "M2\nSurf.")
      p_diff <- pheatmap(d,
               color = colorRampPalette(c("#2166AC", "#92C5DE", "white",
                                          "#F4A582", "#B2182B"))(100),
               breaks = seq(-max_diff, max_diff, length.out = 101),
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.2f",
               fontsize_number = 8, fontsize_row = 10, fontsize_col = 9,
               main = paste0(cond_display[cond], "\nvs Baseline"),
               border_color = "grey70", silent = TRUE)
      grob_list <- c(grob_list, list(p_diff[[4]]))
    }

    grid.arrange(grobs = grob_list, nrow = 1,
                 top = paste(direction, "Signal"))
    dev.off()
  }
  cat("  Saved fig3B_outgoing_heatmap.pdf and fig3B_incoming_heatmap.pdf\n")
}

cat("\n10b CellChat figures complete.\n")
cat("Output:", fig_dir, "\n")
