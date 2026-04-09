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

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Male_HFD", "KO_Female_HFD",
                     "KO_Female_HFDInfected", "WT_Female_HFDInfected")
cond_display <- c("WT_Female_HFD" = "WT Female HFD\n(Baseline)",
                  "WT_Male_HFD" = "WT Male HFD",
                  "KO_Male_HFD" = "KO Male HFD",
                  "KO_Female_HFD" = "KO Female HFD",
                  "KO_Female_HFDInfected" = "KO Female HFD\nInfected",
                  "WT_Female_HFDInfected" = "WT Female HFD\nInfected")

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
# Fig 3A: Chord Diagram — ALL pathways, ALL conditions aggregated
# Shows total Eos <-> Mac communication (not just CCL)
# ==============================================================================
cat("\nGenerating Fig 3A: All-pathway Chord Diagram (all conditions aggregated)...\n")

unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
agg_total <- matrix(0, nrow = 4, ncol = 4, dimnames = list(unique_groups, unique_groups))

for (cond in condition_order) {
  cc <- cellchat_objects[[cond]]
  if (is.null(cc)) next

  # Sum across ALL pathways
  for (pw in cc@netP$pathways) {
    prob_mat <- cc@netP$prob[,, pw]
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
}
# Average across conditions
agg_total <- agg_total / length(cellchat_objects)

pdf(file.path(fig_dir, "fig3A_chord_all_pathways.pdf"), width = 9, height = 9)
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
circos.clear()
dev.off()
cat("  Saved fig3A_chord_all_pathways.pdf\n")

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

# Use CellChat's native netAnalysis_signalingRole_network
# This computes sender/receiver/mediator/influencer correctly from the full network
baseline_cond <- "WT_Female_HFD"
cc_base <- cellchat_objects[[baseline_cond]]

# Panel 1: All eos (aggregated) + Mac subtypes
# Create a copy with eos subclusters merged into "Eosinophils"
cat("  Panel 1: Aggregated Eos + Mac subtypes...\n")
cc_agg <- cc_base
idents_agg <- as.character(cc_agg@idents)
idents_agg[grepl("^Eos_", idents_agg)] <- "Eosinophils"
idents_agg[idents_agg == "Mac_M1"] <- "M1_Mac"
idents_agg[idents_agg == "Mac_M2_Func"] <- "M2_Functional"
idents_agg[idents_agg == "Mac_M2_Surf"] <- "M2_Surface"
cc_agg@idents <- factor(idents_agg,
                         levels = c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface"))

# Re-run communication on aggregated identities
tryCatch({
  cc_agg <- computeCommunProb(cc_agg, type = "triMean")
  cc_agg <- filterCommunication(cc_agg, min.cells = 10)
  cc_agg <- computeCommunProbPathway(cc_agg)
  cc_agg <- aggregateNet(cc_agg)
  cc_agg <- netAnalysis_computeCentrality(cc_agg)

  pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles.pdf"), width = 8, height = 5)
  netAnalysis_signalingRole_network(cc_agg, signaling = "CCL",
                                    width = 8, height = 4, font.size = 10)
  dev.off()
  cat("  Saved fig3C_ccl_signaling_roles.pdf\n")
}, error = function(e) cat(sprintf("  Panel 1 failed: %s\n", e$message)))

# Panel 2: Eos_C3 + Mac subtypes
cat("  Panel 2: Eos_C3 + Mac subtypes...\n")
cc_c3 <- cc_base
idents_c3 <- as.character(cc_c3@idents)
# Keep only Eos_C3 and Mac subtypes
keep_cells <- grepl("^Eos_C3$|^Mac_", idents_c3)

tryCatch({
  cc_c3 <- subsetCellChat(cc_c3, idents.use = c("Eos_C3", "Mac_M1", "Mac_M2_Func", "Mac_M2_Surf"))

  # Rename for display
  idents_c3 <- as.character(cc_c3@idents)
  idents_c3[idents_c3 == "Mac_M1"] <- "M1_Mac"
  idents_c3[idents_c3 == "Mac_M2_Func"] <- "M2_Functional"
  idents_c3[idents_c3 == "Mac_M2_Surf"] <- "M2_Surface"
  cc_c3@idents <- factor(idents_c3,
                          levels = c("Eos_C3", "M1_Mac", "M2_Functional", "M2_Surface"))

  cc_c3 <- computeCommunProb(cc_c3, type = "triMean")
  cc_c3 <- filterCommunication(cc_c3, min.cells = 10)
  cc_c3 <- computeCommunProbPathway(cc_c3)
  cc_c3 <- aggregateNet(cc_c3)
  cc_c3 <- netAnalysis_computeCentrality(cc_c3)

  pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles_C3_focus.pdf"), width = 8, height = 5)
  netAnalysis_signalingRole_network(cc_c3, signaling = "CCL",
                                    width = 8, height = 4, font.size = 10)
  dev.off()
  cat("  Saved fig3C_ccl_signaling_roles_C3_focus.pdf\n")
}, error = function(e) cat(sprintf("  Panel 2 failed: %s\n", e$message)))

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

  diff_conds <- setdiff(condition_order, baseline_cond)
  col_labels <- c("Eos", "M1", "M2\nFunct.", "M2\nSurf.")
  pathway_order <- rev(selected_pathways)  # bottom to top

  library(cowplot)

  # Build data for BOTH directions in one combined figure
  all_df <- data.frame()
  for (direction in c("Outgoing", "Incoming")) {
    sig_list <- if (direction == "Outgoing") out_list2 else in_list2
    sig_base <- if (direction == "Outgoing") base_out else base_in

    # Baseline
    base_mat <- sig_list[[baseline_cond]]
    for (i in seq_along(selected_pathways)) {
      for (j in seq_along(broad_groups)) {
        all_df <- rbind(all_df, data.frame(
          direction = direction, pathway = selected_pathways[i],
          celltype = col_labels[j], value = base_mat[i, j],
          panel = cond_display[baseline_cond],
          type = "baseline", stringsAsFactors = FALSE))
      }
    }
    # Differences
    for (cond in diff_conds) {
      d <- sig_list[[cond]] - sig_base
      for (i in seq_along(selected_pathways)) {
        for (j in seq_along(broad_groups)) {
          all_df <- rbind(all_df, data.frame(
            direction = direction, pathway = selected_pathways[i],
            celltype = col_labels[j], value = d[i, j],
            panel = cond_display[cond],
            type = "diff", stringsAsFactors = FALSE))
        }
      }
    }
  }

  all_df$pathway <- factor(all_df$pathway, levels = pathway_order)
  all_df$celltype <- factor(all_df$celltype, levels = col_labels)
  panel_order <- c(cond_display[baseline_cond], cond_display[diff_conds])
  all_df$panel <- factor(all_df$panel, levels = panel_order)

  # Global max for consistent scales
  max_base_val <- max(all_df$value[all_df$type == "baseline"], na.rm = TRUE)
  max_diff_val <- max(abs(all_df$value[all_df$type == "diff"]), na.rm = TRUE)
  if (max_diff_val == 0) max_diff_val <- 0.5

  # Helper to build one row (outgoing or incoming)
  make_row <- function(dir_name, show_strip = TRUE) {
    df_base <- all_df[all_df$type == "baseline" & all_df$direction == dir_name, ]
    df_diff <- all_df[all_df$type == "diff" & all_df$direction == dir_name, ]

    p_base <- ggplot(df_base, aes(x = celltype, y = pathway, fill = value)) +
      geom_tile(color = "grey40", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 2.8, color = "white") +
      scale_fill_viridis_c(option = "viridis", limits = c(0, max_base_val * 1.1)) +
      facet_wrap(~panel, nrow = 1) +
      labs(x = "", y = "") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 9, face = "bold"),
            strip.text = if (show_strip) element_text(size = 8, face = "bold") else element_blank(),
            panel.grid = element_blank(), legend.position = "none",
            plot.margin = margin(2, 2, 2, 2))

    p_diff <- ggplot(df_diff, aes(x = celltype, y = pathway, fill = value)) +
      geom_tile(color = "grey40", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
      scale_fill_gradientn(colors = c("#2166AC", "#4393C3", "#92C5DE",
                                      "#D1E5F0", "white", "#FDDBC7",
                                      "#F4A582", "#D6604D", "#B2182B"),
                           limits = c(-max_diff_val, max_diff_val)) +
      facet_wrap(~panel, nrow = 1) +
      labs(x = "", y = "") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_blank(),
            strip.text = if (show_strip) element_text(size = 8, face = "bold") else element_blank(),
            panel.grid = element_blank(), legend.position = "none",
            plot.margin = margin(2, 2, 2, 2))

    plot_grid(p_base, p_diff, nrow = 1, rel_widths = c(1.3, 5), align = "h", axis = "tb")
  }

  row_out <- make_row("Outgoing", show_strip = TRUE)
  row_in  <- make_row("Incoming", show_strip = FALSE)

  # Shared legends
  leg_base <- cowplot::get_legend(
    ggplot(data.frame(x = 1, y = 1, v = 0), aes(x, y, fill = v)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Baseline\n(Raw)", option = "viridis",
                           limits = c(0, max_base_val * 1.1)) +
      theme(legend.position = "right", legend.key.height = unit(1.2, "cm"),
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8)))

  leg_diff <- cowplot::get_legend(
    ggplot(data.frame(x = 1, y = 1, v = 0), aes(x, y, fill = v)) +
      geom_tile() +
      scale_fill_gradientn(name = "Diff\n(vs Base)",
                           colors = c("#2166AC", "#4393C3", "#92C5DE",
                                      "#D1E5F0", "white", "#FDDBC7",
                                      "#F4A582", "#D6604D", "#B2182B"),
                           limits = c(-max_diff_val, max_diff_val)) +
      theme(legend.position = "right", legend.key.height = unit(1.2, "cm"),
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8)))

  legends <- plot_grid(leg_base, leg_diff, ncol = 1, rel_heights = c(1, 1))

  # Stack outgoing + incoming, add shared legends on right
  rows_stacked <- plot_grid(
    row_out + draw_label("Outgoing Signal:", x = 0.01, y = 0.95,
                         hjust = 0, vjust = 1, fontface = "bold", size = 11),
    row_in  + draw_label("Incoming Signal:", x = 0.01, y = 0.95,
                         hjust = 0, vjust = 1, fontface = "bold", size = 11),
    ncol = 1, rel_heights = c(1, 1))

  final <- plot_grid(rows_stacked, legends, nrow = 1, rel_widths = c(6, 0.7))

  pdf(file.path(fig_dir, "fig3B_difference_heatmaps.pdf"), width = 20, height = 7.5)
  print(final)
  dev.off()
  cat("  Saved fig3B_difference_heatmaps.pdf\n")
}

# ==============================================================================
# Fig 3D: CCL Ligand Senders & Receivers per cell type
# Left: red heatmap (senders), Right: blue heatmap (receivers)
# Averaged across all conditions
# ==============================================================================
cat("\nGenerating Fig 3D: CCL Ligand Senders & Receivers...\n")

tryCatch({
  # Cell type display order: C3, C5, C0, C1, C2, C4, M1, M2 Funct., M2 Surf.
  ct_order <- c("Eos_C3", "Eos_C5", "Eos_C0", "Eos_C1", "Eos_C2", "Eos_C4",
                "Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")
  ct_display <- c("Eos_C3" = "C3", "Eos_C5" = "C5", "Eos_C0" = "C0",
                  "Eos_C1" = "C1", "Eos_C2" = "C2", "Eos_C4" = "C4",
                  "Mac_M1" = "M1", "Mac_M2_Func" = "M2\nFunct.",
                  "Mac_M2_Surf" = "M2\nSurf.")

  # Collect ligand-level sender/receiver probabilities across conditions
  sender_all <- list()
  receiver_all <- list()

  for (cond in names(cellchat_objects)) {
    cc <- cellchat_objects[[cond]]
    if (is.null(cc)) next

    # Get CCL communications via standard CellChat API
    ccl_comm <- tryCatch(
      subsetCommunication(cc, signaling = "CCL"),
      error = function(e) NULL
    )
    if (is.null(ccl_comm) || nrow(ccl_comm) == 0) {
      cat(sprintf("    No CCL communications in %s\n", cond))
      next
    }

    for (idx in seq_len(nrow(ccl_comm))) {
      ligand      <- ccl_comm$ligand[idx]
      source_ct   <- ccl_comm$source[idx]
      target_ct   <- ccl_comm$target[idx]
      prob_val    <- ccl_comm$prob[idx]

      # Sender: the cell type sending the ligand
      if (source_ct %in% ct_order) {
        sender_all <- c(sender_all, list(data.frame(
          ligand = ligand, celltype = source_ct, value = prob_val,
          condition = cond, stringsAsFactors = FALSE)))
      }
      # Receiver: the cell type receiving the ligand
      if (target_ct %in% ct_order) {
        receiver_all <- c(receiver_all, list(data.frame(
          ligand = ligand, celltype = target_ct, value = prob_val,
          condition = cond, stringsAsFactors = FALSE)))
      }
    }
  }

  if (length(sender_all) == 0) stop("No CCL sender data found across any condition")

  sender_df <- do.call(rbind, sender_all)
  receiver_df <- do.call(rbind, receiver_all)

  # Average across conditions
  avg_sender <- sender_df %>%
    group_by(ligand, celltype) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  avg_receiver <- receiver_df %>%
    group_by(ligand, celltype) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  # Keep top ligands (by total signal)
  top_ligands_s <- avg_sender %>% group_by(ligand) %>%
    summarise(total = sum(value)) %>% arrange(desc(total)) %>%
    filter(total > 0) %>% head(7) %>% pull(ligand)
  top_ligands_r <- avg_receiver %>% group_by(ligand) %>%
    summarise(total = sum(value)) %>% arrange(desc(total)) %>%
    filter(total > 0) %>% head(7) %>% pull(ligand)

  avg_sender <- avg_sender %>% filter(ligand %in% top_ligands_s)
  avg_receiver <- avg_receiver %>% filter(ligand %in% top_ligands_r)

  # Factor ordering
  avg_sender$ligand <- factor(avg_sender$ligand, levels = rev(top_ligands_s))
  avg_sender$celltype <- factor(avg_sender$celltype, levels = ct_order)
  avg_receiver$ligand <- factor(avg_receiver$ligand, levels = rev(top_ligands_r))
  avg_receiver$celltype <- factor(avg_receiver$celltype, levels = ct_order)

  ct_labels_vec <- ct_display[ct_order]

  # Sender heatmap (red)
  p_sender <- ggplot(avg_sender, aes(x = celltype, y = ligand, fill = value)) +
    geom_tile(color = "grey70", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = "") +
    scale_x_discrete(labels = ct_labels_vec) +
    labs(title = "CCL Ligand Senders", x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 10, face = "italic"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          panel.grid = element_blank(),
          legend.key.height = unit(1.5, "cm"))

  # Receiver heatmap (blue)
  p_receiver <- ggplot(avg_receiver, aes(x = celltype, y = ligand, fill = value)) +
    geom_tile(color = "grey70", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#2166AC", name = "") +
    scale_x_discrete(labels = ct_labels_vec) +
    labs(title = "CCL Ligand Receivers", x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 10, face = "italic"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          panel.grid = element_blank(),
          legend.key.height = unit(1.5, "cm"))

  pdf(file.path(fig_dir, "fig3D_ccl_ligand_sender_receiver.pdf"), width = 16, height = 5)
  print(plot_grid(p_sender, p_receiver, nrow = 1, rel_widths = c(1, 1)))
  dev.off()
  cat("  Saved fig3D_ccl_ligand_sender_receiver.pdf\n")

}, error = function(e) cat(sprintf("  Fig 3D failed: %s\n", e$message)))

cat("\n10b CellChat figures complete.\n")
cat("Output:", fig_dir, "\n")
