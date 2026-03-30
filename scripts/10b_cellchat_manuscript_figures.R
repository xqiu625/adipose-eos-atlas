#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=cellchat_figs
#SBATCH --output=cellchat_figs_%j.out

# ==============================================================================
# 10b_cellchat_manuscript_figures.R
# Publication-quality CellChat figures for manuscript Figure 3 & Figure 4.
# Generates panels from atlas-level CellChat objects (eos + mac, 6 conditions).
#
# Input:  CellChat objects from 10_cellchat.R ($BASE/cellchat/cellchat_*.rds)
# Output: Publication PDFs in $BASE/cellchat/figures/
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(pheatmap)
  library(cowplot)
  library(RColorBrewer)
  library(viridis)
  library(grid)
  library(grDevices)
})

# Prevent X11 display issues in headless environments
options(bitmapType = 'cairo')
Sys.setenv(DISPLAY = "")
options(device = pdf)
grDevices::graphics.off()

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
fig_dir <- file.path(BASE, "cellchat", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Focus pathways for Fig 3B
selected_pathways <- c("APP", "CCL", "FN1", "CD45", "THBS")

# Mac types to include in figures (exclude Other_Mac for manuscript consistency)
mac_display <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")

# Condition display order and labels (6 conditions for atlas)
condition_order <- c("WT_Female_HFD", "WT_Male_HFD",
                     "KO_Female_HFD", "KO_Male_HFD",
                     "WT_Female_HFDInfected", "KO_Female_HFDInfected")

# Use cond_labels from 00_theme.R (has sex symbols)

# ==============================================================================
# 1. Load all CellChat objects
# ==============================================================================
cat("=== Loading CellChat objects ===\n")

cc_files <- list.files(cc_dir, pattern = "^cellchat_.*\\.rds$", full.names = TRUE)
cc_files <- cc_files[!grepl("merged", cc_files)]

cellchat_objects <- list()
for (f in cc_files) {
  cond <- gsub("cellchat_|\\.rds", "", basename(f))
  if (cond %in% condition_order) {
    cat(sprintf("  Loading %s ...\n", cond))
    cc <- readRDS(f)
    # Ensure centrality is computed
    tryCatch({
      if (is.null(cc@netP$centr)) {
        cc <- netAnalysis_computeCentrality(cc)
      }
    }, error = function(e) {
      cat(sprintf("    centrality computation failed for %s: %s\n", cond, e$message))
    })
    cellchat_objects[[cond]] <- cc
  }
}

cat(sprintf("Loaded %d conditions: %s\n", length(cellchat_objects),
            paste(names(cellchat_objects), collapse = ", ")))

# Identify available eos clusters dynamically
all_idents <- unique(unlist(lapply(cellchat_objects, function(cc) levels(cc@idents))))
eos_clusters <- sort(all_idents[grepl("^Eos_C", all_idents)])
mac_types_present <- sort(all_idents[grepl("^Mac_", all_idents)])

cat(sprintf("Eos clusters: %s\n", paste(eos_clusters, collapse = ", ")))
cat(sprintf("Mac types: %s\n", paste(mac_types_present, collapse = ", ")))

# Show which conditions have both eos + mac
for (cond in names(cellchat_objects)) {
  cc_tmp <- cellchat_objects[[cond]]
  idents_tmp <- levels(cc_tmp@idents)
  has_eos <- any(grepl("^Eos_C", idents_tmp))
  has_mac <- any(grepl("^Mac_", idents_tmp))
  cat(sprintf("  %s: %d idents (%s eos, %s mac)\n", cond, length(idents_tmp),
              ifelse(has_eos, "YES", "no"), ifelse(has_mac, "YES", "no")))
}

conds_with_both <- names(cellchat_objects)[sapply(cellchat_objects, function(cc) {
  ids <- levels(cc@idents)
  any(grepl("^Eos_C", ids)) && any(grepl("^Mac_", ids))
})]
cat(sprintf("Conditions with both Eos + Mac: %s\n", paste(conds_with_both, collapse = ", ")))

# ==============================================================================
# 2. Helper functions
# ==============================================================================

# --- Map internal cell type names to display names ---
display_name <- function(x) {
  mapping <- c(
    "Mac_M1"       = "M1",
    "Mac_M2_Func"  = "M2 Funct.",
    "Mac_M2_Surf"  = "M2 Surf.",
    "Mac_Other_Mac" = "Other Mac"
  )
  ifelse(x %in% names(mapping), mapping[x], x)
}

# --- Map to 4 chord groups ---
chord_group <- function(x) {
  dplyr::case_when(
    grepl("^Eos_", x)       ~ "Eosinophils",
    x == "Mac_M1"           ~ "M1_Mac",
    x == "Mac_M2_Func"      ~ "M2_Functional",
    x == "Mac_M2_Surf"      ~ "M2_Surface",
    TRUE                     ~ NA_character_
  )
}

# --- Extract pathway-level signaling strength ---
extract_signaling_strength <- function(cc) {
  pathway_names <- cc@netP$pathways
  cell_groups <- levels(cc@idents)
  prob_pathways <- cc@netP$prob

  outgoing <- matrix(0, nrow = length(pathway_names), ncol = length(cell_groups),
                     dimnames = list(pathway_names, cell_groups))
  incoming <- matrix(0, nrow = length(pathway_names), ncol = length(cell_groups),
                     dimnames = list(pathway_names, cell_groups))

  for (pw in pathway_names) {
    if (pw %in% dimnames(prob_pathways)[[3]]) {
      pw_prob <- prob_pathways[, , pw]
      outgoing[pw, ] <- rowSums(pw_prob)
      incoming[pw, ] <- colSums(pw_prob)
    }
  }

  available_pw <- intersect(selected_pathways, pathway_names)
  if (length(available_pw) == 0) return(NULL)

  list(
    outgoing = outgoing[available_pw, , drop = FALSE],
    incoming = incoming[available_pw, , drop = FALSE]
  )
}

# --- Aggregate cell types to 4 display groups for Fig 3B ---
aggregate_to_display <- function(mat, eos_clusters, mac_display) {
  display_cols <- c("Eos", "M1", "M2 Funct.", "M2 Surf.")
  agg <- matrix(0, nrow = nrow(mat), ncol = length(display_cols),
                dimnames = list(rownames(mat), display_cols))

  eos_present <- intersect(eos_clusters, colnames(mat))
  if (length(eos_present) > 0) {
    agg[, "Eos"] <- rowSums(mat[, eos_present, drop = FALSE])
  }
  if ("Mac_M1" %in% colnames(mat))       agg[, "M1"]         <- mat[, "Mac_M1"]
  if ("Mac_M2_Func" %in% colnames(mat))  agg[, "M2 Funct."]  <- mat[, "Mac_M2_Func"]
  if ("Mac_M2_Surf" %in% colnames(mat))  agg[, "M2 Surf."]   <- mat[, "Mac_M2_Surf"]

  agg
}

# --- Single heatmap panel for Fig 3B ---
create_single_heatmap <- function(data_matrix, title, color_limits,
                                   show_y_axis = TRUE, show_x_axis = TRUE,
                                   use_viridis = FALSE) {
  df <- melt(data_matrix)
  colnames(df) <- c("Pathway", "CellType", "Value")

  pathway_order <- rev(selected_pathways[selected_pathways %in% rownames(data_matrix)])
  df$Pathway <- factor(df$Pathway, levels = pathway_order)
  df$CellType <- factor(df$CellType, levels = colnames(data_matrix))

  if (use_viridis) {
    color_palette <- viridis_palette  # from 00_theme.R
    actual_limits <- c(0, max(color_limits))
    text_colors <- ifelse(df$Value > max(color_limits) * 0.6, "black", "white")
  } else {
    color_palette <- diverging_palette  # from 00_theme.R
    actual_limits <- color_limits
    text_colors <- ifelse(abs(df$Value) > max(abs(color_limits)) * 0.5, "white", "black")
  }

  p <- ggplot(df, aes(x = CellType, y = Pathway, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              size = 3, fontface = "bold", color = text_colors) +
    scale_fill_gradientn(colors = color_palette, limits = actual_limits,
                         na.value = "grey50") +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank())
  }
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank())
  }
  return(p)
}

# --- Extract CCL L-R pair probabilities from cc@net$prob ---
extract_ccl_lr_probs <- function(cc) {
  net_prob <- cc@net$prob
  if (is.null(net_prob)) return(NULL)

  lr_names <- dimnames(net_prob)[[3]]
  ccl_lr <- lr_names[grepl("CCL|Ccl", lr_names)]

  cat(sprintf("  net$prob L-R pairs matching CCL/Ccl: %d of %d total\n",
              length(ccl_lr), length(lr_names)))
  if (length(ccl_lr) > 0) cat("    Examples:", paste(head(ccl_lr, 5), collapse = ", "), "\n")

  if (length(ccl_lr) == 0) return(NULL)

  result <- list()
  for (lr in ccl_lr) {
    result[[lr]] <- net_prob[, , lr]
  }
  result
}


# ==============================================================================
# 3. Generate Figures
# ==============================================================================

# =============================================================================
# Fig 3A — CCL Chord Diagram (WT_Female_HFD baseline)
# =============================================================================
cat("\n=== Fig 3A: CCL Chord Diagram ===\n")

baseline_cond <- "WT_Female_HFD"
if (baseline_cond %in% names(cellchat_objects)) {
  cc <- cellchat_objects[[baseline_cond]]

  tryCatch({
    if ("CCL" %in% cc@netP$pathways) {
      ccl_prob <- cc@netP$prob[, , "CCL"]
      cell_types <- rownames(ccl_prob)

      # Map cell types to 4 groups
      group_map <- chord_group(cell_types)
      keep <- !is.na(group_map)
      ccl_prob_sub <- ccl_prob[keep, keep]
      groups <- group_map[keep]

      # Aggregate to 4x4 matrix
      unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
      agg_mat <- matrix(0, nrow = length(unique_groups), ncol = length(unique_groups),
                        dimnames = list(unique_groups, unique_groups))

      for (i in seq_along(groups)) {
        for (j in seq_along(groups)) {
          gi <- groups[i]
          gj <- groups[j]
          agg_mat[gi, gj] <- agg_mat[gi, gj] + ccl_prob_sub[i, j]
        }
      }

      # Draw chord diagram
      pdf(file.path(fig_dir, "fig3A_ccl_chord_diagram.pdf"), width = 10, height = 10)
      circos.clear()
      circos.par(start.degree = 90, gap.after = c(rep(2, nrow(agg_mat) - 1), 10))

      colors_ordered <- group_colors[rownames(agg_mat)]

      chordDiagram(
        agg_mat,
        grid.col = colors_ordered,
        transparency = 0.3,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        link.arr.type = "big.arrow",
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = 0.15),
        reduce = 0
      )

      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                    CELL_META$sector.index,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 1.0, font = 2)
      }, bg.border = NA)

      title(main = paste0("CCL Signaling Network \u2014 ", cond_labels[baseline_cond]),
            cex.main = 1.5, font.main = 2)
      circos.clear()
      dev.off()

      cat("  Saved: fig3A_ccl_chord_diagram.pdf\n")
    } else {
      cat("  CCL pathway not found in baseline condition\n")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  Fig 3A failed: %s\n", e$message))
  })
}


# =============================================================================
# Fig 3B — Outgoing/Incoming Signaling Difference Heatmaps
# =============================================================================
cat("\n=== Fig 3B: Difference Heatmaps ===\n")

tryCatch({
  sig_data <- list()
  for (cond in names(cellchat_objects)) {
    cc <- cellchat_objects[[cond]]
    strength <- extract_signaling_strength(cc)
    if (!is.null(strength)) {
      sig_data[[cond]] <- list(
        outgoing = aggregate_to_display(strength$outgoing, eos_clusters, mac_display),
        incoming = aggregate_to_display(strength$incoming, eos_clusters, mac_display)
      )
    }
  }

  if (length(sig_data) >= 2 && baseline_cond %in% names(sig_data)) {
    baseline <- sig_data[[baseline_cond]]
    comparison_conds <- setdiff(intersect(condition_order, names(sig_data)), baseline_cond)

    # Compute differences
    diff_out <- list()
    diff_in <- list()
    all_diff_values <- c()

    for (cond in comparison_conds) {
      common_pw <- intersect(rownames(sig_data[[cond]]$outgoing), rownames(baseline$outgoing))
      if (length(common_pw) > 0) {
        diff_out[[cond]] <- sig_data[[cond]]$outgoing[common_pw, , drop = FALSE] -
                            baseline$outgoing[common_pw, , drop = FALSE]
        diff_in[[cond]] <- sig_data[[cond]]$incoming[common_pw, , drop = FALSE] -
                           baseline$incoming[common_pw, , drop = FALSE]
        all_diff_values <- c(all_diff_values, c(diff_out[[cond]], diff_in[[cond]]))
      }
    }

    max_abs_diff <- max(abs(all_diff_values), na.rm = TRUE)
    global_limits <- c(-max_abs_diff, max_abs_diff)
    baseline_raw_max <- max(c(baseline$outgoing, baseline$incoming), na.rm = TRUE)

    cat(sprintf("  Difference limits: [%.3f, %.3f]\n", global_limits[1], global_limits[2]))

    # Build plot grid: 2 rows (outgoing, incoming) x (1 baseline + N comparisons)
    out_plots <- list()
    in_plots <- list()

    # Baseline — viridis raw values
    out_plots[[1]] <- create_single_heatmap(
      baseline$outgoing,
      paste0(cond_labels[baseline_cond], "\n(Baseline - Raw)"),
      c(0, baseline_raw_max),
      show_y_axis = TRUE, show_x_axis = FALSE, use_viridis = TRUE
    )
    in_plots[[1]] <- create_single_heatmap(
      baseline$incoming,
      paste0(cond_labels[baseline_cond], "\n(Baseline - Raw)"),
      c(0, baseline_raw_max),
      show_y_axis = TRUE, show_x_axis = TRUE, use_viridis = TRUE
    )

    # Comparisons — diverging differences
    for (i in seq_along(comparison_conds)) {
      cond <- comparison_conds[i]
      cond_label <- ifelse(cond %in% names(cond_labels), cond_labels[cond], cond)

      out_plots[[i + 1]] <- create_single_heatmap(
        diff_out[[cond]],
        paste0(cond_label, "\nvs Baseline"),
        global_limits,
        show_y_axis = FALSE, show_x_axis = FALSE, use_viridis = FALSE
      )
      in_plots[[i + 1]] <- create_single_heatmap(
        diff_in[[cond]],
        paste0(cond_label, "\nvs Baseline"),
        global_limits,
        show_y_axis = FALSE, show_x_axis = TRUE, use_viridis = FALSE
      )
    }

    n_cols <- length(out_plots)

    # Create legends
    legend_df_div <- data.frame(x = 1,
                                y = seq(global_limits[1], global_limits[2], length.out = 100),
                                fill = seq(global_limits[1], global_limits[2], length.out = 100))
    legend_plot_div <- ggplot(legend_df_div, aes(x = x, y = y, fill = fill)) +
      geom_tile() +
      scale_fill_gradientn(colors = diverging_palette, limits = global_limits,
                           name = "Diff\n(vs Base)") +
      theme_void() +
      theme(legend.position = "right",
            legend.key.height = unit(2, "cm"),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 8))

    legend_df_vir <- data.frame(x = 1,
                                y = seq(0, baseline_raw_max, length.out = 100),
                                fill = seq(0, baseline_raw_max, length.out = 100))
    legend_plot_vir <- ggplot(legend_df_vir, aes(x = x, y = y, fill = fill)) +
      geom_tile() +
      scale_fill_gradientn(colors = viridis_palette, limits = c(0, baseline_raw_max),
                           name = "Baseline\n(Raw)") +
      theme_void() +
      theme(legend.position = "right",
            legend.key.height = unit(2, "cm"),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 8))

    # Extract legends
    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp); legend_grob_div <- cowplot::get_legend(legend_plot_div); dev.off(); unlink(tmp)
    pdf(tmp); legend_grob_vir <- cowplot::get_legend(legend_plot_vir); dev.off(); unlink(tmp)

    combined_legends <- plot_grid(legend_grob_vir, legend_grob_div, ncol = 1,
                                  rel_heights = c(1, 1))

    # Assemble
    outgoing_row <- plot_grid(plotlist = out_plots, nrow = 1,
                              rel_widths = c(1.3, rep(1, n_cols - 1)))
    incoming_row <- plot_grid(plotlist = in_plots, nrow = 1,
                              rel_widths = c(1.3, rep(1, n_cols - 1)))

    # Row labels
    out_label <- ggdraw() + draw_label("Outgoing Signal", fontface = "bold", size = 12,
                                        angle = 90)
    in_label <- ggdraw() + draw_label("Incoming Signal", fontface = "bold", size = 12,
                                       angle = 90)

    out_with_label <- plot_grid(out_label, outgoing_row, nrow = 1, rel_widths = c(0.04, 1))
    in_with_label <- plot_grid(in_label, incoming_row, nrow = 1, rel_widths = c(0.04, 1))

    combined_rows <- plot_grid(out_with_label, in_with_label, ncol = 1,
                               rel_heights = c(1, 1.1))
    final_3B <- plot_grid(combined_rows, combined_legends, nrow = 1,
                           rel_widths = c(1, 0.12))

    ggsave(file.path(fig_dir, "fig3B_difference_heatmaps.pdf"),
           final_3B, width = 16, height = 5.9, units = "cm",
           dpi = 300, device = cairo_pdf, scale = 3)

    while (!is.null(dev.list())) dev.off()
    cat("  Saved: fig3B_difference_heatmaps.pdf\n")
  }
}, error = function(e) {
  while (!is.null(dev.list())) dev.off()
  cat(sprintf("  Fig 3B failed: %s\n", e$message))
})


# =============================================================================
# Fig 3C — CCL Signaling Role Heatmap (centrality)
# =============================================================================
cat("\n=== Fig 3C: CCL Signaling Role Heatmap ===\n")

if (baseline_cond %in% names(cellchat_objects)) {
  tryCatch({
    cc <- cellchat_objects[[baseline_cond]]
    centr <- cc@netP$centr

    if (!is.null(centr) && "CCL" %in% names(centr)) {
      ccl_centr <- centr[["CCL"]]

      role_names <- c("outdeg", "indeg", "flowbet", "info")
      role_labels <- c("Sender", "Receiver", "Mediator", "Influencer")

      cell_types <- names(ccl_centr[[role_names[1]]])
      centr_mat <- matrix(0, nrow = length(cell_types), ncol = length(role_names),
                          dimnames = list(cell_types, role_labels))
      for (k in seq_along(role_names)) {
        vals <- ccl_centr[[role_names[k]]]
        centr_mat[names(vals), role_labels[k]] <- vals
      }

      # --- Panel 1: Aggregated 4 groups ---
      groups <- chord_group(cell_types)
      keep <- !is.na(groups)
      unique_grps <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
      agg_centr <- matrix(0, nrow = length(unique_grps), ncol = ncol(centr_mat),
                          dimnames = list(unique_grps, role_labels))
      for (g in unique_grps) {
        idx <- which(groups == g)
        if (length(idx) > 0) {
          agg_centr[g, ] <- colMeans(centr_mat[idx, , drop = FALSE])
        }
      }

      # --- Panel 2: Replace Eos with Eos_C3 specifically ---
      panel2_types <- c("Eos_C3", intersect(mac_display, cell_types))
      panel2_types <- panel2_types[panel2_types %in% cell_types]
      panel2_mat <- centr_mat[panel2_types, , drop = FALSE]
      rownames(panel2_mat) <- display_name(rownames(panel2_mat))

      pdf(file.path(fig_dir, "fig3C_ccl_signaling_roles.pdf"), width = 8, height = 6)

      pheatmap(agg_centr,
               main = "CCL Signaling Roles (Aggregated Groups)",
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.3f",
               fontsize_number = 9, fontsize = 11,
               color = colorRampPalette(c("white", "#2ca02c", "#7b2d8e"))(100),
               border_color = "grey80")

      pheatmap(panel2_mat,
               main = "CCL Signaling Roles (Eos_C3 Focus)",
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.3f",
               fontsize_number = 9, fontsize = 11,
               color = colorRampPalette(c("white", "#2ca02c", "#7b2d8e"))(100),
               border_color = "grey80")

      dev.off()
      cat("  Saved: fig3C_ccl_signaling_roles.pdf\n")
    } else {
      cat("  CCL centrality data not available\n")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  Fig 3C failed: %s\n", e$message))
  })
}


# =============================================================================
# Fig 3D — CCL Ligand Sender/Receiver Heatmaps
# =============================================================================
cat("\n=== Fig 3D: CCL Ligand Sender/Receiver Heatmaps ===\n")

if (baseline_cond %in% names(cellchat_objects)) {
  tryCatch({
    cc <- cellchat_objects[[baseline_cond]]
    ccl_probs <- extract_ccl_lr_probs(cc)

    if (!is.null(ccl_probs)) {
      cell_types <- dimnames(cc@net$prob)[[1]]
      eos_present <- intersect(eos_clusters, cell_types)
      mac_present <- intersect(mac_display, cell_types)
      display_types <- c(eos_present, mac_present)

      ligands <- unique(sub("_.*", "", names(ccl_probs)))

      sender_mat <- matrix(0, nrow = length(ligands), ncol = length(display_types),
                           dimnames = list(ligands, display_name(display_types)))
      receiver_mat <- matrix(0, nrow = length(ligands), ncol = length(display_types),
                             dimnames = list(ligands, display_name(display_types)))

      for (lr in names(ccl_probs)) {
        lig <- sub("_.*", "", lr)
        prob_mat <- ccl_probs[[lr]]
        for (ct in display_types) {
          if (ct %in% rownames(prob_mat)) {
            sender_mat[lig, display_name(ct)] <- sender_mat[lig, display_name(ct)] +
                                                  sum(prob_mat[ct, ])
          }
          if (ct %in% colnames(prob_mat)) {
            receiver_mat[lig, display_name(ct)] <- receiver_mat[lig, display_name(ct)] +
                                                    sum(prob_mat[, ct])
          }
        }
      }

      # Remove all-zero rows
      sender_mat <- sender_mat[rowSums(sender_mat) > 0, , drop = FALSE]
      receiver_mat <- receiver_mat[rowSums(receiver_mat) > 0, , drop = FALSE]

      red_scale <- colorRampPalette(c("white", "red3"))(100)

      pdf(file.path(fig_dir, "fig3D_ccl_sender_receiver_heatmaps.pdf"),
          width = 18, height = 6)
      par(mfrow = c(1, 2))

      pheatmap(sender_mat,
               main = "CCL Ligand Senders",
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.3f",
               fontsize_number = 8, fontsize = 10,
               color = red_scale,
               border_color = "grey80")

      pheatmap(receiver_mat,
               main = "CCL Ligand Receivers",
               cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.3f",
               fontsize_number = 8, fontsize = 10,
               color = red_scale,
               border_color = "grey80")

      dev.off()
      cat("  Saved: fig3D_ccl_sender_receiver_heatmaps.pdf\n")
    } else {
      cat("  No CCL L-R probabilities found\n")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  Fig 3D failed: %s\n", e$message))
  })
}


# =============================================================================
# Fig 3E — Eos Subclusters to Mac Subtypes CCL Heatmap
# =============================================================================
cat("\n=== Fig 3E: Eos -> Mac CCL Heatmap ===\n")

if (baseline_cond %in% names(cellchat_objects)) {
  tryCatch({
    cc <- cellchat_objects[[baseline_cond]]

    cell_types <- levels(cc@idents)
    eos_present <- intersect(eos_clusters, cell_types)
    mac_present <- intersect(mac_display, cell_types)

    cat(sprintf("  Eos present: %s\n", paste(eos_present, collapse = ", ")))
    cat(sprintf("  Mac present: %s\n", paste(mac_present, collapse = ", ")))

    comm <- subsetCommunication(cc)
    ccl_comm <- comm[grepl("CCL", comm$pathway_name) &
                       comm$source %in% eos_present &
                       comm$target %in% mac_present, ]

    cat(sprintf("  CCL Eos->Mac interactions from subsetCommunication: %d\n", nrow(ccl_comm)))

    if (nrow(ccl_comm) > 0) {
      ccl_comm$pair_label <- paste0(gsub("Eos_", "", ccl_comm$source), "\u2192",
                                     display_name(ccl_comm$target))
      eos_to_mac_mat <- reshape2::acast(ccl_comm, interaction_name ~ pair_label,
                                         value.var = "prob", fun.aggregate = sum, fill = 0)

      eos_to_mac_mat <- eos_to_mac_mat[rowSums(eos_to_mac_mat) > 0, , drop = FALSE]
      eos_to_mac_mat <- eos_to_mac_mat[, colSums(eos_to_mac_mat) > 0, drop = FALSE]

      cat(sprintf("  After filtering: %d L-R pairs x %d cell pairs\n",
                  nrow(eos_to_mac_mat), ncol(eos_to_mac_mat)))

      if (nrow(eos_to_mac_mat) > 0 && ncol(eos_to_mac_mat) > 0) {
        green_scale <- colorRampPalette(c("white", "#238b45", "#00441b"))(100)

        pdf(file.path(fig_dir, "fig3E_eos_to_mac_ccl_heatmap.pdf"), width = 14, height = 6)
        pheatmap(eos_to_mac_mat,
                 main = "Eos Subclusters to Mac Subtypes \u2014 CCL",
                 cluster_rows = nrow(eos_to_mac_mat) >= 2,
                 cluster_cols = ncol(eos_to_mac_mat) >= 2,
                 display_numbers = TRUE, number_format = "%.3f",
                 fontsize_number = 7, fontsize = 10,
                 color = green_scale,
                 border_color = "grey80")
        dev.off()
        cat("  Saved: fig3E_eos_to_mac_ccl_heatmap.pdf\n")
      }
    } else {
      cat("  No CCL Eos->Mac interactions found\n")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  Fig 3E failed: %s\n", e$message))
  })
}


# =============================================================================
# Fig 3F — Mac Subtypes to Eos Subclusters CCL Heatmap
# =============================================================================
cat("\n=== Fig 3F: Mac -> Eos CCL Heatmap ===\n")

if (baseline_cond %in% names(cellchat_objects)) {
  tryCatch({
    cc <- cellchat_objects[[baseline_cond]]

    cell_types <- levels(cc@idents)
    eos_present <- intersect(eos_clusters, cell_types)
    mac_present <- intersect(mac_display, cell_types)

    comm <- subsetCommunication(cc)
    ccl_comm <- comm[grepl("CCL", comm$pathway_name) &
                       comm$source %in% mac_present &
                       comm$target %in% eos_present, ]

    cat(sprintf("  CCL Mac->Eos interactions from subsetCommunication: %d\n", nrow(ccl_comm)))

    if (nrow(ccl_comm) > 0) {
      ccl_comm$pair_label <- paste0(display_name(ccl_comm$source), "\u2192",
                                     gsub("Eos_", "", ccl_comm$target))
      mac_to_eos_mat <- reshape2::acast(ccl_comm, interaction_name ~ pair_label,
                                         value.var = "prob", fun.aggregate = sum, fill = 0)
      mac_to_eos_mat <- mac_to_eos_mat[rowSums(mac_to_eos_mat) > 0, , drop = FALSE]
      mac_to_eos_mat <- mac_to_eos_mat[, colSums(mac_to_eos_mat) > 0, drop = FALSE]

      cat(sprintf("  After filtering: %d L-R pairs x %d cell pairs\n",
                  nrow(mac_to_eos_mat), ncol(mac_to_eos_mat)))

      if (nrow(mac_to_eos_mat) > 0 && ncol(mac_to_eos_mat) > 0) {
        blue_red_scale <- colorRampPalette(c("white", "#4393C3", "#2166AC",
                                              "#D6604D", "#B2182B"))(100)

        pdf(file.path(fig_dir, "fig3F_mac_to_eos_ccl_heatmap.pdf"), width = 14, height = 6)
        pheatmap(mac_to_eos_mat,
                 main = "Mac Subtypes to Eos Subclusters \u2014 CCL",
                 cluster_rows = nrow(mac_to_eos_mat) >= 2,
                 cluster_cols = ncol(mac_to_eos_mat) >= 2,
                 display_numbers = TRUE, number_format = "%.3f",
                 fontsize_number = 7, fontsize = 10,
                 color = blue_red_scale,
                 border_color = "grey80")
        dev.off()
        cat("  Saved: fig3F_mac_to_eos_ccl_heatmap.pdf\n")
      }
    } else {
      cat("  No CCL Mac->Eos interactions found\n")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("  Fig 3F failed: %s\n", e$message))
  })
}


# =============================================================================
# Fig 4A — Protective Axis: Mac -> Eos via CCL6-CCR3 / CCL24-CCR3
# =============================================================================
cat("\n=== Fig 4A: Protective Axis (Mac -> Eos CCR3) ===\n")

protective_lr <- c("CCL6_CCR3", "CCL24_CCR3")

tryCatch({
  protective_rows <- list()
  for (cond in names(cellchat_objects)) {
    cc <- cellchat_objects[[cond]]
    comm <- subsetCommunication(cc)
    hits <- comm[comm$interaction_name %in% protective_lr &
                   grepl("^Mac_", comm$source) &
                   grepl("^Eos_C", comm$target) &
                   comm$source %in% mac_display, ]
    if (nrow(hits) > 0) {
      hits$condition <- cond
      protective_rows[[cond]] <- hits
    }
  }

  if (length(protective_rows) > 0) {
    prot_df <- do.call(rbind, protective_rows)
    prot_df$pair_label <- paste0(display_name(prot_df$source), " \u2192 ",
                                  prot_df$target, "\n", prot_df$interaction_name)
    prot_df$condition <- factor(prot_df$condition, levels = condition_order)
    prot_df$cond_label <- cond_labels[as.character(prot_df$condition)]
    prot_df$cond_label <- factor(prot_df$cond_label,
                                  levels = cond_labels[condition_order])

    pair_order <- prot_df %>%
      group_by(pair_label) %>%
      summarise(mean_prob = mean(prob), .groups = "drop") %>%
      arrange(mean_prob) %>%
      pull(pair_label)
    prot_df$pair_label <- factor(prot_df$pair_label, levels = pair_order)

    p4a <- ggplot(prot_df, aes(x = cond_label, y = pair_label,
                                size = prob, color = prob)) +
      geom_point() +
      scale_color_viridis_c(option = "viridis", name = "Probability") +
      scale_size_continuous(range = c(1, 8), name = "Probability") +
      labs(title = "Protective Axis: Macrophage \u2192 Eosinophil\n(CCL6-CCR3 / CCL24-CCR3)",
           x = "", y = "") +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 9),
        legend.position = "right"
      ) +
      guides(size = guide_legend(order = 1),
             color = guide_colorbar(order = 2))

    ggsave(file.path(fig_dir, "fig4A_protective_mac_to_eos_CCR3.pdf"),
           p4a, width = 10, height = 8, device = cairo_pdf)
    while (!is.null(dev.list())) dev.off()
    cat("  Saved: fig4A_protective_mac_to_eos_CCR3.pdf\n")
  } else {
    cat("  No protective CCR3 interactions found\n")
  }
}, error = function(e) {
  while (!is.null(dev.list())) dev.off()
  cat(sprintf("  Fig 4A failed: %s\n", e$message))
})


# =============================================================================
# Fig 4B — Pathogenic Axis: Eos -> Mac via CCL6-CCR2 (+ CCL6-CCR1, CCL9-CCR1)
# =============================================================================
cat("\n=== Fig 4B: Pathogenic Axis (Eos -> Mac CCR2) ===\n")

pathogenic_lr <- c("CCL6_CCR2", "CCL6_CCR1", "CCL9_CCR1")

tryCatch({
  pathogenic_rows <- list()
  for (cond in names(cellchat_objects)) {
    cc <- cellchat_objects[[cond]]
    comm <- subsetCommunication(cc)
    hits <- comm[comm$interaction_name %in% pathogenic_lr &
                   grepl("^Eos_C", comm$source) &
                   grepl("^Mac_", comm$target) &
                   comm$target %in% mac_display, ]
    if (nrow(hits) > 0) {
      hits$condition <- cond
      pathogenic_rows[[cond]] <- hits
    }
  }

  if (length(pathogenic_rows) > 0) {
    path_df <- do.call(rbind, pathogenic_rows)
    path_df$pair_label <- paste0(path_df$source, " \u2192 ",
                                  display_name(path_df$target), "\n",
                                  path_df$interaction_name)
    path_df$condition <- factor(path_df$condition, levels = condition_order)
    path_df$cond_label <- cond_labels[as.character(path_df$condition)]
    path_df$cond_label <- factor(path_df$cond_label,
                                  levels = cond_labels[condition_order])

    pair_order <- path_df %>%
      group_by(pair_label) %>%
      summarise(mean_prob = mean(prob), .groups = "drop") %>%
      arrange(mean_prob) %>%
      pull(pair_label)
    path_df$pair_label <- factor(path_df$pair_label, levels = pair_order)

    p4b <- ggplot(path_df, aes(x = cond_label, y = pair_label,
                                size = prob, color = prob)) +
      geom_point() +
      scale_color_viridis_c(option = "viridis", name = "Probability") +
      scale_size_continuous(range = c(1, 8), name = "Probability") +
      labs(title = "Pathogenic Axis: Eosinophil \u2192 Macrophage\n(CCL6-CCR2 / CCL6-CCR1 / CCL9-CCR1)",
           x = "", y = "") +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 9),
        legend.position = "right"
      ) +
      guides(size = guide_legend(order = 1),
             color = guide_colorbar(order = 2))

    ggsave(file.path(fig_dir, "fig4B_pathogenic_eos_to_mac_CCR2.pdf"),
           p4b, width = 10, height = 8, device = cairo_pdf)
    while (!is.null(dev.list())) dev.off()
    cat("  Saved: fig4B_pathogenic_eos_to_mac_CCR2.pdf\n")
  } else {
    cat("  No pathogenic CCR2 interactions found\n")
  }
}, error = function(e) {
  while (!is.null(dev.list())) dev.off()
  cat(sprintf("  Fig 4B failed: %s\n", e$message))
})


# =============================================================================
# Fig 4C — The Axis Switch: CCL6-CCR3 (protective) vs CCL6-CCR2 (pathogenic)
# =============================================================================
cat("\n=== Fig 4C: Axis Switch (CCR3 vs CCR2) ===\n")

tryCatch({
  avail_conds <- intersect(condition_order, names(cellchat_objects))

  ccr3_mat <- matrix(0, nrow = length(avail_conds), ncol = length(eos_clusters),
                     dimnames = list(avail_conds, eos_clusters))
  ccr2_mat <- matrix(0, nrow = length(avail_conds), ncol = length(eos_clusters),
                     dimnames = list(avail_conds, eos_clusters))

  for (cond in avail_conds) {
    cc <- cellchat_objects[[cond]]
    comm <- subsetCommunication(cc)
    cell_types <- levels(cc@idents)
    mac_here <- intersect(mac_display, cell_types)
    eos_here <- intersect(eos_clusters, cell_types)

    # Protective: Mac -> Eos via CCL6_CCR3
    ccr3_hits <- comm[comm$interaction_name == "CCL6_CCR3" &
                        comm$source %in% mac_here &
                        comm$target %in% eos_here, ]
    if (nrow(ccr3_hits) > 0) {
      for (e in eos_here) {
        sub <- ccr3_hits[ccr3_hits$target == e, ]
        if (nrow(sub) > 0) ccr3_mat[cond, e] <- sum(sub$prob)
      }
    }

    # Pathogenic: Eos -> Mac via CCL6_CCR2
    ccr2_hits <- comm[comm$interaction_name == "CCL6_CCR2" &
                        comm$source %in% eos_here &
                        comm$target %in% mac_here, ]
    if (nrow(ccr2_hits) > 0) {
      for (e in eos_here) {
        sub <- ccr2_hits[ccr2_hits$source == e, ]
        if (nrow(sub) > 0) ccr2_mat[cond, e] <- sum(sub$prob)
      }
    }
  }

  ccr3_df <- melt(ccr3_mat, varnames = c("Condition", "Cluster"), value.name = "Prob")
  ccr3_df$Axis <- "Protective: Mac \u2192 Eos\n(CCL6-CCR3)"
  ccr2_df <- melt(ccr2_mat, varnames = c("Condition", "Cluster"), value.name = "Prob")
  ccr2_df$Axis <- "Pathogenic: Eos \u2192 Mac\n(CCL6-CCR2)"

  switch_df <- rbind(ccr3_df, ccr2_df)
  switch_df$Condition <- factor(switch_df$Condition, levels = avail_conds)
  switch_df$cond_label <- cond_labels[as.character(switch_df$Condition)]
  switch_df$cond_label <- factor(switch_df$cond_label,
                                  levels = cond_labels[avail_conds])
  switch_df$Axis <- factor(switch_df$Axis,
                            levels = c("Protective: Mac \u2192 Eos\n(CCL6-CCR3)",
                                       "Pathogenic: Eos \u2192 Mac\n(CCL6-CCR2)"))

  max_val <- max(switch_df$Prob, na.rm = TRUE)
  switch_df$text_col <- ifelse(switch_df$Prob > max_val * 0.6, "white", "black")

  p4c <- ggplot(switch_df, aes(x = Cluster, y = cond_label, fill = Prob)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(Prob > 0, sprintf("%.3f", Prob), ""),
                  color = text_col),
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    facet_wrap(~ Axis, ncol = 2) +
    scale_fill_gradient2(
      low = "#08519c", mid = "#6baed6", high = "#d73027",
      midpoint = max_val * 0.35,
      limits = c(0, max_val),
      name = "Probability"
    ) +
    labs(title = "CCL6 Axis Switch by Eos Cluster",
         x = "Eosinophil Cluster", y = "") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.text.y = element_text(face = "bold", size = 10),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "right"
    )

  ggsave(file.path(fig_dir, "fig4C_axis_switch_CCR3_vs_CCR2.pdf"),
         p4c, width = 14, height = 6, device = cairo_pdf)
  while (!is.null(dev.list())) dev.off()
  cat("  Saved: fig4C_axis_switch_CCR3_vs_CCR2.pdf\n")
}, error = function(e) {
  while (!is.null(dev.list())) dev.off()
  cat(sprintf("  Fig 4C failed: %s\n", e$message))
})


# =============================================================================
# Fig 4D — Eos Cluster as CCL6 Hub (bar chart)
# =============================================================================
cat("\n=== Fig 4D: Eos Cluster Hub Analysis ===\n")

tryCatch({
  avail_conds <- intersect(condition_order, names(cellchat_objects))

  hub_rows <- list()
  for (cond in avail_conds) {
    cc <- cellchat_objects[[cond]]
    comm <- subsetCommunication(cc)
    cell_types <- levels(cc@idents)
    mac_here <- intersect(mac_display, cell_types)
    eos_here <- intersect(eos_clusters, cell_types)

    for (e in eos_here) {
      ccr3 <- comm[comm$interaction_name == "CCL6_CCR3" &
                      comm$source %in% mac_here & comm$target == e, ]
      ccr2 <- comm[comm$interaction_name == "CCL6_CCR2" &
                      comm$source == e & comm$target %in% mac_here, ]

      hub_rows[[length(hub_rows) + 1]] <- data.frame(
        condition = cond,
        cluster = e,
        CCR3_received = if (nrow(ccr3) > 0) sum(ccr3$prob) else 0,
        CCR2_sent     = if (nrow(ccr2) > 0) sum(ccr2$prob) else 0,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(hub_rows) > 0) {
    hub_df <- do.call(rbind, hub_rows)

    hub_avg <- hub_df %>%
      group_by(cluster) %>%
      summarise(
        CCR3_received = mean(CCR3_received),
        CCR2_sent     = mean(CCR2_sent),
        .groups = "drop"
      )

    hub_long <- melt(hub_avg, id.vars = "cluster",
                     variable.name = "Role", value.name = "Mean_Prob")
    hub_long$Role <- factor(hub_long$Role,
                            levels = c("CCR3_received", "CCR2_sent"),
                            labels = c("CCL6-CCR3 Received\n(Protective)",
                                       "CCL6-CCR2 Sent\n(Pathogenic)"))
    hub_long$cluster <- factor(hub_long$cluster, levels = sort(eos_clusters))

    cluster_fills <- c("CCL6-CCR3 Received\n(Protective)" = "#2166AC",
                       "CCL6-CCR2 Sent\n(Pathogenic)"     = "#B2182B")

    p4d <- ggplot(hub_long, aes(x = cluster, y = Mean_Prob, fill = Role)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      scale_fill_manual(values = cluster_fills, name = "") +
      labs(title = "Eos Clusters as CCL6 Signaling Hubs",
           subtitle = "Mean probability across all conditions",
           x = "Eosinophil Cluster", y = "Mean Interaction Probability") +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.position = "top",
        legend.text = element_text(size = 10, face = "bold")
      )

    ggsave(file.path(fig_dir, "fig4D_eos_cluster_hub_CCL6.pdf"),
           p4d, width = 10, height = 6, device = cairo_pdf)
    while (!is.null(dev.list())) dev.off()
    cat("  Saved: fig4D_eos_cluster_hub_CCL6.pdf\n")

    # --- Condition-resolved version ---
    hub_cond_long <- melt(hub_df, id.vars = c("condition", "cluster"),
                          variable.name = "Role", value.name = "Prob")
    hub_cond_long$Role <- factor(hub_cond_long$Role,
                                  levels = c("CCR3_received", "CCR2_sent"),
                                  labels = c("CCL6-CCR3 Received\n(Protective)",
                                             "CCL6-CCR2 Sent\n(Pathogenic)"))
    hub_cond_long$condition <- factor(hub_cond_long$condition, levels = condition_order)
    hub_cond_long$cond_label <- cond_labels[as.character(hub_cond_long$condition)]
    hub_cond_long$cond_label <- factor(hub_cond_long$cond_label,
                                        levels = cond_labels[condition_order])
    hub_cond_long$cluster <- factor(hub_cond_long$cluster, levels = sort(eos_clusters))

    p4d_cond <- ggplot(hub_cond_long, aes(x = cluster, y = Prob, fill = Role)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      facet_wrap(~ cond_label, nrow = 1) +
      scale_fill_manual(values = cluster_fills, name = "") +
      labs(title = "CCL6 Axis Switch by Eos Cluster and Condition",
           x = "Eosinophil Cluster", y = "Interaction Probability") +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 9),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "top",
        legend.text = element_text(size = 9, face = "bold")
      )

    ggsave(file.path(fig_dir, "fig4D_eos_cluster_hub_by_condition.pdf"),
           p4d_cond, width = 16, height = 5, device = cairo_pdf)
    while (!is.null(dev.list())) dev.off()
    cat("  Saved: fig4D_eos_cluster_hub_by_condition.pdf\n")
  } else {
    cat("  No CCL6 hub data found\n")
  }
}, error = function(e) {
  while (!is.null(dev.list())) dev.off()
  cat(sprintf("  Fig 4D failed: %s\n", e$message))
})


# ==============================================================================
# Verification: check output files
# ==============================================================================
cat("\n=== Verification ===\n")

expected_files <- c(
  "fig3A_ccl_chord_diagram.pdf",
  "fig3B_difference_heatmaps.pdf",
  "fig3C_ccl_signaling_roles.pdf",
  "fig3D_ccl_sender_receiver_heatmaps.pdf",
  "fig3E_eos_to_mac_ccl_heatmap.pdf",
  "fig3F_mac_to_eos_ccl_heatmap.pdf",
  "fig4A_protective_mac_to_eos_CCR3.pdf",
  "fig4B_pathogenic_eos_to_mac_CCR2.pdf",
  "fig4C_axis_switch_CCR3_vs_CCR2.pdf",
  "fig4D_eos_cluster_hub_CCL6.pdf",
  "fig4D_eos_cluster_hub_by_condition.pdf"
)

for (f in expected_files) {
  fpath <- file.path(fig_dir, f)
  if (file.exists(fpath)) {
    sz <- file.size(fpath)
    status <- if (sz > 5000) "OK" else "WARNING: small file"
    cat(sprintf("  [%s] %s (%s)\n", status, f, format(sz, big.mark = ",")))
  } else {
    cat(sprintf("  [MISSING] %s\n", f))
  }
}

cat(sprintf("\nDone! 10b_cellchat_manuscript_figures.R completed.\n"))
cat(sprintf("Output directory: %s\n", fig_dir))
sessionInfo()
