library(CellChat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(RColorBrewer)

# Prevent X11 display issues in headless environments
options(bitmapType = 'cairo')
Sys.setenv(DISPLAY = "")

# Disable interactive graphics
options(device = pdf)
grDevices::graphics.off()

# Set working directory
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat")

# Define color scheme for cell types
group_colors <- c(
  "Eosinophils" = "#777B7E",
  "M1_Mac" = "#cb416b",
  "M2_Functional" = "#8FC50E",
  "M2_Surface" = "#5633AA"
)

# Define specific pathways to include
selected_pathways <- c("APP", "CCL", "FN1", "CD45", "THBS")

# Define blue-white-red diverging color palette
diverging_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                        "white",
                                        "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

# Define viridis color palette for baseline (raw values)
viridis_colors <- colorRampPalette(c("#440154", "#482878", "#3e4a89", "#31688e",
                                      "#26838e", "#1f9e89", "#35b779", "#6ece58",
                                      "#b5de2b", "#fde725"))(100)

# Function to extract data from CellChat object
extract_signaling_data <- function(cellchat_path, output_prefix) {
  
  cat("\n=== Extracting data from:", output_prefix, "===\n")
  
  cellchat <- readRDS(cellchat_path)
  
  # Get pathway activity for each cell type
  pathway_names <- cellchat@netP$pathways
  cell_groups <- levels(cellchat@idents)
  
  cat("Found", length(pathway_names), "pathways\n")
  
  # Check if netP$prob exists
  if (is.null(cellchat@netP$prob)) {
    cat("ERROR: netP$prob is NULL\n")
    return(NULL)
  }
  
  # Initialize matrices
  outgoing_strength <- matrix(0, nrow = length(pathway_names), ncol = length(cell_groups))
  incoming_strength <- matrix(0, nrow = length(pathway_names), ncol = length(cell_groups))
  rownames(outgoing_strength) <- pathway_names
  rownames(incoming_strength) <- pathway_names
  colnames(outgoing_strength) <- cell_groups
  colnames(incoming_strength) <- cell_groups
  
  # Get pathway-level probabilities
  prob_pathways <- cellchat@netP$prob
  
  # Calculate outgoing and incoming for each pathway
  for (pathway in pathway_names) {
    if (pathway %in% dimnames(prob_pathways)[[3]]) {
      pathway_prob <- prob_pathways[,,pathway]
      outgoing_strength[pathway, ] <- rowSums(pathway_prob)
      incoming_strength[pathway, ] <- colSums(pathway_prob)
    }
  }
  
  # Filter for selected pathways only
  available_pathways <- intersect(selected_pathways, pathway_names)
  
  if (length(available_pathways) == 0) {
    cat("Warning: None of the selected pathways found\n")
    return(NULL)
  }
  
  outgoing_data <- outgoing_strength[available_pathways, , drop = FALSE]
  incoming_data <- incoming_strength[available_pathways, , drop = FALSE]
  
  return(list(
    outgoing = outgoing_data,
    incoming = incoming_data,
    name = output_prefix
  ))
}

# Function to create a single heatmap panel using ggplot2
create_single_heatmap <- function(data_matrix, title, color_limits, show_y_axis = TRUE, show_x_axis = TRUE,
                                   use_magma = FALSE) {

  # Convert matrix to long format
  df <- melt(data_matrix)
  colnames(df) <- c("Pathway", "CellType", "Value")

  # Ensure pathway order matches selected_pathways (reversed for display)
  pathway_order <- rev(selected_pathways[selected_pathways %in% rownames(data_matrix)])
  df$Pathway <- factor(df$Pathway, levels = pathway_order)

  # Choose color palette
  if (use_magma) {
    color_palette <- viridis_colors
    # For viridis, use limits from 0 to max for raw values
    actual_limits <- c(0, max(color_limits))
    # Text color based on value position in viridis scale
    text_colors <- ifelse(df$Value > max(color_limits) * 0.6, "black", "white")
  } else {
    color_palette <- diverging_colors
    actual_limits <- color_limits
    text_colors <- ifelse(abs(df$Value) > max(abs(color_limits)) * 0.5, "white", "black")
  }

  # Create heatmap
  p <- ggplot(df, aes(x = CellType, y = Pathway, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              size = 3, fontface = "bold",
              color = text_colors) +
    scale_fill_gradientn(colors = color_palette,
                         limits = actual_limits,
                         na.value = "grey50") +
    labs(title = title, x = "Cell Type", y = "Pathway") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  # Optionally hide y-axis
  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }

  # Optionally hide x-axis title
  if (!show_x_axis) {
    p <- p + theme(axis.title.x = element_blank())
  }

  return(p)
}

# Function to create difference comparison figure (relative to baseline)
create_difference_comparison_figure <- function(data_list, baseline_condition, comparison_conditions, 
                                                 output_filename, show_baseline_raw = FALSE) {
  
  baseline_data <- data_list[[baseline_condition]]
  
  # Prepare data for each condition
  all_conditions <- c(baseline_condition, comparison_conditions)
  n_conditions <- length(all_conditions)
  
  # Calculate differences and collect all values for global limits
  diff_outgoing_list <- list()
  diff_incoming_list <- list()
  all_diff_values <- c()
  
  for (cond in all_conditions) {
    data <- data_list[[cond]]
    
    if (cond == baseline_condition) {
      if (show_baseline_raw) {
        # Show raw values for baseline
        diff_outgoing_list[[cond]] <- data$outgoing
        diff_incoming_list[[cond]] <- data$incoming
      } else {
        # Show zeros for baseline (as reference point)
        diff_outgoing_list[[cond]] <- data$outgoing * 0
        diff_incoming_list[[cond]] <- data$incoming * 0
      }
    } else {
      # Calculate difference: comparison - baseline
      diff_outgoing_list[[cond]] <- data$outgoing - baseline_data$outgoing
      diff_incoming_list[[cond]] <- data$incoming - baseline_data$incoming
      all_diff_values <- c(all_diff_values, 
                           c(diff_outgoing_list[[cond]], diff_incoming_list[[cond]]))
    }
  }
  
  # Calculate symmetric limits based on differences
  if (length(all_diff_values) > 0) {
    max_abs_diff <- max(abs(all_diff_values), na.rm = TRUE)
    global_limits <- c(-max_abs_diff, max_abs_diff)
  } else {
    global_limits <- c(-1, 1)
  }
  
  cat("Difference limits:", global_limits[1], "to", global_limits[2], "\n")
  
  # Create lists to store plots
  outgoing_plots <- list()
  incoming_plots <- list()
  
  # Calculate baseline raw limits for magma color scale
  baseline_raw_limits <- c(0, max(c(baseline_data$outgoing, baseline_data$incoming), na.rm = TRUE))

  for (i in seq_along(all_conditions)) {
    cond <- all_conditions[i]
    show_y <- (i == 1)  # Only show y-axis for first column
    is_baseline <- (cond == baseline_condition)

    # Create title
    if (is_baseline) {
      if (show_baseline_raw) {
        out_title <- paste0(cond, "\n(Baseline - Raw)")
        in_title <- paste0(cond, "\n(Baseline - Raw)")
      } else {
        out_title <- paste0(cond, "\n(Baseline)")
        in_title <- paste0(cond, "\n(Baseline)")
      }
    } else {
      out_title <- paste0(cond, "\nvs ", baseline_condition, " (Outgoing)")
      in_title <- paste0(cond, "\nvs ", baseline_condition, " (Incoming)")
    }

    # Use magma for baseline raw values, diverging for differences
    use_magma_color <- is_baseline && show_baseline_raw
    color_lim <- if (use_magma_color) baseline_raw_limits else global_limits

    # Outgoing heatmap
    outgoing_plots[[i]] <- create_single_heatmap(
      diff_outgoing_list[[cond]],
      out_title,
      color_lim,
      show_y_axis = show_y,
      show_x_axis = FALSE,
      use_magma = use_magma_color
    )

    # Incoming heatmap
    incoming_plots[[i]] <- create_single_heatmap(
      diff_incoming_list[[cond]],
      in_title,
      color_lim,
      show_y_axis = show_y,
      show_x_axis = TRUE,
      use_magma = use_magma_color
    )
  }
  
  # Create diverging color legend (for differences)
  legend_df <- data.frame(x = 1, y = seq(global_limits[1], global_limits[2], length.out = 100),
                          fill = seq(global_limits[1], global_limits[2], length.out = 100))

  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = fill)) +
    geom_tile() +
    scale_fill_gradientn(colors = diverging_colors, limits = global_limits,
                         name = "Difference\n(vs Baseline)") +
    theme_void() +
    theme(legend.position = "right",
          legend.key.height = unit(2, "cm"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9)) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 8))

  # Extract just the legend (use temp pdf device to avoid X11)
  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  legend_grob <- suppressWarnings(cowplot::get_legend(legend_plot))
  dev.off()
  unlink(temp_pdf)

  # Create viridis legend for baseline raw values (if applicable)
  if (show_baseline_raw) {
    viridis_legend_df <- data.frame(x = 1, y = seq(0, baseline_raw_limits[2], length.out = 100),
                                     fill = seq(0, baseline_raw_limits[2], length.out = 100))

    viridis_legend_plot <- ggplot(viridis_legend_df, aes(x = x, y = y, fill = fill)) +
      geom_tile() +
      scale_fill_gradientn(colors = viridis_colors, limits = c(0, baseline_raw_limits[2]),
                           name = "Baseline\n(Raw Value)") +
      theme_void() +
      theme(legend.position = "right",
            legend.key.height = unit(2, "cm"),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 8))

    temp_pdf2 <- tempfile(fileext = ".pdf")
    pdf(temp_pdf2)
    viridis_legend_grob <- suppressWarnings(cowplot::get_legend(viridis_legend_plot))
    dev.off()
    unlink(temp_pdf2)

    # Combine both legends vertically
    combined_legends <- plot_grid(viridis_legend_grob, legend_grob, ncol = 1, rel_heights = c(1, 1))
  } else {
    combined_legends <- legend_grob
  }

  # Combine outgoing plots in one row
  outgoing_row <- plot_grid(plotlist = outgoing_plots, nrow = 1,
                            rel_widths = c(1.3, rep(1, n_conditions - 1)))

  # Combine incoming plots in one row
  incoming_row <- plot_grid(plotlist = incoming_plots, nrow = 1,
                            rel_widths = c(1.3, rep(1, n_conditions - 1)))

  # Combine rows
  combined_plots <- plot_grid(outgoing_row, incoming_row, ncol = 1, rel_heights = c(1, 1.1))

  # Add legend(s)
  final_plot <- plot_grid(combined_plots, combined_legends, nrow = 1, rel_widths = c(1, 0.12))
  
  # Save figure using cairo_pdf to avoid X11 issues
  # Dimensions: 16cm width x 5.9cm height
  ggsave(output_filename, final_plot, width = 16, height = 5.9, units = "cm", dpi = 300, device = cairo_pdf)
  
  # Ensure all graphics devices are closed to prevent X11 fallback
  while (!is.null(dev.list())) dev.off()
  
  cat("✓ Saved:", output_filename, "\n")
  
  invisible(final_plot)
}

# Process all CellChat objects
cellchat_files <- list(
  "WT_Female_HFD" = "WT_Female_HFD_cellchat.rds",
  "WT_Male_HFD" = "WT_Male_HFD_cellchat.rds",
  "WT_Female_HFDInfected" = "WT_Female_HFDInfected_cellchat.rds",
  "KO_Female_HFD" = "KO_Female_HFD_cellchat.rds",
  "KO_Male_HFD" = "KO_Male_HFD_cellchat.rds",
  "KO_Female_HFDInfected" = "KO_Female_HFDInfected_cellchat.rds",
  "eos_mac" = "eos_mac_cellchat.rds"
)

# FIRST PASS: Extract all data
cat("\n=== FIRST PASS: Extracting all data ===\n")
all_data <- list()

for (name in names(cellchat_files)) {
  file_path <- cellchat_files[[name]]
  if (file.exists(file_path)) {
    data <- extract_signaling_data(file_path, name)
    if (!is.null(data)) {
      all_data[[name]] <- data
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# SECOND PASS: Create difference comparison figures

cat("\n=== Creating difference heatmaps (relative to baselines) ===\n")

# Comparison 1: KO_Female_HFD and WT_Male_HFD vs WT_Female_HFD (as in screenshot)
if (all(c("WT_Female_HFD", "KO_Female_HFD", "WT_Male_HFD") %in% names(all_data))) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Female_HFD",
    comparison_conditions = c("KO_Female_HFD", "WT_Male_HFD"),
    output_filename = "Difference_vs_WT_Female_HFD.pdf",
    show_baseline_raw = TRUE  # Show raw values for baseline
  )
}

# Comparison 2: Same but show baseline as zeros
if (all(c("WT_Female_HFD", "KO_Female_HFD", "WT_Male_HFD") %in% names(all_data))) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Female_HFD",
    comparison_conditions = c("KO_Female_HFD", "WT_Male_HFD"),
    output_filename = "Difference_vs_WT_Female_HFD_baseline_zero.pdf",
    show_baseline_raw = FALSE  # Show zeros for baseline
  )
}

# Comparison 3: All KO vs their WT counterparts (Female HFD)
if (all(c("WT_Female_HFD", "KO_Female_HFD") %in% names(all_data))) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Female_HFD",
    comparison_conditions = c("KO_Female_HFD"),
    output_filename = "Difference_KO_vs_WT_Female_HFD.pdf",
    show_baseline_raw = TRUE
  )
}

# Comparison 4: Male HFD - KO vs WT
if (all(c("WT_Male_HFD", "KO_Male_HFD") %in% names(all_data))) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Male_HFD",
    comparison_conditions = c("KO_Male_HFD"),
    output_filename = "Difference_KO_vs_WT_Male_HFD.pdf",
    show_baseline_raw = TRUE
  )
}

# Comparison 5: Female HFDInfected - KO vs WT
if (all(c("WT_Female_HFDInfected", "KO_Female_HFDInfected") %in% names(all_data))) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Female_HFDInfected",
    comparison_conditions = c("KO_Female_HFDInfected"),
    output_filename = "Difference_KO_vs_WT_Female_HFDInfected.pdf",
    show_baseline_raw = TRUE
  )
}

# Comparison 6: All conditions vs WT_Female_HFD (specific order)
# Order: WT Female HFD (Baseline), WT Male HFD, KO Male HFD, KO Female HFD, KO Female Infected, WT Female Infected
ordered_comparisons <- c("WT_Male_HFD", "KO_Male_HFD", "KO_Female_HFD", "KO_Female_HFDInfected", "WT_Female_HFDInfected")
available_ordered <- ordered_comparisons[ordered_comparisons %in% names(all_data)]
if ("WT_Female_HFD" %in% names(all_data) && length(available_ordered) > 0) {
  create_difference_comparison_figure(
    all_data,
    baseline_condition = "WT_Female_HFD",
    comparison_conditions = available_ordered,
    output_filename = "Difference_All_vs_WT_Female_HFD.pdf",
    show_baseline_raw = TRUE
  )
}

# Print summary
cat("\n\n=== SUMMARY ===\n")
cat("Processed", length(all_data), "conditions\n")
cat("Selected pathways:", paste(selected_pathways, collapse = ", "), "\n")
cat("Color scheme: Blue (decreased vs baseline) - White (no change) - Red (increased vs baseline)\n")

cat("\n=== All difference heatmaps completed ===\n")
cat("Generated files:\n")
cat("  - Difference_vs_WT_Female_HFD.pdf (KO_Female_HFD & WT_Male_HFD vs baseline)\n")
cat("  - Difference_vs_WT_Female_HFD_baseline_zero.pdf (same, baseline shown as zeros)\n")
cat("  - Difference_KO_vs_WT_*.pdf (KO vs WT comparisons)\n")
cat("  - Difference_All_vs_WT_Female_HFD.pdf (all conditions vs WT_Female_HFD)\n")