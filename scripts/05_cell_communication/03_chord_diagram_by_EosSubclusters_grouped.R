# Combined Chord Diagram for All Eosinophil Clusters
library(CellChat)
library(circlize)
library(ggplot2)
library(RColorBrewer)

# Set input/output paths
input_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R03_draft_v2/Figure2/CellChat"
output_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/05_EosSub_CellChat_Chord"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Eosinophil clusters to process
eos_clusters <- 0:6

cat("Creating combined chord diagram with all eosinophil subclusters...\n")

# Load all CellChat objects and combine their interaction data
all_interactions <- list()
all_cell_types <- c()

for (eos_clust in eos_clusters) {
  rds_path <- file.path(input_dir, paste0("EosC", eos_clust, "_mac_cellchat.rds"))
  if (file.exists(rds_path)) {
    cat("Loading cluster", eos_clust, "...\n")
    cellchat <- readRDS(rds_path)
    
    # Get weight matrix and standardize cell type names
    weight_matrix <- cellchat@net$weight
    cell_types <- levels(cellchat@idents)
    
    cat("  Original cell types:", paste(cell_types, collapse = ", "), "\n")
    
    # Standardize eosinophil cell type names to match cluster number
    new_names <- cell_types
    eos_idx <- grepl("Eos|eosinophil", cell_types, ignore.case = TRUE)
    if (any(eos_idx)) {
      new_names[eos_idx] <- paste0("Eosinophil_C", eos_clust)
    }
    
    cat("  Standardized names:", paste(new_names, collapse = ", "), "\n")
    
    # Update matrix with new names
    dimnames(weight_matrix) <- list(new_names, new_names)
    all_interactions[[paste0("EosC", eos_clust)]] <- weight_matrix
    all_cell_types <- c(all_cell_types, new_names)
  } else {
    cat("File not found:", rds_path, "\n")
  }
}

if (length(all_interactions) == 0) {
  stop("No CellChat files found!")
}

# Get unique cell types
unique_cell_types <- unique(all_cell_types)
cat("\nUnique cell types found:", paste(unique_cell_types, collapse = ", "), "\n")

# Create combined interaction matrix
combined_matrix <- matrix(0, 
                         nrow = length(unique_cell_types), 
                         ncol = length(unique_cell_types),
                         dimnames = list(unique_cell_types, unique_cell_types))

# Sum interactions from all clusters
for (i in seq_along(all_interactions)) {
  interaction_matrix <- all_interactions[[i]]
  current_types <- rownames(interaction_matrix)
  
  cat("Processing interaction matrix", i, "with types:", paste(current_types, collapse = ", "), "\n")
  
  for (from_cell in current_types) {
    for (to_cell in current_types) {
      if (from_cell %in% unique_cell_types && to_cell %in% unique_cell_types) {
        combined_matrix[from_cell, to_cell] <- combined_matrix[from_cell, to_cell] + interaction_matrix[from_cell, to_cell]
      }
    }
  }
}

cat("\nCombined matrix dimensions:", dim(combined_matrix), "\n")
cat("Combined matrix row sums:", rowSums(combined_matrix), "\n")

# Remove zero-sum rows/columns
non_zero_idx <- rowSums(combined_matrix) > 0 | colSums(combined_matrix) > 0
if (sum(non_zero_idx) > 0) {
  combined_matrix <- combined_matrix[non_zero_idx, non_zero_idx, drop = FALSE]
  cat("After removing zero rows/cols, dimensions:", dim(combined_matrix), "\n")
  
  cell_types_final <- rownames(combined_matrix)
  cat("Final cell types:", paste(cell_types_final, collapse = ", "), "\n")
  
  # Method 1: Detailed view - all individual cell types
  cat("\nCreating detailed combined chord diagram...\n")
  
  pdf(file.path(output_dir, "All_EosClusters_combined_detailed.pdf"), width = 12, height = 10)
  
  # Create color scheme
  eos_types <- cell_types_final[grepl("Eosinophil_C", cell_types_final)]
  mac_types <- cell_types_final[!grepl("Eosinophil_C", cell_types_final)]
  
  # Colors for eosinophil clusters (different colors for each cluster)
  eos_colors <- rainbow(length(eos_types), alpha = 0.8)
  names(eos_colors) <- eos_types
  
  # Colors for macrophages (green palette)
  mac_colors <- brewer.pal(min(length(mac_types), 9), "Greens")[3:(2+length(mac_types))]
  names(mac_colors) <- mac_types
  
  # Combine colors
  sector_colors <- c(eos_colors, mac_colors)
  
  # Clear any previous plots
  circos.clear()
  
  # Set parameters
  circos.par(start.degree = 90, gap.after = c(rep(2, length(cell_types_final)-1), 10))
  
  # Create chord diagram
  chordDiagram(
    combined_matrix,
    grid.col = sector_colors,
    transparency = 0.3,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(combined_matrix)))))
  )
  
  # Add sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
               CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
  }, bg.border = NA)
  
  title("Combined Eosinophil Subclusters (C0-C6) Interactions with Macrophages", cex.main = 1.2)
  
  # Add legend
  legend("bottomright", 
         legend = c("Eosinophil Subclusters", "Macrophages"),
         fill = c("red", "darkgreen"),
         cex = 0.8, bty = "n")
  
  circos.clear()
  dev.off()
  
  cat("✅ Detailed combined chord diagram saved to: All_EosClusters_combined_detailed.pdf\n")
  
  # Method 2: Simplified grouped view
  cat("\nCreating simplified grouped chord diagram...\n")
  
  pdf(file.path(output_dir, "All_EosClusters_combined_grouped.pdf"), width = 10, height = 8)
  
  # Extract cluster numbers from eosinophil names
  eos_clusters_found <- unique(gsub("Eosinophil_C(\\d+)", "C\\1", eos_types))
  
  # Create simplified matrix grouped by cluster
  group_names <- c(eos_clusters_found, mac_types)
  group_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names),
                        dimnames = list(group_names, group_names))
  
  # Aggregate interactions by group
  for (i in 1:nrow(combined_matrix)) {
    for (j in 1:ncol(combined_matrix)) {
      from_cell <- rownames(combined_matrix)[i]
      to_cell <- colnames(combined_matrix)[j]
      
      # Determine group names
      if (grepl("Eosinophil_C", from_cell)) {
        from_group <- gsub("Eosinophil_C(\\d+)", "C\\1", from_cell)
      } else {
        from_group <- from_cell
      }
      
      if (grepl("Eosinophil_C", to_cell)) {
        to_group <- gsub("Eosinophil_C(\\d+)", "C\\1", to_cell)
      } else {
        to_group <- to_cell
      }
      
      if (from_group %in% group_names && to_group %in% group_names) {
        group_matrix[from_group, to_group] <- group_matrix[from_group, to_group] + combined_matrix[i, j]
      }
    }
  }
  
  # Remove zero rows/columns from grouped matrix
  non_zero_groups <- rowSums(group_matrix) > 0 | colSums(group_matrix) > 0
  if (sum(non_zero_groups) > 0) {
    group_matrix <- group_matrix[non_zero_groups, non_zero_groups, drop = FALSE]
    
    # Color scheme for groups
    final_groups <- rownames(group_matrix)
    eos_group_types <- final_groups[grepl("^C\\d+$", final_groups)]
    mac_group_types <- final_groups[!grepl("^C\\d+$", final_groups)]
    
    group_colors <- rep("gray", length(final_groups))
    if (length(eos_group_types) > 0) {
      group_colors[grepl("^C\\d+$", final_groups)] <- rainbow(length(eos_group_types), alpha = 0.8)
    }
    if (length(mac_group_types) > 0) {
      group_colors[!grepl("^C\\d+$", final_groups)] <- brewer.pal(min(length(mac_group_types), 9), "Greens")[3:(2+length(mac_group_types))]
    }
    names(group_colors) <- final_groups
    
    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(5, length(final_groups)-1), 15))
    
    chordDiagram(
      group_matrix,
      grid.col = group_colors,
      transparency = 0.3,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(group_matrix)))))
    )
    
    # Add labels
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                 facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
    }, bg.border = NA)
    
    title("Eosinophil Subclusters (C0-C6) Grouped Interactions with Macrophages", cex.main = 1.2)
    
    # Add legend
    legend("bottomright", 
           legend = c("Eosinophil Clusters", "Macrophages"),
           fill = c("red", "darkgreen"),
           cex = 0.8, bty = "n")
    
    circos.clear()
  }
  
  dev.off()
  
  cat("✅ Grouped chord diagram saved to: All_EosClusters_combined_grouped.pdf\n")
  
  # Save combined interaction matrix
  write.csv(combined_matrix, file.path(output_dir, "combined_interaction_matrix.csv"))
  cat("✅ Combined interaction matrix saved to: combined_interaction_matrix.csv\n")
  
  # Create summary statistics
  summary_stats <- data.frame(
    Cell_Type = rownames(combined_matrix),
    Total_Outgoing = rowSums(combined_matrix),
    Total_Incoming = colSums(combined_matrix),
    Total_Interactions = rowSums(combined_matrix) + colSums(combined_matrix),
    Is_Eosinophil = grepl("Eosinophil_C", rownames(combined_matrix))
  )
  
  write.csv(summary_stats, file.path(output_dir, "interaction_summary_stats.csv"), row.names = FALSE)
  
  # Create summary plot
  p <- ggplot(summary_stats, aes(x = Cell_Type, y = Total_Interactions, fill = Is_Eosinophil)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = "darkgreen", "TRUE" = "red"),
                     labels = c("Macrophages", "Eosinophils")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Total Interaction Strength by Cell Type",
         x = "Cell Type",
         y = "Total Interaction Strength",
         fill = "Cell Type Category") +
    geom_text(aes(label = round(Total_Interactions, 1)), vjust = -0.5, size = 3)
  
  ggsave(file.path(output_dir, "interaction_summary_plot.pdf"), plot = p, width = 10, height = 6)
  
  cat("✅ Summary statistics and plot saved.\n")
  
} else {
  cat("❌ No interactions found in combined matrix.\n")
}

cat("\n🎉 Combined chord diagram analysis complete!\n")
cat("📁 Output files saved to:", output_dir, "\n")
cat("📄 Files created:\n")
cat("  - All_EosClusters_combined_detailed.pdf\n")
cat("  - All_EosClusters_combined_grouped.pdf\n")
cat("  - combined_interaction_matrix.csv\n")
cat("  - interaction_summary_stats.csv\n")
cat("  - interaction_summary_plot.pdf\n")