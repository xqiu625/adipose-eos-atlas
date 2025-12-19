# scripts/R05_20250529/05_EosSub_CellChat_Chord/02_chord_diagram_by_EosSubclusters.R

library(CellChat)
library(circlize)
library(ggplot2)

# Set input/output paths
input_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R03_draft_v2/Figure2/CellChat"
output_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R05_20250529/05_EosSub_CellChat_Chord"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Eosinophil clusters to process
eos_clusters <- 0:6  # Cluster labels in file names are 0 to 6

for (eos_clust in eos_clusters) {
  cat("Generating chord diagram for Eosinophil Cluster", eos_clust, "...\n")
  
  # Load CellChat object
  rds_path <- file.path(input_dir, paste0("EosC", eos_clust, "_mac_cellchat.rds"))
  if (!file.exists(rds_path)) {
    warning("File not found:", rds_path)
    next
  }
  
  cellchat <- readRDS(rds_path)
  
  # Debug: Check the structure of the CellChat object
  cat("  Cell types in CellChat object:", paste(levels(cellchat@idents), collapse = ", "), "\n")
  cat("  Dimensions of weight matrix:", dim(cellchat@net$weight), "\n")
  
  # Get cell types from the CellChat object
  cell_types <- levels(cellchat@idents)
  
  # Force consistent labeling: replace any eosinophil cell type with standardized name
  # that matches the file naming convention (C0-C6)
  standardized_cell_types <- cell_types
  eos_idx <- grepl("Eos|eosinophil", cell_types, ignore.case = TRUE)
  if (any(eos_idx)) {
    standardized_cell_types[eos_idx] <- paste0("Eosinophil_C", eos_clust)
  }
  
  # Update the cellchat object with corrected labels
  levels(cellchat@idents) <- standardized_cell_types
  dimnames(cellchat@net$weight) <- list(standardized_cell_types, standardized_cell_types)
  
  # Create named vector for grouping
  # Method 1: Simple approach - each cell type is its own group
  group_vector <- setNames(standardized_cell_types, standardized_cell_types)
  
  # Alternative Method 2: Group by cell type categories
  # Uncomment this section if you want to group Eosinophil subclusters together
  # group_vector <- setNames(
  #   ifelse(grepl("Eosinophil_C", standardized_cell_types), "Eosinophils", "Macrophages"),
  #   standardized_cell_types
  # )
  
  cat("  Group vector:", paste(names(group_vector), "=", group_vector, collapse = "; "), "\n")
  
  # Create and save chord diagram
  tryCatch({
    pdf(file.path(output_dir, paste0("EosC", eos_clust, "_mac_chord.pdf")), 
        width = 10, height = 8)
    
    # Method 1: Using netVisual_chord_cell with proper group vector
    netVisual_chord_cell(
      net = cellchat@net$weight,
      group = group_vector,
      title.name = paste0("Eosinophil Cluster ", eos_clust, " Interactions with Macrophages"),
      small.gap = 1.5,
      big.gap = 10
    )
    
    dev.off()
    
    cat("  ✅ Saved chord diagram to:", file.path(output_dir, paste0("EosC", eos_clust, "_mac_chord.pdf")), "\n")
    
  }, error = function(e) {
    dev.off()  # Make sure to close PDF device even if error occurs
    cat("  ❌ Error creating chord diagram for cluster", eos_clust, ":", e$message, "\n")
    
    # Try alternative approach using base circlize
    tryCatch({
      cat("  Trying alternative approach with circlize...\n")
      
      pdf(file.path(output_dir, paste0("EosC", eos_clust, "_mac_chord_alt.pdf")), 
          width = 10, height = 8)
      
      # Get the weight matrix
      weight_matrix <- cellchat@net$weight
      
      # Create chord diagram using circlize directly
      circos.clear()
      circos.par(start.degree = 90)
      
      chordDiagram(
        weight_matrix,
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(weight_matrix))))),
        transparency = 0.3
      )
      
      # Add labels
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                   facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      }, bg.border = NA)
      
      title(paste0("Eosinophil Cluster ", eos_clust, " Interactions with Macrophages"))
      
      circos.clear()
      dev.off()
      
      cat("  ✅ Saved alternative chord diagram to:", 
          file.path(output_dir, paste0("EosC", eos_clust, "_mac_chord_alt.pdf")), "\n")
      
    }, error = function(e2) {
      dev.off()
      cat("  ❌ Alternative approach also failed:", e2$message, "\n")
    })
  })
}

cat("✅ All chord diagrams processing complete.\n")

# Create combined chord diagram with all subclusters
cat("Creating combined chord diagram with all eosinophil subclusters...\n")

tryCatch({
  # Load all CellChat objects and combine their interaction data
  all_interactions <- list()
  all_cell_types <- c()
  
  for (eos_clust in eos_clusters) {
    rds_path <- file.path(input_dir, paste0("EosC", eos_clust, "_mac_cellchat.rds"))
    if (file.exists(rds_path)) {
      cellchat <- readRDS(rds_path)
      
      # Get weight matrix and standardize eosinophil cell type names
      weight_matrix <- cellchat@net$weight
      cell_types <- levels(cellchat@idents)
      
      # Force consistent labeling: replace any eosinophil cell type with standardized name
      new_names <- cell_types
      eos_idx <- grepl("Eos|eosinophil", cell_types, ignore.case = TRUE)
      if (any(eos_idx)) {
        new_names[eos_idx] <- paste0("Eosinophil_C", eos_clust)
      }
      
      dimnames(weight_matrix) <- list(new_names, new_names)
      all_interactions[[paste0("EosC", eos_clust)]] <- weight_matrix
      all_cell_types <- c(all_cell_types, new_names)
    }
  }
  
  if (length(all_interactions) > 0) {
    # Get unique cell types
    unique_cell_types <- unique(all_cell_types)
    
    # Create combined interaction matrix
    combined_matrix <- matrix(0, 
                             nrow = length(unique_cell_types), 
                             ncol = length(unique_cell_types),
                             dimnames = list(unique_cell_types, unique_cell_types))
    
    # Sum interactions from all clusters
    for (interaction_matrix in all_interactions) {
      current_types <- rownames(interaction_matrix)
      for (i in current_types) {
        for (j in current_types) {
          if (i %in% unique_cell_types && j %in% unique_cell_types) {
            combined_matrix[i, j] <- combined_matrix[i, j] + interaction_matrix[i, j]
          }
        }
      }
    }
    
    # Remove zero-sum rows/columns
    non_zero_idx <- rowSums(combined_matrix) > 0 | colSums(combined_matrix) > 0
    if (sum(non_zero_idx) > 0) {
      combined_matrix <- combined_matrix[non_zero_idx, non_zero_idx, drop = FALSE]
      
      # Create group vector for combined plot
      cell_types_final <- rownames(combined_matrix)
      group_vector_combined <- setNames(
        ifelse(grepl("Eosinophil_C", cell_types_final), "Eosinophils", "Macrophages"),
        cell_types_final
      )
      
      # Method 1: CellChat style
      tryCatch({
        pdf(file.path(output_dir, "All_EosClusters_combined_chord.pdf"), width = 12, height = 10)
        
        # Use circlize directly for better control
        circos.clear()
        circos.par(start.degree = 90, gap.after = c(rep(2, length(cell_types_final)-1), 10))
        
        # Create color scheme
        eos_colors <- rainbow(sum(grepl("Eosinophil_C", cell_types_final)), alpha = 0.7)
        mac_colors <- c("darkgreen", "forestgreen", "green3")[1:sum(!grepl("Eosinophil_C", cell_types_final))]
        
        sector_colors <- rep("gray", length(cell_types_final))
        sector_colors[grepl("Eosinophil_C", cell_types_final)] <- eos_colors
        sector_colors[!grepl("Eosinophil_C", cell_types_final)] <- mac_colors
        names(sector_colors) <- cell_types_final
        
        chordDiagram(
          combined_matrix,
          grid.col = sector_colors,
          transparency = 0.3,
          annotationTrack = "grid",
          preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(combined_matrix)))))
        )
        
        # Add sector labels (show full Eosinophil_C# labels)
        circos.track(track.index = 1, panel.fun = function(x, y) {
          circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                     CELL_META$sector.index, 
                     facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
        }, bg.border = NA)
        
        title("Combined Eosinophil Subclusters (C0-C6) Interactions with Macrophages", cex.main = 1.2)
        
        # Add legend
        legend("topright", 
               legend = c("Eosinophil Subclusters", "Macrophages"),
               fill = c("red", "darkgreen"),
               cex = 0.8, bty = "n")
        
        circos.clear()
        dev.off()
        
        cat("✅ Combined chord diagram saved to: All_EosClusters_combined_chord.pdf\n")
        
      }, error = function(e) {
        dev.off()
        cat("❌ Error creating combined chord diagram:", e$message, "\n")
      })
      
      # Method 2: Simplified version with cluster grouping
      tryCatch({
        pdf(file.path(output_dir, "All_EosClusters_grouped_chord.pdf"), width = 10, height = 8)
        
        # Aggregate by cluster groups (extract C0-C6 from Eosinophil_C# labels)
        eos_clusters_found <- unique(gsub("Eosinophil_C(\\d+)", "C\\1", cell_types_final[grepl("Eosinophil_C", cell_types_final)]))
        mac_types <- cell_types_final[!grepl("Eosinophil_C", cell_types_final)]
        
        # Create simplified matrix grouped by cluster (display as C0-C6)
        group_names <- c(eos_clusters_found, mac_types)
        group_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names),
                              dimnames = list(group_names, group_names))
        
        # Aggregate interactions by group
        for (i in 1:nrow(combined_matrix)) {
          for (j in 1:ncol(combined_matrix)) {
            from_cell <- rownames(combined_matrix)[i]
            to_cell <- colnames(combined_matrix)[j]
            
            # Determine group names (extract cluster labels C0-C6)
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
        
        # Remove zero rows/columns
        non_zero_groups <- rowSums(group_matrix) > 0 | colSums(group_matrix) > 0
        if (sum(non_zero_groups) > 0) {
          group_matrix <- group_matrix[non_zero_groups, non_zero_groups, drop = FALSE]
          
          circos.clear()
          circos.par(start.degree = 90)
          
          # Color scheme for groups
          final_groups <- rownames(group_matrix)
          group_colors <- rep("gray", length(final_groups))
          group_colors[grepl("^C\\d+$", final_groups)] <- rainbow(sum(grepl("^C\\d+$", final_groups)), alpha = 0.8)
          group_colors[!grepl("^C\\d+$", final_groups)] <- c("darkgreen", "forestgreen")[1:sum(!grepl("^C\\d+$", final_groups))]
          names(group_colors) <- final_groups
          
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
                       facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
          }, bg.border = NA)
          
          title("Eosinophil Subclusters (C0-C6) Interactions with Macrophages")rophages")rophages")
          
          circos.clear()
        }
        
        dev.off()
        
        cat("✅ Grouped chord diagram saved to: All_EosClusters_grouped_chord.pdf\n")
        
      }, error = function(e) {
        dev.off()
        cat("❌ Error creating grouped chord diagram:", e$message, "\n")
      })
      
      # Save combined interaction matrix
      write.csv(combined_matrix, file.path(output_dir, "combined_interaction_matrix.csv"))
      cat("✅ Combined interaction matrix saved to: combined_interaction_matrix.csv\n")
    }
  }
  
}, error = function(e) {
  cat("❌ Error creating combined chord diagrams:", e$message, "\n")
})

# Optional: Create a summary plot showing all clusters
cat("Creating summary visualization...\n")

tryCatch({
  # Load all successful CellChat objects and create summary
  summary_data <- list()
  
  for (eos_clust in eos_clusters) {
    rds_path <- file.path(input_dir, paste0("EosC", eos_clust, "_mac_cellchat.rds"))
    if (file.exists(rds_path)) {
      cellchat <- readRDS(rds_path)
      
      # Extract total interaction strength
      total_strength <- sum(cellchat@net$weight, na.rm = TRUE)
      summary_data[[paste0("EosC", eos_clust)]] <- data.frame(
        Cluster = paste0("EosC", eos_clust),
        Total_Interactions = total_strength,
        Num_Cell_Types = length(levels(cellchat@idents))
      )
    }
  }
  
  if (length(summary_data) > 0) {
    summary_df <- do.call(rbind, summary_data)
    
    # Create summary plot
    p <- ggplot(summary_df, aes(x = Cluster, y = Total_Interactions)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Total Interaction Strength by Eosinophil Cluster",
           x = "Eosinophil Cluster",
           y = "Total Interaction Strength") +
      geom_text(aes(label = round(Total_Interactions, 2)), vjust = -0.5)
    
    ggsave(file.path(output_dir, "interaction_summary.pdf"), plot = p, width = 8, height = 6)
    
    # Save summary table
    write.csv(summary_df, file.path(output_dir, "interaction_summary.csv"), row.names = FALSE)
    
    cat("✅ Summary visualization saved.\n")
  }
  
}, error = function(e) {
  cat("Warning: Could not create summary visualization:", e$message, "\n")
})



