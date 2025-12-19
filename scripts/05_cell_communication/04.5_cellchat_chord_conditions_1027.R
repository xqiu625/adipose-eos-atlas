library(CellChat)
library(circlize)

# Set working directory
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat")

# Define color scheme
group_colors <- c(
  "Eosinophils" = "#777B7E",
  "M1_Mac" = "#cb416b",
  "M2_Functional" = "#8FC50E",
  "M2_Surface" = "#5633AA"
)

# Function to create chord diagrams using circlize directly
create_chord_diagram_circlize <- function(cellchat_path, output_prefix) {
  
  cat("\n=== Processing:", output_prefix, "===\n")
  
  # Load CellChat object
  cellchat <- readRDS(cellchat_path)
  
  # Get interaction matrices
  mat_count <- cellchat@net$count
  mat_weight <- cellchat@net$weight
  
  cat("Groups:", paste(rownames(mat_count), collapse = ", "), "\n")
  cat("Total interactions (count):", sum(mat_count), "\n")
  cat("Total interactions (weight):", round(sum(mat_weight), 2), "\n")
  
  # Ensure colors match the matrix row/col names
  colors_ordered <- group_colors[rownames(mat_count)]
  
  # Create combined chord diagram
  pdf(paste0(output_prefix, "_chord_diagram.pdf"), width = 16, height = 8)
  par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))
  
  # COUNT diagram - using circlize directly
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(mat_count,
               grid.col = colors_ordered,
               transparency = 0.4,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = 1, font = 2)
  }, bg.border = NA)
  
  title(main = "Number of Interactions", cex.main = 1.3)
  circos.clear()
  
  # WEIGHT diagram
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(mat_weight,
               grid.col = colors_ordered,
               transparency = 0.4,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = 1, font = 2)
  }, bg.border = NA)
  
  title(main = "Interaction Strength", cex.main = 1.3)
  circos.clear()
  
  mtext(output_prefix, outer = TRUE, cex = 1.5, line = -1.5, font = 2)
  dev.off()
  
  # Create separate high-res count diagram
  pdf(paste0(output_prefix, "_chord_count.pdf"), width = 10, height = 10)
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(mat_count,
               grid.col = colors_ordered,
               transparency = 0.4,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.15))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = 1.2, font = 2)
  }, bg.border = NA)
  
  title(main = paste0(output_prefix, "\nNumber of Interactions"), 
        cex.main = 1.5, font.main = 2)
  circos.clear()
  dev.off()
  
  # Create separate high-res weight diagram
  pdf(paste0(output_prefix, "_chord_weight.pdf"), width = 10, height = 10)
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4)
  chordDiagram(mat_weight,
               grid.col = colors_ordered,
               transparency = 0.4,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.15))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, 
                adj = c(0, 0.5), cex = 1.2, font = 2)
  }, bg.border = NA)
  
  title(main = paste0(output_prefix, "\nInteraction Strength"), 
        cex.main = 1.5, font.main = 2)
  circos.clear()
  dev.off()
  
  cat("✓ Saved chord diagrams for", output_prefix, "\n")
  
  # Return summary statistics
  return(data.frame(
    sample = output_prefix,
    total_count = sum(mat_count),
    total_weight = round(sum(mat_weight), 2),
    mean_weight = round(mean(mat_weight[mat_weight > 0]), 3),
    max_weight = round(max(mat_weight), 3),
    n_interactions = sum(mat_count > 0),
    stringsAsFactors = FALSE
  ))
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

summary_stats <- list()
for (name in names(cellchat_files)) {
  file_path <- cellchat_files[[name]]
  if (file.exists(file_path)) {
    result <- create_chord_diagram_circlize(file_path, name)
    summary_stats[[name]] <- result
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Print summary table
cat("\n\n=== SUMMARY TABLE ===\n")
summary_df <- do.call(rbind, summary_stats)
rownames(summary_df) <- NULL
print(summary_df)

# Save summary table
write.csv(summary_df, "chord_diagram_summary.csv", row.names = FALSE)
cat("\nSummary saved to: chord_diagram_summary.csv\n")