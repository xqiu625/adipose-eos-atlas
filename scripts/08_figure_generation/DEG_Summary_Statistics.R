# =====================================
# DEG SUMMARY STATISTICS FOR MANUSCRIPT
# Author: Xinru Qiu
# Date: December 2024
#
# Purpose: Calculate for each comparison:
#   - Number of significant DEGs (up and down)
#   - Average log2FC for up and down regulated genes
#
# Output: Summary table for manuscript text
# =====================================

library(dplyr)

# =====================================
# CONFIGURATION
# =====================================

# Input directory
deg_input_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R13_20250910/DEG_results"

# Output directory
output_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R21_20251211"

# DEG thresholds (same as used in original analysis)
logfc_threshold <- 0.25
pval_threshold <- 0.05

# =====================================
# MAIN ANALYSIS
# =====================================

# Find all DEG files
deg_files <- list.files(deg_input_dir, pattern = "_DEG_results\\.csv$", full.names = TRUE)

cat("========================================\n")
cat("DEG SUMMARY STATISTICS FOR MANUSCRIPT\n")
cat("========================================\n")
cat("Input directory:", deg_input_dir, "\n")
cat("Log2FC threshold:", logfc_threshold, "\n")
cat("P-value threshold:", pval_threshold, "\n")
cat("Found", length(deg_files), "DEG files\n\n")

# Initialize results dataframe
results <- data.frame(
  Cell_Type = character(),
  Comparison = character(),
  Group1 = character(),
  Group2 = character(),
  Total_Sig_DEGs = integer(),
  Upregulated_Count = integer(),
  Upregulated_Avg_Log2FC = numeric(),
  Downregulated_Count = integer(),
  Downregulated_Avg_Log2FC = numeric(),
  stringsAsFactors = FALSE
)

# Process each file
for (deg_file in deg_files) {

  filename <- basename(deg_file)

  # Parse filename: CellType_Comparison_DEG_results.csv
  # Examples:
  #   Macrophages_WT_F_vs_WT_M_DEG_results.csv
  #   Eosinophils_F_NotInfected_vs_F_Infected_DEG_results.csv

  filename_clean <- gsub("_DEG_results\\.csv$", "", filename)
  parts <- strsplit(filename_clean, "_")[[1]]

  cell_type <- parts[1]
  comparison_parts <- parts[2:length(parts)]

  # Find "vs" position to split groups
  vs_pos <- which(comparison_parts == "vs")
  if (length(vs_pos) == 0) {
    cat("Warning: Could not parse comparison in", filename, "\n")
    next
  }

  group1 <- paste(comparison_parts[1:(vs_pos-1)], collapse = "_")
  group2 <- paste(comparison_parts[(vs_pos+1):length(comparison_parts)], collapse = "_")
  comparison <- paste(group1, "vs", group2)

  # Read DEG data
  deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)

  # Filter significant DEGs
  sig_degs <- deg_data[abs(deg_data$avg_log2FC) >= logfc_threshold &
                        deg_data$p_val_adj < pval_threshold, ]

  # Separate up and down regulated
  upregulated <- sig_degs[sig_degs$avg_log2FC > 0, ]
  downregulated <- sig_degs[sig_degs$avg_log2FC < 0, ]

  # Calculate statistics
  total_sig <- nrow(sig_degs)
  up_count <- nrow(upregulated)
  down_count <- nrow(downregulated)

  up_avg_log2fc <- ifelse(up_count > 0, mean(upregulated$avg_log2FC), NA)
  down_avg_log2fc <- ifelse(down_count > 0, mean(downregulated$avg_log2FC), NA)

  # Add to results
  results <- rbind(results, data.frame(
    Cell_Type = cell_type,
    Comparison = comparison,
    Group1 = group1,
    Group2 = group2,
    Total_Sig_DEGs = total_sig,
    Upregulated_Count = up_count,
    Upregulated_Avg_Log2FC = round(up_avg_log2fc, 3),
    Downregulated_Count = down_count,
    Downregulated_Avg_Log2FC = round(down_avg_log2fc, 3),
    stringsAsFactors = FALSE
  ))
}

# Sort results
results <- results %>%
  arrange(Cell_Type, Comparison)

# =====================================
# PRINT RESULTS
# =====================================

cat("\n========================================\n")
cat("SUMMARY TABLE\n")
cat("========================================\n\n")

print(results, row.names = FALSE)

# =====================================
# GENERATE MANUSCRIPT TEXT
# =====================================

cat("\n\n========================================\n")
cat("MANUSCRIPT TEXT (Copy-paste ready)\n")
cat("========================================\n\n")

# Separate by cell type
for (ct in unique(results$Cell_Type)) {

  cat("---", ct, "---\n\n")

  ct_results <- results[results$Cell_Type == ct, ]

  for (i in 1:nrow(ct_results)) {
    row <- ct_results[i, ]

    # Generate text
    # Note: In DEG analysis, positive log2FC means higher in Group1 vs Group2
    # So "upregulated" = higher in Group1

    cat("Comparison:", row$Comparison, "\n")
    cat(sprintf("  Total significant DEGs: %d\n", row$Total_Sig_DEGs))
    cat(sprintf("  Higher in %s: %d genes (avg log2FC = %.2f)\n",
                row$Group1, row$Upregulated_Count, row$Upregulated_Avg_Log2FC))
    cat(sprintf("  Higher in %s: %d genes (avg log2FC = %.2f)\n",
                row$Group2, row$Downregulated_Count, abs(row$Downregulated_Avg_Log2FC)))
    cat("\n")

    # Alternative phrasing for manuscript
    cat("  Manuscript text option:\n")
    cat(sprintf("  \"%s showed %d DEGs enriched (avg log2FC = %.2f), while %s showed %d DEGs enriched (avg log2FC = %.2f).\"\n\n",
                row$Group1, row$Upregulated_Count, row$Upregulated_Avg_Log2FC,
                row$Group2, row$Downregulated_Count, abs(row$Downregulated_Avg_Log2FC)))
  }

  cat("\n")
}

# =====================================
# SAVE RESULTS
# =====================================

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save CSV
output_file <- file.path(output_dir, "DEG_Summary_Statistics_for_Manuscript.csv")
write.csv(results, output_file, row.names = FALSE)
cat("\nResults saved to:", output_file, "\n")

# Save formatted text file
text_output <- file.path(output_dir, "DEG_Summary_Manuscript_Text.txt")

sink(text_output)

cat("DEG SUMMARY STATISTICS FOR MANUSCRIPT\n")
cat("=====================================\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log2FC threshold:", logfc_threshold, "\n")
cat("P-value (adjusted) threshold:", pval_threshold, "\n\n")

for (ct in unique(results$Cell_Type)) {

  cat("\n===", toupper(ct), "===\n\n")

  ct_results <- results[results$Cell_Type == ct, ]

  for (i in 1:nrow(ct_results)) {
    row <- ct_results[i, ]

    cat(row$Comparison, ":\n")
    cat(sprintf("  - Total significant DEGs: %d\n", row$Total_Sig_DEGs))
    cat(sprintf("  - %s enriched: %d genes, avg log2FC = %.2f\n",
                row$Group1, row$Upregulated_Count, row$Upregulated_Avg_Log2FC))
    cat(sprintf("  - %s enriched: %d genes, avg log2FC = %.2f\n\n",
                row$Group2, row$Downregulated_Count, abs(row$Downregulated_Avg_Log2FC)))
  }
}

cat("\n\nREADY-TO-USE SENTENCES:\n")
cat("========================\n\n")

for (ct in unique(results$Cell_Type)) {

  cat("---", ct, "---\n\n")

  ct_results <- results[results$Cell_Type == ct, ]

  for (i in 1:nrow(ct_results)) {
    row <- ct_results[i, ]

    # Clean up group names for readability
    g1_clean <- gsub("_", " ", row$Group1)
    g2_clean <- gsub("_", " ", row$Group2)
    g1_clean <- gsub("NotInfected", "uninfected", g1_clean)
    g2_clean <- gsub("NotInfected", "uninfected", g2_clean)
    g1_clean <- gsub("Infected", "infected", g1_clean)
    g2_clean <- gsub("Infected", "infected", g2_clean)

    cat(sprintf("%s vs %s: %s had %d DEGs with average log2FC of %.2f, while %s had %d DEGs with average log2FC of %.2f.\n\n",
                g1_clean, g2_clean,
                g1_clean, row$Upregulated_Count, row$Upregulated_Avg_Log2FC,
                g2_clean, row$Downregulated_Count, abs(row$Downregulated_Avg_Log2FC)))
  }
}

sink()

cat("Text file saved to:", text_output, "\n")

cat("\n========================================\n")
cat("DONE!\n")
cat("========================================\n")
