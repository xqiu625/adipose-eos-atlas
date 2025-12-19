# =====================================
# GO PATHWAY ANALYSIS - ALL TERMS (NOT TOP 10)
# For Macrophage Comparisons
# Author: Xinru Qiu
# Date: December 2024
#
# Comparisons:
#   1. WT Female vs WT Male (Macrophages)
#   2. WT Female vs KO Female (Macrophages)
#   3. WT Male vs KO Male (Macrophages)
#   4. Uninfected vs Infected (Macrophages)
#
# Output: ALL significant GO terms with associated genes
# =====================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scales)
library(tidyr)

# =====================================
# CONFIGURATION
# =====================================

# Input/Output directories - UPDATE THESE PATHS FOR CLUSTER
deg_input_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R13_20250910/DEG_results"
output_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R21_20251211/go_all_terms_macrophages"

# DEG thresholds
logfc_threshold <- 0.25
pval_threshold <- 0.05

# GO analysis parameters - relaxed to get ALL terms
go_pvalue_cutoff <- 0.05
go_qvalue_cutoff <- 0.2

# Target comparisons for macrophages
# Based on actual file names in DEG_results:
#   Macrophages_WT_F_vs_WT_M_DEG_results.csv
#   Macrophages_WT_F_vs_KO_F_DEG_results.csv
#   Macrophages_WT_M_vs_KO_M_DEG_results.csv
#   Macrophages_F_NotInfected_vs_F_Infected_DEG_results.csv
target_comparisons <- c(
  "WT_F_vs_WT_M",              # WT Female vs WT Male
  "WT_F_vs_KO_F",              # WT Female vs KO Female
  "WT_M_vs_KO_M",              # WT Male vs KO Male
  "F_NotInfected_vs_F_Infected" # Uninfected vs Infected (Female)
)

# =====================================
# HELPER FUNCTIONS
# =====================================

# Convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  return(entrez_ids[!is.na(entrez_ids)])
}

# =====================================
# MAIN ANALYSIS FUNCTION
# =====================================

process_deg_file_all_terms <- function(deg_file_path, output_dir) {

  # Extract comparison info from filename
  filename <- basename(deg_file_path)
  filename_parts <- strsplit(gsub(".csv$", "", filename), "_")[[1]]

  cell_type <- filename_parts[1]
  comparison_parts <- filename_parts[2:(length(filename_parts)-2)]
  comparison_name <- paste(comparison_parts, collapse = "_")

  # Create readable comparison name
  readable_comparison <- gsub("_", " ", comparison_name)
  readable_comparison <- gsub("NotInfected", "Uninfected", readable_comparison)

  print(paste(rep("=", 60), collapse = ""))
  print(paste("Processing:", cell_type, "-", readable_comparison))
  print(paste(rep("=", 60), collapse = ""))

  # Read DEG results
  deg_data <- read.csv(deg_file_path, stringsAsFactors = FALSE)

  # Filter significant DEGs
  significant_degs <- deg_data[
    abs(deg_data$avg_log2FC) >= logfc_threshold &
    deg_data$p_val_adj < pval_threshold,
  ]

  if (nrow(significant_degs) < 5) {
    print(paste("  Warning: Only", nrow(significant_degs), "significant DEGs. Skipping."))
    return(NULL)
  }

  print(paste("  Total significant DEGs:", nrow(significant_degs)))

  # Separate upregulated and downregulated genes
  upregulated_genes <- significant_degs[significant_degs$avg_log2FC > logfc_threshold, "gene"]
  downregulated_genes <- significant_degs[significant_degs$avg_log2FC < -logfc_threshold, "gene"]

  print(paste("  Upregulated genes:", length(upregulated_genes)))
  print(paste("  Downregulated genes:", length(downregulated_genes)))

  # Convert to Entrez IDs
  up_entrez <- convert_to_entrez(upregulated_genes)
  down_entrez <- convert_to_entrez(downregulated_genes)

  print(paste("  Upregulated genes with Entrez ID:", length(up_entrez)))
  print(paste("  Downregulated genes with Entrez ID:", length(down_entrez)))

  results <- list()

  # Process upregulated genes for all GO categories
  if (length(up_entrez) >= 5) {
    for (go_ont in c("BP", "MF", "CC")) {
      tryCatch({
        go_result <- enrichGO(gene = up_entrez,
                             OrgDb = org.Mm.eg.db,
                             ont = go_ont,
                             pAdjustMethod = "BH",
                             pvalueCutoff = go_pvalue_cutoff,
                             qvalueCutoff = go_qvalue_cutoff,
                             readable = TRUE)

        if (!is.null(go_result) && nrow(go_result@result) > 0) {
          results[[paste0("Up_", go_ont)]] <- go_result
          print(paste("  Upregulated", go_ont, "pathways:", nrow(go_result@result)))
        }
      }, error = function(e) {
        print(paste("  Error in upregulated", go_ont, "analysis:", e$message))
      })
    }
  }

  # Process downregulated genes for all GO categories
  if (length(down_entrez) >= 5) {
    for (go_ont in c("BP", "MF", "CC")) {
      tryCatch({
        go_result <- enrichGO(gene = down_entrez,
                             OrgDb = org.Mm.eg.db,
                             ont = go_ont,
                             pAdjustMethod = "BH",
                             pvalueCutoff = go_pvalue_cutoff,
                             qvalueCutoff = go_qvalue_cutoff,
                             readable = TRUE)

        if (!is.null(go_result) && nrow(go_result@result) > 0) {
          results[[paste0("Down_", go_ont)]] <- go_result
          print(paste("  Downregulated", go_ont, "pathways:", nrow(go_result@result)))
        }
      }, error = function(e) {
        print(paste("  Error in downregulated", go_ont, "analysis:", e$message))
      })
    }
  }

  # Save ALL terms (not just top 10)
  if (length(results) > 0) {
    save_all_go_terms(results, readable_comparison, cell_type, comparison_name, output_dir)
  }

  return(results)
}

# =====================================
# SAVE ALL GO TERMS FUNCTION
# =====================================

save_all_go_terms <- function(go_results, comparison_name, cell_type, comparison_code, output_dir) {

  # Clean names for file saving
  comparison_clean <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
  cell_type_clean <- gsub("[^A-Za-z0-9_]", "_", cell_type)

  for (analysis_key in names(go_results)) {
    go_data <- go_results[[analysis_key]]

    if (is.null(go_data) || nrow(go_data@result) == 0) {
      next
    }

    # Parse analysis key
    parts <- strsplit(analysis_key, "_")[[1]]
    dir_val <- parts[1]
    go_cat <- parts[2]

    category_names <- list(
      "BP" = "Biological Process",
      "MF" = "Molecular Function",
      "CC" = "Cellular Component"
    )
    full_category_name <- category_names[[go_cat]]

    # Get ALL pathways (not just top 10)
    all_pathways <- go_data@result %>%
      arrange(p.adjust)

    n_pathways <- nrow(all_pathways)
    print(paste("  Saving", n_pathways, dir_val, go_cat, "pathways"))

    # =========================================
    # OUTPUT 1: Summary table (one row per pathway)
    # =========================================
    summary_table <- all_pathways
    summary_table$cell_type <- cell_type
    summary_table$comparison <- comparison_name
    summary_table$direction <- dir_val
    summary_table$go_category <- go_cat
    summary_table$go_full_name <- full_category_name
    summary_table$neg_log10_padj <- -log10(summary_table$p.adjust)
    summary_table$gene_symbols <- gsub("/", ", ", summary_table$geneID)

    summary_table <- summary_table %>%
      dplyr::select(
        cell_type, comparison, direction, go_category, go_full_name,
        ID, Description, GeneRatio, BgRatio,
        pvalue, p.adjust, qvalue, Count, neg_log10_padj,
        gene_symbols
      ) %>%
      arrange(p.adjust)

    summary_filename <- paste0(cell_type_clean, "_", comparison_clean, "_",
                               dir_val, "_", go_cat, "_ALL_Terms_Summary.csv")
    summary_filepath <- file.path(output_dir, summary_filename)

    tryCatch({
      write.csv(summary_table, summary_filepath, row.names = FALSE)
      print(paste("    Saved summary:", summary_filename))
    }, error = function(e) {
      print(paste("    Error saving summary:", e$message))
    })

    # =========================================
    # OUTPUT 2: Expanded table (one row per gene-pathway pair)
    # =========================================
    expanded_table <- all_pathways
    expanded_table$cell_type <- cell_type
    expanded_table$comparison <- comparison_name
    expanded_table$direction <- dir_val
    expanded_table$go_category <- go_cat
    expanded_table$go_full_name <- full_category_name

    expanded_table <- expanded_table %>%
      # Split geneID into separate rows
      separate_rows(geneID, sep = "/") %>%
      dplyr::rename(gene_symbol = geneID) %>%
      dplyr::select(
        cell_type, comparison, direction, go_category, go_full_name,
        ID, Description, p.adjust, Count,
        gene_symbol
      ) %>%
      arrange(p.adjust, gene_symbol)

    expanded_filename <- paste0(cell_type_clean, "_", comparison_clean, "_",
                                dir_val, "_", go_cat, "_ALL_Terms_Gene_Expanded.csv")
    expanded_filepath <- file.path(output_dir, expanded_filename)

    tryCatch({
      write.csv(expanded_table, expanded_filepath, row.names = FALSE)
      print(paste("    Saved expanded:", expanded_filename,
                  "(", nrow(expanded_table), "gene-pathway pairs)"))
    }, error = function(e) {
      print(paste("    Error saving expanded:", e$message))
    })
  }

  # =========================================
  # OUTPUT 3: Combined file with all GO categories
  # =========================================
  all_combined <- data.frame()

  category_names_map <- list(
    "BP" = "Biological Process",
    "MF" = "Molecular Function",
    "CC" = "Cellular Component"
  )

  for (analysis_key in names(go_results)) {
    go_data <- go_results[[analysis_key]]

    if (is.null(go_data) || nrow(go_data@result) == 0) next

    parts <- strsplit(analysis_key, "_")[[1]]
    dir_val2 <- parts[1]
    go_cat2 <- parts[2]

    temp_df <- go_data@result
    temp_df$cell_type <- cell_type
    temp_df$comparison <- comparison_name
    temp_df$direction <- dir_val2
    temp_df$go_category <- go_cat2
    temp_df$go_full_name <- category_names_map[[go_cat2]]
    temp_df$neg_log10_padj <- -log10(temp_df$p.adjust)
    temp_df$gene_symbols <- gsub("/", ", ", temp_df$geneID)

    all_combined <- rbind(all_combined, temp_df)
  }

  if (nrow(all_combined) > 0) {
    combined_filename <- paste0(cell_type_clean, "_", comparison_clean, "_ALL_GO_Terms_Combined.csv")
    combined_filepath <- file.path(output_dir, combined_filename)

    all_combined <- all_combined %>%
      dplyr::select(
        cell_type, comparison, direction, go_category, go_full_name,
        ID, Description, GeneRatio, BgRatio,
        pvalue, p.adjust, qvalue, Count, neg_log10_padj,
        gene_symbols
      ) %>%
      arrange(direction, go_category, p.adjust)

    tryCatch({
      write.csv(all_combined, combined_filepath, row.names = FALSE)
      print(paste("  Saved combined file:", combined_filename,
                  "(", nrow(all_combined), "total pathways)"))
    }, error = function(e) {
      print(paste("  Error saving combined:", e$message))
    })
  }
}

# =====================================
# MAIN EXECUTION
# =====================================

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print(paste("Created output directory:", output_dir))
}

# Find all DEG CSV files for Macrophages
deg_files <- list.files(deg_input_dir, pattern = "^Macrophages.*_DEG_results\\.csv$", full.names = TRUE)

print(paste(rep("=", 80), collapse = ""))
print("GO PATHWAY ANALYSIS - ALL TERMS FOR MACROPHAGES")
print(paste(rep("=", 80), collapse = ""))
print(paste("Input directory:", deg_input_dir))
print(paste("Output directory:", output_dir))
print("")
print("Target comparisons:")
for (comp in target_comparisons) {
  print(paste("  -", comp))
}
print("")

# Filter to target comparisons only
filtered_deg_files <- c()
for (deg_file in deg_files) {
  filename <- basename(deg_file)
  for (target in target_comparisons) {
    if (grepl(target, filename)) {
      filtered_deg_files <- c(filtered_deg_files, deg_file)
      break
    }
  }
}

print(paste("Found", length(filtered_deg_files), "matching DEG files:"))
for (file in filtered_deg_files) {
  print(paste("  ", basename(file)))
}

print(paste(rep("=", 80), collapse = ""))

# Process each DEG file
analysis_summary <- data.frame(
  File = character(),
  Cell_Type = character(),
  Comparison = character(),
  Up_BP = integer(),
  Up_MF = integer(),
  Up_CC = integer(),
  Down_BP = integer(),
  Down_MF = integer(),
  Down_CC = integer(),
  Total_Pathways = integer(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (deg_file in filtered_deg_files) {

  filename <- basename(deg_file)
  filename_parts <- strsplit(gsub(".csv$", "", filename), "_")[[1]]
  cell_type <- filename_parts[1]
  comparison_parts <- filename_parts[2:(length(filename_parts)-2)]
  comparison_name <- paste(comparison_parts, collapse = "_")

  tryCatch({
    result <- process_deg_file_all_terms(deg_file, output_dir)

    # Count pathways
    up_bp <- 0; up_mf <- 0; up_cc <- 0
    down_bp <- 0; down_mf <- 0; down_cc <- 0

    if (!is.null(result)) {
      if ("Up_BP" %in% names(result)) up_bp <- nrow(result[["Up_BP"]]@result)
      if ("Up_MF" %in% names(result)) up_mf <- nrow(result[["Up_MF"]]@result)
      if ("Up_CC" %in% names(result)) up_cc <- nrow(result[["Up_CC"]]@result)
      if ("Down_BP" %in% names(result)) down_bp <- nrow(result[["Down_BP"]]@result)
      if ("Down_MF" %in% names(result)) down_mf <- nrow(result[["Down_MF"]]@result)
      if ("Down_CC" %in% names(result)) down_cc <- nrow(result[["Down_CC"]]@result)
    }

    total <- up_bp + up_mf + up_cc + down_bp + down_mf + down_cc
    status <- ifelse(total == 0, "No significant pathways", "Completed")

    analysis_summary <- rbind(analysis_summary, data.frame(
      File = filename,
      Cell_Type = cell_type,
      Comparison = comparison_name,
      Up_BP = up_bp,
      Up_MF = up_mf,
      Up_CC = up_cc,
      Down_BP = down_bp,
      Down_MF = down_mf,
      Down_CC = down_cc,
      Total_Pathways = total,
      Status = status,
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    print(paste("Error processing", filename, ":", e$message))

    analysis_summary <- rbind(analysis_summary, data.frame(
      File = filename,
      Cell_Type = cell_type,
      Comparison = comparison_name,
      Up_BP = 0, Up_MF = 0, Up_CC = 0,
      Down_BP = 0, Down_MF = 0, Down_CC = 0,
      Total_Pathways = 0,
      Status = paste("Error:", e$message),
      stringsAsFactors = FALSE
    ))
  })
}

# Save analysis summary
write.csv(analysis_summary, file.path(output_dir, "GO_Analysis_Summary.csv"), row.names = FALSE)

# Print final summary
print(paste(rep("=", 80), collapse = ""))
print("ANALYSIS COMPLETED!")
print(paste(rep("=", 80), collapse = ""))
print("")
print("SUMMARY:")
print(analysis_summary)
print("")

# Count output files
csv_files <- list.files(output_dir, pattern = "\\.csv$")
print(paste("Generated", length(csv_files), "CSV files"))
print("")
print("Output files include:")
print("  - *_ALL_Terms_Summary.csv: One row per GO term with gene list")
print("  - *_ALL_Terms_Gene_Expanded.csv: One row per gene-pathway pair")
print("  - *_ALL_GO_Terms_Combined.csv: All GO categories in one file")
print("")
print(paste("Results saved to:", output_dir))
