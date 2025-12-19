# cd /bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/03_EOS_GO_enrichment

# Load required libraries
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse annotation (change to org.Hs.eg.db for human)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

# Read the marker data
eos_markers <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8_all_markers.csv")

# Function to perform GO:BP enrichment for each cluster
perform_go_enrichment <- function(markers_df, pval_cutoff = 0.05, logfc_cutoff = 0.25) {
  
  # Get unique clusters
  clusters <- unique(markers_df$cluster)
  
  # Initialize list to store results
  go_results <- list()
  
  for(cluster in clusters) {
    cat("Processing cluster", cluster, "\n")
    
    # Filter genes for current cluster
    cluster_genes <- markers_df %>%
      filter(cluster == !!cluster, 
             p_val_adj < pval_cutoff,
             avg_log2FC > logfc_cutoff) %>%
      pull(gene)
    
    if(length(cluster_genes) < 5) {
      cat("Cluster", cluster, "has fewer than 5 significant genes, skipping...\n")
      next
    }
    
    # Convert gene symbols to Entrez IDs
    gene_ids <- bitr(cluster_genes, 
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db,
                     drop = TRUE)
    
    if(nrow(gene_ids) < 3) {
      cat("Cluster", cluster, "has fewer than 3 mappable genes, skipping...\n")
      next
    }
    
    # Perform GO:BP enrichment
    ego <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
    
    if(nrow(ego@result) > 0) {
      go_results[[paste0("Cluster_", cluster)]] <- ego@result
    }
  }
  
  return(go_results)
}

# Perform GO enrichment analysis
cat("Starting GO:BP enrichment analysis...\n")
go_enrichment_results <- perform_go_enrichment(eos_markers)

# Function to create heatmap data - exactly 6 terms per cluster
create_heatmap_data <- function(go_results, terms_per_cluster = 6) {
  
  # Collect exactly 6 top terms from each cluster
  all_terms <- list()
  
  for(cluster_name in names(go_results)) {
    cluster_results <- go_results[[cluster_name]]
    
    # Get exactly top 6 significant terms
    top_results <- cluster_results %>%
      arrange(p.adjust) %>%
      head(terms_per_cluster)
    
    cat("Cluster", cluster_name, "- selected", nrow(top_results), "terms\n")
    
    for(i in 1:nrow(top_results)) {
      term_id <- top_results$ID[i]
      term_desc <- top_results$Description[i]
      gene_ratio <- eval(parse(text = top_results$GeneRatio[i]))
      
      if(!term_id %in% names(all_terms)) {
        all_terms[[term_id]] <- list(
          description = term_desc,
          clusters = list()
        )
      }
      
      all_terms[[term_id]]$clusters[[cluster_name]] <- gene_ratio
    }
  }
  
  # Create matrix for heatmap - include ALL selected terms
  # Order cluster columns from 0 to 6 (as character: "Cluster_0", ..., "Cluster_6")
  cluster_order <- paste0("Cluster_", 0:6)
  cluster_names <- intersect(cluster_order, names(go_results))
  term_ids <- names(all_terms)
  
  # Create matrix with all terms (no filtering)
  heatmap_matrix <- matrix(0, 
                          nrow = length(all_terms), 
                          ncol = length(cluster_names),
                          dimnames = list(names(all_terms), cluster_names))
  
  # Fill matrix
  for(i in 1:nrow(heatmap_matrix)) {
    term_id <- rownames(heatmap_matrix)[i]
    term_data <- all_terms[[term_id]]
    
    for(cluster_name in cluster_names) {
      if(cluster_name %in% names(term_data$clusters)) {
        heatmap_matrix[i, cluster_name] <- term_data$clusters[[cluster_name]]
      }
    }
  }
  
  # Add term descriptions as row names
  term_descriptions <- sapply(all_terms, function(x) x$description)
  rownames(heatmap_matrix) <- term_descriptions
  
  return(heatmap_matrix)
}

# Create heatmap data - exactly 6 terms per cluster
cat("Creating heatmap data with 6 terms per cluster...\n")
heatmap_data <- create_heatmap_data(go_enrichment_results, terms_per_cluster = 6)

# Function to create the heatmap with gene ratio coloring (NO GRID LINES)
create_go_heatmap <- function(heatmap_matrix, filename = "GO_BP_enrichment_heatmap.pdf") {
  
  # Define color palette for gene ratio (white to dark red)
  colors <- colorRampPalette(c("white", "#FFE5E5", "#FFB3B3", "#FF8080", 
                              "#FF4D4D", "#FF1A1A", "#E60000", "#B30000", "#800000"))(100)
  
  # Print matrix range for color scale reference
  cat("Gene ratio range:", round(min(heatmap_matrix), 3), "to", round(max(heatmap_matrix), 3), "\n")
  
  # Create the heatmap: no clustering, clusters ordered 0-6, no grid lines
  p <- pheatmap(heatmap_matrix,
           color = colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 9,
           fontsize_col = 11,
           angle_col = 45,
           main = "GO:BP Enrichment Analysis\n(Gene Ratio - 6 Terms per Cluster)",
           legend_title = "Gene Ratio",
           filename = filename,
           width = 14,
           height = max(8, nrow(heatmap_matrix) * 0.3),
           cellwidth = 30,
           cellheight = 15,
           border_color = NA  # This removes all cell borders/grid lines
           )
  
  return(p)
}

# Generate the heatmap
cat("Generating heatmap with gene ratio coloring (no grid lines)...\n")
heatmap_plot <- create_go_heatmap(heatmap_data, "EOS_GO_BP_enrichment_heatmap.pdf")

# Print summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total clusters analyzed:", length(go_enrichment_results), "\n")
cat("Terms per cluster: 6\n")
cat("Total GO terms in heatmap:", nrow(heatmap_data), "\n")
cat("Clusters in heatmap:", ncol(heatmap_data), "\n")
cat("Gene ratio range:", round(min(heatmap_data), 3), "to", round(max(heatmap_data), 3), "\n")

# Display the matrix for inspection
cat("\nHeatmap data preview (showing gene ratios):\n")
print(round(heatmap_data[1:min(10, nrow(heatmap_data)), ], 3))

# Optional: Save detailed results
save_detailed_results <- function(go_results, filename = "GO_enrichment_detailed_results.xlsx") {
  require(openxlsx)
  
  wb <- createWorkbook()
  
  for(cluster_name in names(go_results)) {
    addWorksheet(wb, cluster_name)
    writeData(wb, cluster_name, go_results[[cluster_name]])
  }
  
  saveWorkbook(wb, filename, overwrite = TRUE)
  cat("Detailed results saved to:", filename, "\n")
}

# Uncomment to save detailed results
# save_detailed_results(go_enrichment_results)

cat("\nAnalysis complete! Check the generated PDF files for heatmaps.\n")