# (for WT female v WT female Inf) and (KO female vs KO female inf).
library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)

eos_subset <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")


# 2. Set up the variables for comparison
eos_subset$Sex <- factor(eos_subset$Sex)
eos_subset$condition <- factor(eos_subset$condition)
eos_subset$infection_status <- factor(eos_subset$infection_status)
eos_subset$genotype <- factor(eos_subset$genotype)

eos_subset_f <- subset(eos_subset, Sex == "Female")
eos_subset_f_WT <- subset(eos_subset_f, genotype == "WT")
eos_subset_f_KO <- subset(eos_subset_f, genotype == "KO")

# 3. Perform differential expression analysis for each variable
Idents(eos_subset_f_WT) <- "infection_status"
deg_wt_f_inf <- FindMarkers(eos_subset_f_WT, ident.1 = "Infected", ident.2 = "Not_Infected", test.use = "MAST")

Idents(eos_subset_f_KO) <- "infection_status"
deg_ko_f_inf <- FindMarkers(eos_subset_f_KO, ident.1 = "Infected", ident.2 = "Not_Infected", test.use = "MAST")


# 4. Add gene names and sort results
add_gene_names <- function(deg_results) {
  deg_results$gene <- rownames(deg_results)
  deg_results <- deg_results %>% arrange(p_val_adj)
  return(deg_results)
}

deg_wt_f_inf <- add_gene_names(deg_wt_f_inf)
deg_ko_f_inf <- add_gene_names(deg_ko_f_inf)



# 5. Save results
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data")
write.csv(deg_wt_f_inf, "DEG_eosinophils_WT_female_v_WT_female_Inf.csv", row.names = FALSE)
write.csv(deg_ko_f_inf, "DEG_eosinophils_KO_female_vs_KO_female_inf.csv", row.names = FALSE)


# 5. enrichment

# Function to perform enrichment analysis dot plot and cnet plot
# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(enrichplot)

# Function to perform enrichment analysis dot plot and emap plot
perform_enrichment <- function(deg_results, comparison_name) {
  # Check if required packages are installed and loaded
  required_packages <- c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "dplyr", "ggplot2", "enrichplot")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it using install.packages('", pkg, "')"))
    }
    library(pkg, character.only = TRUE)
  }

  # Function to process genes and perform enrichment
  process_genes <- function(genes, direction) {
    if (length(genes) == 0) {
      message(paste("No", direction, "genes found for", comparison_name))
      return(NULL)
    }
    
    # Convert gene symbols to ENTREZ IDs
    entrez_ids <- tryCatch({
      mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
    }, error = function(e) {
      stop(paste("Error in mapping genes to ENTREZ IDs:", e$message))
    })
    
    # Check for genes that didn't map to ENTREZ IDs
    unmapped_genes <- genes[is.na(entrez_ids)]
    if (length(unmapped_genes) > 0) {
      message(paste("Warning:", length(unmapped_genes), "genes could not be mapped to ENTREZ IDs for", direction, "genes in", comparison_name))
      message("First few unmapped genes: ", paste(head(unmapped_genes), collapse = ", "))
    }
    
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    if (length(entrez_ids) == 0) {
      message(paste("No valid ENTREZ IDs found for", direction, "genes in", comparison_name))
      return(NULL)
    }
    
    # GO Enrichment (Biological Process)
    go_bp <- tryCatch({
      enrichGO(gene = entrez_ids,
               OrgDb = org.Mm.eg.db,
               ont = "BP",
               pvalueCutoff = 0.05)
    }, error = function(e) {
      message(paste("GO BP enrichment failed for", direction, "genes in", comparison_name, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      write.csv(as.data.frame(go_bp), paste0("GO_BP_enrichment_", comparison_name, "_", direction, ".csv"), row.names = FALSE)
      
      # Dotplot
      pdf(paste0("GO_BP_dotplot_", comparison_name, "_", direction, ".pdf"), width = 10, height = 16)
      print(dotplot(go_bp, showCategory = 20, title = paste("GO Biological Process -", comparison_name, "-", direction)))
      dev.off()
      
    } else {
      message(paste("No significant GO BP terms found for", direction, "genes in", comparison_name))
    }
    
    return(go_bp)
  }
  
  # Select upregulated and downregulated genes
  upregulated_genes <- deg_results %>%
    filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
    pull(gene)
  
  downregulated_genes <- deg_results %>%
    filter(avg_log2FC < -0.25 & p_val_adj < 0.05) %>%
    pull(gene)
  
  # Perform enrichment for upregulated and downregulated genes
  up_enrichment <- process_genes(upregulated_genes, "upregulated")
  down_enrichment <- process_genes(downregulated_genes, "downregulated")
  
  return(list(upregulated = up_enrichment, downregulated = down_enrichment))
}
# Perform enrichment for each comparison
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/figures/grant_figures/enrichment")
enrichment_results_1 <- perform_enrichment(deg_wt_f_inf, "WT_female_Inf_vs_NonInf")
enrichment_results_2 <- perform_enrichment(deg_ko_f_inf, "KO_female_Inf_vs_NonInf")






