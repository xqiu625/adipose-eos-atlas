# Load required additional libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # For mouse genes. Use org.Hs.eg.db for human genes
library(enrichplot)
library(dplyr)

setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data")
# Function to perform enrichment analysis
perform_enrichment <- function(deg_results, comparison_name) {
  # Function to process genes and perform enrichment
  process_genes <- function(genes, direction) {
    if (length(genes) == 0) {
      message(paste("No", direction, "genes found for", comparison_name))
      return(NULL)
    }
    
    # Convert gene symbols to ENTREZ IDs
    entrez_ids <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
    
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
    
    # KEGG Pathway Enrichment
    kegg_enrichment <- tryCatch({
      enrichKEGG(gene = entrez_ids,
                 organism = 'mmu',  # 'mmu' for mouse, 'hsa' for human
                 pvalueCutoff = 0.05)
    }, error = function(e) {
      message(paste("KEGG enrichment failed for", direction, "genes in", comparison_name, ":", e$message))
      return(NULL)
    })
    
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
    
    # Save results (only if enrichment was successful)
    if (!is.null(kegg_enrichment) && nrow(kegg_enrichment) > 0) {
      write.csv(as.data.frame(kegg_enrichment), paste0("KEGG_enrichment_", comparison_name, "_", direction, ".csv"), row.names = FALSE)
      pdf(paste0("KEGG_dotplot_", comparison_name, "_", direction, ".pdf"), width = 10, height = 16)
      print(dotplot(kegg_enrichment, showCategory = 20, title = paste("KEGG Pathways -", comparison_name, "-", direction)))
      dev.off()
    } else {
      message(paste("No significant KEGG pathways found for", direction, "genes in", comparison_name))
    }
    
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      write.csv(as.data.frame(go_bp), paste0("GO_BP_enrichment_", comparison_name, "_", direction, ".csv"), row.names = FALSE)
      pdf(paste0("GO_BP_dotplot_", comparison_name, "_", direction, ".pdf"), width = 10, height = 16)
      print(dotplot(go_bp, showCategory = 20, title = paste("GO Biological Process -", comparison_name, "-", direction)))
      dev.off()
    } else {
      message(paste("No significant GO BP terms found for", direction, "genes in", comparison_name))
    }
    
    return(list(kegg = kegg_enrichment, go_bp = go_bp))
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
enrichment_results_1 <- perform_enrichment(deg_condition_1, "KO_Female_vs_KO_Male_HFD")
enrichment_results_2 <- perform_enrichment(deg_condition_2, "WT_Female_vs_WT_Male_HFD")
enrichment_results_3 <- perform_enrichment(deg_condition_3, "KO_Female_vs_WT_Female_HFD")
enrichment_results_4 <- perform_enrichment(deg_condition_4, "KO_Male_vs_WT_Male_HFD")

# Perform enrichment for other comparisons
deg_sex <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/DEG_eosinophils_sex.csv")
deg_infected <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/DEG_eosinophils_infected.csv")
deg_genotype <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/DEG_eosinophils_genotype.csv")

enrichment_results_sex <- perform_enrichment(deg_sex, "Sex")
enrichment_results_infected <- perform_enrichment(deg_infected, "Infection_Status")
enrichment_results_genotype <- perform_enrichment(deg_genotype, "Genotype")
