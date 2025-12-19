# === Packages ===
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(dplyr)

# === 创建输出文件夹 ===
dir.create(file.path(out_dir, "pathway_by_cluster"), showWarnings = FALSE)

# === 定义输出目录 ===
out_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R06_20250606/"

# === 从已保存的 DEG 结果 CSV 中读取 ===
deg <- read.csv(file.path(out_dir, "macrophage_cluster_DEGs.csv"))

# === 获取 cluster 列表 ===
clusters <- unique(deg$cluster)

for (clust in clusters) {
  message("Processing cluster ", clust, "...")
  
  # Subset DEG
  deg_sub <- deg %>% filter(cluster == clust) %>% arrange(desc(avg_log2FC))
  gene_symbols <- deg_sub$gene
  
  # Map SYMBOL -> ENTREZID
  mapping <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  gene_entrez <- unique(mapping$ENTREZID)
  
  # === Enrichment: GO BP ===
  ego <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  write.csv(as.data.frame(ego), file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_GO_BP.csv")), row.names = FALSE)
  saveRDS(ego, file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_GO_BP.rds")))

  # === Enrichment: KEGG ===
  ekegg <- enrichKEGG(gene = gene_entrez,
                      organism = "mmu",
                      pvalueCutoff = 0.05)
  write.csv(as.data.frame(ekegg), file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_KEGG.csv")), row.names = FALSE)
  saveRDS(ekegg, file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_KEGG.rds")))

  # === Enrichment: REACTOME ===
  epath <- enrichPathway(gene = gene_entrez,
                         organism = "mouse",
                         pvalueCutoff = 0.05,
                         readable = TRUE)
  write.csv(as.data.frame(epath), file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_REACTOME.csv")), row.names = FALSE)
  saveRDS(epath, file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_REACTOME.rds")))

  # === GSEA input: ranked gene list ===
  gsea_ranks <- deg_sub %>%
    left_join(mapping, by = c("gene" = "SYMBOL")) %>%
    dplyr::select(ENTREZID, avg_log2FC) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    arrange(desc(avg_log2FC))

  gene_list <- gsea_ranks$avg_log2FC
  names(gene_list) <- gsea_ranks$ENTREZID

  egsea <- gseGO(geneList = gene_list,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 keyType = "ENTREZID",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)
  write.csv(as.data.frame(egsea), file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_GSEA_GO_BP.csv")), row.names = FALSE)
  saveRDS(egsea, file = file.path(out_dir, "pathway_by_cluster", paste0("cluster", clust, "_GSEA_GO_BP.rds")))
}




library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# === 设置路径 ===
path_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R06_20250606/pathway_by_cluster"
plot_dir <- file.path(path_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# === cluster 列表 ===
clusters <- 0:5

for (clust in clusters) {
  message("Generating barplots for cluster ", clust)

  # === 文件路径 ===
  rds_files <- list(
    GO_BP = file.path(path_dir, paste0("cluster", clust, "_GO_BP.rds")),
    KEGG = file.path(path_dir, paste0("cluster", clust, "_KEGG.rds")),
    REACTOME = file.path(path_dir, paste0("cluster", clust, "_REACTOME.rds")),
    GSEA = file.path(path_dir, paste0("cluster", clust, "_GSEA_GO_BP.rds"))
  )

  # === 绘图函数 ===
  for (type in names(rds_files)) {
    rds_path <- rds_files[[type]]
    if (file.exists(rds_path)) {
      enr_obj <- readRDS(rds_path)
      # Check for empty enrichment results
      if (nrow(as.data.frame(enr_obj)) == 0) {
        message(paste("No enrichment results for cluster", clust, type, "- skipping."))
        next
      }
      out_pdf <- file.path(plot_dir, paste0("cluster", clust, "_", type, "_barplot.pdf"))
      
      # Create plot with adjusted parameters for better readability
      p <- enrichplot:::barplot.enrichResult(enr_obj, showCategory = 15) +
        ggtitle(paste0(type, " - Cluster ", clust)) +
        theme(
          axis.text.y = element_text(size = 8, face = "plain"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 120, unit = "pt")
        )
      
      # Calculate height based on number of categories
      n_categories <- min(15, nrow(as.data.frame(enr_obj)))
      plot_height <- max(6, n_categories * 0.4)  # Minimum height of 6, scales with categories
      
      pdf(out_pdf, width = 12, height = plot_height)  # Increased width for better label visibility
      print(p)
      dev.off()
    }
  }
}