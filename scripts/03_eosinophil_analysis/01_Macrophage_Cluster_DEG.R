# === Load packages ===
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ReactomePA)
library(ggplot2)

# === Define input/output ===
seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed/seurat_object_labeled.rds"
out_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R06_20250606/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# === Load data and subset Macrophages ===
seurat_object <- readRDS(seurat_file)
macrophages <- subset(seurat_object, subset = cell_type == "Macrophages")

# === Re-process Macrophages only ===
macrophages <- NormalizeData(macrophages)
macrophages <- FindVariableFeatures(macrophages, selection.method = "vst", nfeatures = 2000)
macrophages <- ScaleData(macrophages)
macrophages <- RunPCA(macrophages)
macrophages <- RunUMAP(macrophages, dims = 1:20)
macrophages <- FindNeighbors(macrophages, dims = 1:20)
macrophages <- FindClusters(macrophages, resolution = 0.15)  # ~4–6 clusters

# Save UMAP plot
umap_plot <- DimPlot(macrophages, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Macrophage Subclusters")
ggsave(paste0(out_dir, "macrophage_umap.pdf"), umap_plot, width = 7, height = 6)

# === Differential Expression Analysis ===
deg <- FindAllMarkers(macrophages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(deg, file = paste0(out_dir, "macrophage_cluster_DEGs.csv"), row.names = FALSE)

# === 输出目录 ===
dir.create(file.path(out_dir, "deg_by_cluster"), showWarnings = FALSE)
dir.create(file.path(out_dir, "violin_by_cluster"), showWarnings = FALSE)


# === 设置 ID 为你想分析的分群 ===
Idents(macrophages) <- "RNA_snn_res.0.15"

# === 按 cluster 拆分 DEG 并保存 top gene CSV 与 violin 图 ===
deg <- read.csv(paste0(out_dir, "macrophage_cluster_DEGs.csv"))
clusters <- unique(deg$cluster)

for (clust in clusters) {
  deg_sub <- deg %>% filter(cluster == clust) %>% arrange(desc(avg_log2FC))
  
  # Save CSV
  write.csv(deg_sub, file = file.path(out_dir, "deg_by_cluster", paste0("cluster", clust, "_DEGs.csv")), row.names = FALSE)
  
  # Top 5 genes violin plot
  top5 <- head(deg_sub$gene, 5)
  vln_plot <- VlnPlot(macrophages, features = top5, group.by = "RNA_snn_res.0.15", pt.size = 0.1, ncol = 1) +
    ggtitle(paste("Cluster", clust, "- Top 5 DEGs"))
  ggsave(file.path(out_dir, "violin_by_cluster", paste0("cluster", clust, "_top5_violin.pdf")), vln_plot, height = 10)
}

# === Top 20 genes heatmap（所有 cluster）===
top20_genes <- deg %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

heatmap <- DoHeatmap(macrophages, features = top20_genes, group.by = "RNA_snn_res.0.15") +
  ggtitle("Top 20 DEGs per Cluster")
ggsave(file.path(out_dir, "macrophage_top20_heatmap.pdf"), heatmap, width = 10, height = 12)

# === Top 10 genes dotmap（所有 cluster）===
# Get top 10 genes per cluster
top10_genes <- deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

# Create and save dotplot
dotplot <- DotPlot(macrophages, features = top10_genes, group.by = "RNA_snn_res.0.15") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 DEGs per Cluster")
ggsave(file.path(out_dir, "macrophage_top10_dotplot.pdf"), dotplot, width = 16, height = 8)


