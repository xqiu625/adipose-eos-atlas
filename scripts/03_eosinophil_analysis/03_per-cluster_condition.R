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

# === Save the processed macrophage object with clusters ===
saveRDS(macrophages, file = file.path(out_dir, "macrophages_clustered.rds"))

# === Create cross-tabulations ===
# Table 1: cluster -> condition
cluster_condition_table <- table(macrophages$seurat_clusters, macrophages$condition)
write.csv(as.data.frame.matrix(cluster_condition_table),
          file = file.path(out_dir, "macrophage_cluster_by_condition.csv"))

# Table 2: condition -> cluster
condition_cluster_table <- table(macrophages$condition, macrophages$seurat_clusters)
write.csv(as.data.frame.matrix(condition_cluster_table),
          file = file.path(out_dir, "macrophage_condition_by_cluster.csv"))


# === Plot stacked barplot: cells per condition by cluster ===
df2 <- as.data.frame(table(macrophages$condition, macrophages$seurat_clusters))
colnames(df2) <- c("Condition", "Cluster", "Count")
df2$Condition <- factor(df2$Condition)
df2$Cluster <- factor(df2$Cluster, levels = sort(as.numeric(as.character(unique(df2$Cluster)))))

p2 <- ggplot(df2, aes(x = Condition, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of Cells") +
  xlab("Condition") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  ggtitle("Conditions by Macrophage Cluster (Stacked Proportion Barplot)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(out_dir, "macrophage_condition_by_cluster_barplot.pdf"),
       plot = p2, width = 8, height = 6)


# === Load required package for plotting ===
library(ggplot2)

# === Prepare data for stacked barplot ===
df <- as.data.frame(table(macrophages$seurat_clusters, macrophages$condition))
colnames(df) <- c("Cluster", "Condition", "Count")

# Ensure Cluster is a factor with numeric ordering
df$Cluster <- factor(df$Cluster, levels = sort(as.numeric(as.character(unique(df$Cluster)))))

# === Plot stacked barplot: cells per cluster by condition ===
p <- ggplot(df, aes(x = Cluster, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of Cells") +
  xlab("Macrophage Cluster") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  ggtitle("Macrophage Clusters by Condition (Stacked Proportion Barplot)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === Save plot ===
ggsave(filename = file.path(out_dir, "macrophage_cluster_condition_barplot.pdf"),
       plot = p, width = 8, height = 6)


# === Plot stacked barplot: cells per condition by cluster ===
df2 <- as.data.frame(table(macrophages$condition, macrophages$seurat_clusters))
colnames(df2) <- c("Condition", "Cluster", "Count")
df2$Condition <- factor(df2$Condition)
df2$Cluster <- factor(df2$Cluster, levels = sort(as.numeric(as.character(unique(df2$Cluster)))))

p2 <- ggplot(df2, aes(x = Condition, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Proportion of Cells") +
  xlab("Condition") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  ggtitle("Conditions by Macrophage Cluster (Stacked Proportion Barplot)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(out_dir, "macrophage_condition_by_cluster_barplot.pdf"),
       plot = p2, width = 8, height = 6)

