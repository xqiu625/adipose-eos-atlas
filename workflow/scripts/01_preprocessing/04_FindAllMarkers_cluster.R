
library(Seurat)
library(ggplot2)
library(dplyr)

seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_singleR_clustered.rds"
seurat_object <- readRDS(file = seurat_file)

seurat_object <- JoinLayers(seurat_object)

# Ensure the Seurat object is properly set up
seurat_object <- SetIdent(seurat_object, value = "seurat_clusters")

# Find markers for all clusters
all_markers <- FindAllMarkers(seurat_object, 
                              only.pos = TRUE,  # Only positive markers
                              min.pct = 0.25,   # Expressed in at least 25% of cells
                              logfc.threshold = 0.25)  # At least 0.25 log2 fold change

# Sort the markers by adjusted p-value
all_markers <- all_markers[order(all_markers$cluster, all_markers$p_val_adj),]

# Get top 10 markers per cluster
top10 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Print the top 10 markers for each cluster
print(top10)

# Save the results to a CSV file
write.csv(all_markers, "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/all_cluster_markers.csv", row.names = FALSE)
write.csv(top10, "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/top10_cluster_markers.csv", row.names = FALSE)

# Visualize top markers
library(ggplot2)

top5 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

plot <- ggplot(top5, aes(x = avg_log2FC, y = gene, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cluster, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(x = "Average log2 Fold Change", y = "Gene", title = "Top 5 Marker Genes per Cluster")
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/top5_markers_per_cluster.png"
ggsave(plot_output_file, plot,width = 15, height = 20, units = "in", dpi = 300)


