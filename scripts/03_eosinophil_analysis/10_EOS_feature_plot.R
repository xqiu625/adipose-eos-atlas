library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
eos <- readRDS(eos_subset_file)

eos@meta.data$Clusters <- paste0("C", eos@meta.data$RNA_snn_res.0.8)

# Color by feature genes
# While there is overlap, microglia express some unique genes not found in macrophages, such as P2ry12, Tmem119, and Sall1
plot <- FeaturePlot(eos, 
            reduction = "umap", 
            features = c("Il4", "Ccr3"),
            ncol = 2)
plot_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_UMAP_plot_by_il4-ccr3_feature.png"
ggsave(plot_output_file, plot = plot, width = 8, height = 8)

