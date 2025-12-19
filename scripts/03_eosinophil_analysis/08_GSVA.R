lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
# BiocManager::install(version = "3.18", lib = lib_path)
# BiocManager::install("GSVA", lib = lib_path, dependencies = TRUE)

library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

genesets <- getGmt("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/public_data/mh.all.v2024.1.Mm.symbols.gmt")

# Load data and prepare
eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
eos <- readRDS(eos_subset_file)


# Extract expression data
eos.data <- GetAssayData(object = eos, layer = "data")

# Perform GSVA
gsvascore <- gsva(data.matrix(eos.data), genesets, parallel.sz = 2)
saveRDS(gsvascore, '/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_gsvascore.RDS')

# Prepare data frame for heatmap annotation
eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
eos <- readRDS(eos_subset_file)
eos@meta.data$Clusters <- paste0("C", eos@meta.data$RNA_snn_res.0.8)
DefaultAssay(eos) <- "RNA"


gsvascore <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_gsvascore.RDS")

eos.group <- data.frame(umi = rownames(eos@meta.data),
                       cluster = eos$Clusters,
                       stringsAsFactors = FALSE)

# Create a heatmap annotation based on cluster
ha.t <- HeatmapAnnotation(Cluster = eos.group$cluster)

# Generate heatmap
# Generate heatmap
plot <- Heatmap(as.matrix(gsvascore),
                show_column_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                top_annotation = ha.t,
                column_split = eos.group$cluster,
                row_names_gp = gpar(fontsize = 8),
                row_names_max_width = max_text_width(rownames(gsvascore), gp = gpar(fontsize = 8)))

# Save the plot
pdf("/bigdata/nairlab/shared/XQiu/EOS_project/figures/Eos_cluster_GSVA.pdf", width = 12, height = 8)
draw(plot)
dev.off()


