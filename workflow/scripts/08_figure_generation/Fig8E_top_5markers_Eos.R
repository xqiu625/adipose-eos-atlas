
library(Seurat)
library(dplyr)
library(ggplot2)

eos_markers <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8_all_markers.csv")
top_markers <- eos_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

eos_subset <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")

# Preparing the data for plotting
# Selecting features
features <- c("H2-Aa", "H2-Eb1", "H2-Ab1", "Fn1", "Ctss", "Itgb2", "H2-Q6", "Map2k3",
                 "Il4", "Stk40", "Ddx3y", "Eprs", "Cwc25", "Nop53", "Il4ra", "Mmp9",
                 "Krt80", "Pglyrp1", "Cited2", "P2ry10", "Vps37b", "Wnt11", "Cwc25", "Ptgs2",
                 "Jmy", "Polg2", "Iqcn", "Gstp3", "Tmem230", "Plekhg2", "Dand5", "Tnfrsf25",
                 "Spata25", "Hes2", "Muc2")

# Assuming you have top_markers$gene
gene_list <- top_markers$gene

unique_genes <- unique(gene_list)

# Create the actual vector
features <- unique_genes

# If you want to see the vector in the format you specified:
cat('features <- c(', paste(sprintf('"%s"', features), collapse = ", "), ')\n')

# To use the vector:
print(features)

# Extract expression data
expression_data <- GetAssayData(eos_subset, assay = "RNA", layer = "data")[features, ]

# Update Seurat cluster identities in metadata
eos_subset@meta.data$Clusters <- factor(paste0("C", eos_subset@meta.data$RNA_snn_res.0.8),
                                               levels = paste0("C", 0:6))

# Calculate average expression per cluster
avg_expression <- AverageExpression(eos_subset, features = features, group.by = "Clusters")$RNA

# Convert sparse matrix to regular matrix
avg_expression <- as.matrix(avg_expression)
avg_expression <- avg_expression[features, ]

# Scale the data
scaled_data <- t(scale(t(avg_expression)))

# Create the heatmap
library(pheatmap)

# Define a color palette
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)



# Save the plot
pdf("/bigdata/nairlab/shared/XQiu/EOS_project/figures/grant_figures/Eos_Heatmap_clusterMarker_0.8_top10.pdf", 
    width = 10, 
    height = 24)

heatmap <- pheatmap(scaled_data,
                    color = color_palette,
                    cluster_rows = FALSE,  # Disable row clustering to maintain order
                    cluster_cols = FALSE,
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    main = "Gene Expression Heatmap",
                    fontsize = 14,  # Increased base font size
                    fontsize_row = 12,  # Increased row font size
                    fontsize_col = 12,  # Increased column font size
                    angle_col = 45,
                    cellwidth = 15,  # Adjust cell width to make heatmap denser
                    cellheight = 15,  # Adjust cell height to make heatmap denser
                    border_color = NA,  # Remove cell borders for a denser look
                    legend = TRUE,
                    legend_breaks = c(min(scaled_data), 0, max(scaled_data)),
                    legend_labels = c("Low", "Medium", "High"),
                    treeheight_row = 0,  # Remove row dendrogram
                    treeheight_col = 0)  # Remove column dendrogram

dev.off()
