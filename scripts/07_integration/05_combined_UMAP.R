save_path <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/"
filename <- "harmony_integrated_eos_with_gurtner_posner_jorssen.rds"
full_path <- file.path(save_path, filename)

# Load required libraries
library(Seurat)
library(ggplot2)

# Read the RDS file
merged_seurat <- readRDS(file = full_path)

# Create UMAP plots based on sex, condition, organ, and source_study
umap_sex <- DimPlot(merged_seurat, reduction = "umap", group.by = "sex", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Merged Data Colored by Sex")

umap_condition <- DimPlot(merged_seurat, reduction = "umap", group.by = "condition", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Merged Data Colored by Condition")

umap_organ <- DimPlot(merged_seurat, reduction = "umap", group.by = "organ", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Merged Data Colored by Organ")

umap_source_study <- DimPlot(merged_seurat, reduction = "umap", group.by = "source_study", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of Merged Data Colored by Source Study")

# Save the plots
save_path <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/combined_published_data/"
ggsave(filename = file.path(save_path, "umap_by_sex.pdf"), plot = umap_sex, width = 12, height = 10)
ggsave(filename = file.path(save_path, "umap_by_condition.pdf"), plot = umap_condition, width = 12, height = 10)
ggsave(filename = file.path(save_path, "umap_by_organ.pdf"), plot = umap_organ, width = 12, height = 10)
ggsave(filename = file.path(save_path, "umap_by_source_study.pdf"), plot = umap_source_study, width = 12, height = 10)


# Create separate UMAP plots for each organ
organs <- unique(merged_seurat$organ)
for (organ in organs) {
  umap_organ_specific <- DimPlot(merged_seurat, reduction = "umap", group.by = "organ", label = TRUE, repel = TRUE, cells.highlight = WhichCells(merged_seurat, expression = organ == !!organ)) +
    ggtitle(paste("UMAP of Merged Data for", organ))
  
  # Save the plot
  ggsave(filename = file.path(save_path, paste0("umap_", organ, ".pdf")), plot = umap_organ_specific, width = 12, height = 10)
}


# Create separate UMAP plots for each condition
conditions <- unique(merged_seurat$condition)
for (condition in conditions) {
  umap_condition_specific <- DimPlot(merged_seurat, reduction = "umap", group.by = "condition", label = TRUE, repel = TRUE, cells.highlight = WhichCells(merged_seurat, expression = condition == !!condition)) +
    ggtitle(paste("UMAP of Merged Data for", condition)) +
    scale_color_manual(values = c("#7f7f7f","#002fa7", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",  "#bcbd22", "#17becf"))
  
  # Save the plot
  ggsave(filename = file.path(save_path, paste0("umap_", condition, ".pdf")), plot = umap_condition_specific, width = 12, height = 10)
}

# Create separate UMAP plots for each source study
source_studies <- unique(merged_seurat$source_study)
for (source_study in source_studies) {
  umap_source_study_specific <- DimPlot(merged_seurat, reduction = "umap", group.by = "source_study", label = TRUE, repel = TRUE, cells.highlight = WhichCells(merged_seurat, expression = source_study == !!source_study)) +
    ggtitle(paste("UMAP of Merged Data for", source_study)) +
        scale_color_manual(values = c("#7f7f7f","#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",  "#bcbd22", "#17becf"))
  
  # Save the plot
  ggsave(filename = file.path(save_path, paste0("umap_", source_study, ".pdf")), plot = umap_source_study_specific, width = 12, height = 10)
}

# Create separate UMAP plots for each cluster
clusters <- unique(merged_seurat$cluster)
for (cluster in clusters) {
  umap_cluster_specific <- DimPlot(merged_seurat, reduction = "umap", group.by = "cluster", label = TRUE, repel = TRUE, cells.highlight = WhichCells(merged_seurat, expression = cluster == !!cluster)) +
    ggtitle(paste("UMAP of Merged Data for Cluster", cluster)) +
    scale_color_manual(values = c("#7f7f7f","#9467bd", "#8c564b", "#e377c2",  "#bcbd22", "#17becf"))
  
  # Save the plot
  ggsave(filename = file.path(save_path, paste0("umap_cluster_", cluster, ".pdf")), plot = umap_cluster_specific, width = 12, height = 10)
}

