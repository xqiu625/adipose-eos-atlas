# === Load packages ===
library(Seurat)
library(dplyr)
library(ggplot2)

# === Define input/output ===
macrophages_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/R06_20250606/macrophages_clustered.rds"
out_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R06_20250606/"

# Load the macrophages object
macrophages <- readRDS(macrophages_file)

# === Custom marker genes dot plot ===
# Define the genes of interest
marker_genes <- c("Plac8", "Lyz1", "Lyz2", "Socs6", "Fcgr1", "Ly6c2", "Ccr2", 
                  "Ace", "Dgkg", "Retnla", "Lyve1", "Cd163", "Cd209f", "Cd9", 
                  "Trem2", "Lpl", "Mrc1", "Adgre1", "Maoa", "St3gal6", "Lilra5", 
                  "Siglecf", "Flt3", "Ccr7", "Clec9a", "Itgax", "Cd64", 
                  "Il1b", "Tnf", "Lgals3")

# Check which genes are present in the dataset
available_genes <- marker_genes[marker_genes %in% rownames(macrophages)]
missing_genes <- marker_genes[!marker_genes %in% rownames(macrophages)]

# Print information about gene availability
cat("Available genes:", length(available_genes), "out of", length(marker_genes), "\n")
if(length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

# Create dot plot with available genes
marker_dotplot <- DotPlot(macrophages, 
                         features = available_genes, 
                         group.by = "RNA_snn_res.0.15") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle("Macrophage Marker Genes Expression") +
  xlab("Genes") +
  ylab("Clusters")

# Save the plot
ggsave(file.path(out_dir, "macrophage_marker_genes_dotplot.pdf"), 
       marker_dotplot, width = 18, height = 8)

# Also save as PNG for easier viewing
ggsave(file.path(out_dir, "macrophage_marker_genes_dotplot.png"), 
       marker_dotplot, width = 18, height = 8, dpi = 300)

# Optional: Create a violin plot for these marker genes as well
marker_violin <- VlnPlot(macrophages, 
                        features = available_genes, 
                        group.by = "RNA_snn_res.0.15", 
                        pt.size = 0, 
                        ncol = 5) +
  ggtitle("Macrophage Marker Genes Expression")

ggsave(file.path(out_dir, "macrophage_marker_genes_violin.pdf"), 
       marker_violin, width = 20, height = 16)

# Print summary
cat("Dot plot saved as:", file.path(out_dir, "macrophage_marker_genes_dotplot.pdf"), "\n")
cat("Violin plot saved as:", file.path(out_dir, "macrophage_marker_genes_violin.pdf"), "\n")


# === Check the current cluster information ===
cat("Current clusters:", unique(macrophages$seurat_clusters), "\n")
cat("Cluster class:", class(macrophages$seurat_clusters), "\n")

# === Create simplified cluster labels ===
# Convert seurat_clusters to character first, then map
macrophages$simplified_label <- case_when(
  as.character(macrophages$seurat_clusters) %in% c("0", "1") ~ "Mac_PV",
  as.character(macrophages$seurat_clusters) == "2" ~ "Mac_ISG", 
  as.character(macrophages$seurat_clusters) %in% c("3", "4") ~ "Mac_LAMinflam",
  as.character(macrophages$seurat_clusters) == "5" ~ "Mac_Cycling",
  TRUE ~ paste0("Cluster_", macrophages$seurat_clusters)  # fallback for any unexpected clusters
)

# Verify the mapping worked
cat("Simplified label mapping:\n")
table(macrophages$seurat_clusters, macrophages$simplified_label)

# === Create UMAP plot with simplified labels ===
simplified_umap <- DimPlot(macrophages, 
                          group.by = "simplified_label",
                          label = TRUE,
                          repel = TRUE) +
  ggtitle("Macrophage Clusters - Simplified Labels") +
  theme_minimal()

# Save the plot
ggsave(file.path(out_dir, "macrophage_simplified_labels_umap.pdf"),
       simplified_umap, width = 10, height = 8)

# Also save as PNG
ggsave(file.path(out_dir, "macrophage_simplified_labels_umap.png"),
       simplified_umap, width = 10, height = 8, dpi = 300)

# Print completion message
cat("UMAP with simplified labels saved as:", 
    file.path(out_dir, "macrophage_simplified_labels_umap.pdf"), "\n")






# === Custom marker genes dot plot ===
# Define the genes of interest
marker_genes <- c("Plac8", "Lyz1", "Lyz2", "Socs6", "Fcgr1", "Ly6c2", "Ccr2", 
                  "Ace", "Dgkg", "Retnla", "Lyve1", "Cd163", "Cd209f", "Cd9", 
                  "Trem2", "Lpl", "Mrc1", "Adgre1", "Maoa", "St3gal6", "Lilra5", 
                  "Siglecf", "Flt3", "Ccr7", "Clec9a", "Itgax", "Cd64", 
                  "Il1b", "Tnf", "Lgals3")

# Check which genes are present in the dataset
available_genes <- marker_genes[marker_genes %in% rownames(macrophages)]
missing_genes <- marker_genes[!marker_genes %in% rownames(macrophages)]

# Print information about gene availability
cat("Available genes:", length(available_genes), "out of", length(marker_genes), "\n")
if(length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

# Create dot plot with available genes
marker_dotplot <- DotPlot(macrophages, 
                         features = available_genes, 
                         group.by = "simplified_label") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle("Macrophage Marker Genes Expression") +
  xlab("Genes") +
  ylab("Clusters")

# Save the plot
ggsave(file.path(out_dir, "macrophage_marker_genes_dotplot_simplified_label.pdf"), 
       marker_dotplot, width = 18, height = 8)