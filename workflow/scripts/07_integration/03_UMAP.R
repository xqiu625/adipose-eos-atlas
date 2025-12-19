# Visualize the results
# Plot UMAP and color by orig.ident
library(Seurat)
library(ggplot2)
library(cowplot)

# Set working directory
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/")

# Function to create UMAP plot
create_umap_plot <- function(seurat_obj, group.by = "orig.ident", title = NULL) {
  plot_data <- data.frame(
    UMAP_1 = Embeddings(seurat_obj, "umap")[,1],
    UMAP_2 = Embeddings(seurat_obj, "umap")[,2],
    Group = seurat_obj@meta.data[[group.by]]
  )
  
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_minimal() +  # Use minimal theme as a base
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  
  return(p)
}

# Create UMAP plot
umap_plot <- create_umap_plot(merged_seurat, group.by = "orig.ident", title = "UMAP by Original Identity")

# Save plot
ggsave("imm_UMAP_plot_by_orig.ident.png", plot = umap_plot, width = 10, height = 8, units = "in", dpi = 300)



# Function to create UMAP plot with white background
create_umap_plot <- function(seurat_obj, group.by = "Clusters", title = NULL) {
  plot_data <- data.frame(
    UMAP_1 = Embeddings(seurat_obj, "umap")[,1],
    UMAP_2 = Embeddings(seurat_obj, "umap")[,2],
    Group = seurat_obj@meta.data[[group.by]]
  )
  
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_minimal() +  # Use minimal theme as a base
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  
  # Add labels for clusters
  cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ Group, data = plot_data, FUN = mean)
  p <- p + geom_text(data = cluster_centers, aes(label = Group), size = 4, fontface = "bold", color = "black")
  
  return(p)
}

# Create UMAP plot
umap_plot <- create_umap_plot(merged_seurat, group.by = "Clusters", title = "UMAP by Clusters")

# Save plot
ggsave("imm_UMAP_plot_by_Clusters.png", plot = umap_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")





library(Seurat)
library(ggplot2)
library(cowplot)

# Set working directory
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/")

# Function to create UMAP plot with custom color for GSM7919060
create_umap_plot <- function(seurat_obj, group.by = "condition", title = NULL) {
  plot_data <- data.frame(
    UMAP_1 = Embeddings(seurat_obj, "umap")[,1],
    UMAP_2 = Embeddings(seurat_obj, "umap")[,2],
    Group = seurat_obj@meta.data[[group.by]]
  )
  
  # Get unique conditions
  conditions <- unique(plot_data$Group)
  
  # Create a color palette
  color_palette <- scales::hue_pal()(length(conditions))
  
  # Set GSM7919060 to black
  color_palette[conditions == "GSM7919060"] <- "black"
  
  # Name the color palette
  names(color_palette) <- conditions
  
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Group)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  
  return(p)
}

# Create UMAP plot
umap_plot <- create_umap_plot(merged_seurat, group.by = "condition", title = "UMAP by Condition")

# Save plot
ggsave("imm_UMAP_plot_by_condition.png", plot = umap_plot, width = 10, height = 8, units = "in", dpi = 300)


