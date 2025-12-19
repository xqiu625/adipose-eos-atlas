library(Seurat)
library(ggplot2)
library(dplyr)

# This script is used to label the UCE identity to the seurat object, and then for each cluster, find the UCE identity with the highest percentage
# Load the seurat object
seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(file = seurat_file)

# Add UCE identity to the seurat object
uce_ima <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/predicted_cell_types_metadata_ima.csv")
uce_mouseKidney <- read.csv("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/predicted_cell_types_metadata_mouseKidney.csv")


meta_data <- seurat_object@meta.data
meta_data$cell_barcode <- rownames(meta_data)

# Check the structure of uce_ima and uce_mouseKidney
print(str(uce_ima))
print(str(uce_mouseKidney))

# Add UCE identity to the meta data
meta_data <- meta_data %>%
  left_join(uce_ima, by = "cell_barcode") %>%
  left_join(uce_mouseKidney, by = "cell_barcode")

# Check the structure of meta_data after joins
print(str(meta_data))

# Find the majority cell type for predicted_cell_type, predicted_coarse_cell_type, and predicted_group for each seurat_clusters
uce_identity_by_cluster <- meta_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    majority_predicted_cell_type = names(which.max(table(predicted_cell_type))),
    majority_predicted_coarse_cell_type = names(which.max(table(predicted_coarse_cell_type))),
    majority_predicted_group = names(which.max(table(predicted_group)))
  )

# Print the results
print(uce_identity_by_cluster, n = 46)

# Find the second majority cell type for predicted_cell_type, predicted_coarse_cell_type, and predicted_group for each seurat_clusters
second_uce_identity_by_cluster <- meta_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    second_majority_predicted_cell_type = names(sort(table(predicted_cell_type), decreasing = TRUE)[2]),
    second_majority_predicted_coarse_cell_type = names(sort(table(predicted_coarse_cell_type), decreasing = TRUE)[2]),
    second_majority_predicted_group = names(sort(table(predicted_group), decreasing = TRUE)[2])
  )

# Print the results
print(second_uce_identity_by_cluster, n = 46)

# list the seurat_clusters 28 all cell types and their counts
meta_data %>%
  filter(seurat_clusters == 30) %>%
  select(predicted_cell_type) %>%
  table()

# UMAP plot of the seurat_clusters by predicted_cell_type

seurat_object@meta.data <- meta_data

# Create a stacked bar plot to analyze each cluster's predicted_cell_type composition

# First, ensure we have the necessary libraries
library(ggplot2)
library(dplyr)

# Prepare the data
cluster_composition <- seurat_object@meta.data %>%
  group_by(seurat_clusters, predicted_cell_type) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

# Create a folder to save individual cluster plots
output_folder <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/cluster_composition_plots"
dir.create(output_folder, showWarnings = FALSE)

# Generate a figure for each cluster
for (cluster in unique(cluster_composition$seurat_clusters)) {
  cluster_data <- cluster_composition %>% filter(seurat_clusters == cluster)
  
  cluster_plot <- ggplot(cluster_data, aes(x = reorder(predicted_cell_type, -percentage), y = percentage, fill = predicted_cell_type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Predicted Cell Type", y = "Percentage", fill = "Predicted Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle(paste("Cluster", cluster, "Composition"))
  
  # Save the plot
  output_file <- file.path(output_folder, paste0("Cluster_", cluster, "_Composition.png"))
  ggsave(output_file, plot = cluster_plot, width = 10, height = 6)
  
  # Print a message to confirm the plot has been saved
  print(paste("Cluster", cluster, "plot saved to:", output_file))
}

# Create and save the overall stacked bar plot
stacked_bar_plot <- ggplot(cluster_composition, aes(x = seurat_clusters, y = percentage, fill = predicted_cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Seurat Clusters", y = "Percentage", fill = "Predicted Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Overall Cluster Composition by Predicted Cell Type")

stacked_bar_output_file <- file.path(output_folder, "Overall_Cluster_Composition_Stacked_Bar.png")
ggsave(stacked_bar_output_file, plot = stacked_bar_plot, width = 15, height = 10)

print(paste("Overall stacked bar plot saved to:", stacked_bar_output_file))
