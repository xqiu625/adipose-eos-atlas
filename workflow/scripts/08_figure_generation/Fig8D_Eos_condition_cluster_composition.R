library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
eos <- readRDS(eos_subset_file)
eos@meta.data$Clusters <- paste0("C", eos@meta.data$RNA_snn_res.0.8)

# Set the desired order for conditions
condition_order <- c("WT_Male_HFD", "KO_Male_HFD", "WT_Female_HFD", "KO_Female_HFD", "WT_Female_HFDInfected", "KO_Female_HFDInfected")

# Calculate proportions
cluster_condition_proportions <- eos@meta.data %>%
  group_by(condition, Clusters) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Ensure all conditions are present
all_combinations <- expand.grid(
  condition = factor(condition_order, levels = condition_order),
  Clusters = unique(cluster_condition_proportions$Clusters)
)

cluster_condition_proportions <- merge(
  all_combinations,
  cluster_condition_proportions,
  by = c("condition", "Clusters"),
  all.x = TRUE
)

# Replace NA values with 0
cluster_condition_proportions[is.na(cluster_condition_proportions)] <- 0

# Create stacked barplot with increased font size, Arial font, specified condition order, and no percentage labels
stacked_barplot <- ggplot(cluster_condition_proportions, 
                          aes(x = condition, y = proportion, fill = Clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_discrete(limits = condition_order) +  # Set the order of conditions on x-axis
  labs(title = "Composition of Conditions by Cluster",
       x = "Condition",
       y = "Proportion",
       fill = "Clusters") +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, family = "Arial"),
    axis.title = element_text(size = 16, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, family = "Arial"),
    axis.text.y = element_text(size = 14, family = "Arial"),
    legend.title = element_text(size = 16, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set3")

# Save the plot
ggsave("/bigdata/nairlab/shared/XQiu/EOS_project/figures/grant_figures/Eos_condition_cluster_composition.pdf",
       plot = stacked_barplot,
       width = 12,
       height = 10,
       device = cairo_pdf)
