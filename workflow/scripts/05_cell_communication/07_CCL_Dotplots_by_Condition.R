# CCL Interaction Dotplots by Condition
# Purpose: Create dotplots showing CCL interactions across 6 experimental conditions
# X-axis: Conditions in specified order
# Date: 2024-12-08

library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# SETUP
# =============================================================================

# Input files (local paths)
eos_to_mac_file <- "/Users/xinruqiu/Desktop/302025/302025_git/scripts/R20_20251208/CCL_Eosinophil_to_Macrophage_TOP20.csv"
mac_to_eos_file <- "/Users/xinruqiu/Desktop/302025/302025_git/scripts/R20_20251208/CCL_Macrophage_to_Eosinophil_TOP20.csv"

# Output directory
output_dir <- "/Users/xinruqiu/Desktop/302025/302025_git/scripts/R20_20251208/CCL_dotplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define condition order as specified
condition_order <- c(
  "WT_Male_HFD",
  "KO_Male_HFD",
  "KO_Female_HFD",
  "KO_Female_HFDInfected",
  "WT_Female_HFD",
  "WT_Female_HFDInfected"
)

# Clean condition labels for display
condition_labels <- c(
  "WT_Male_HFD" = "WT Male HFD",
  "KO_Male_HFD" = "KO Male HFD",
  "KO_Female_HFD" = "KO Female HFD",
  "KO_Female_HFDInfected" = "KO Female Infected",
  "WT_Female_HFD" = "WT Female HFD",
  "WT_Female_HFDInfected" = "WT Female Infected"
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

eos_to_mac <- read.csv(eos_to_mac_file)
mac_to_eos <- read.csv(mac_to_eos_file)

cat("Eos to Mac rows:", nrow(eos_to_mac), "\n")
cat("Mac to Eos rows:", nrow(mac_to_eos), "\n")

# Check conditions present
cat("\nConditions in Eos to Mac:\n")
print(unique(eos_to_mac$Condition))
cat("\nConditions in Mac to Eos:\n")
print(unique(mac_to_eos$Condition))

# =============================================================================
# PREPARE DATA
# =============================================================================

# Function to prepare data for dotplot
prepare_dotplot_data <- function(df, direction_label) {

  # Create interaction label (Source -> Target : Pathway)
  df$Interaction_Label <- paste0(df$Source, " -> ", df$Target, "\n(", df$Pathway, ")")

  # Get unique interactions and their max strength for ordering
  interaction_order <- df %>%
    group_by(Interaction_Label) %>%
    summarise(max_strength = max(Strength), .groups = 'drop') %>%
    arrange(desc(max_strength)) %>%
    pull(Interaction_Label)

  # Set factor levels
  df$Interaction_Label <- factor(df$Interaction_Label, levels = rev(interaction_order))
  df$Condition <- factor(df$Condition, levels = condition_order)

  # Add display labels for conditions
  df$Condition_Display <- condition_labels[as.character(df$Condition)]
  df$Condition_Display <- factor(df$Condition_Display,
                                  levels = condition_labels[condition_order])

  return(df)
}

eos_to_mac_plot <- prepare_dotplot_data(eos_to_mac, "Eos to Mac")
mac_to_eos_plot <- prepare_dotplot_data(mac_to_eos, "Mac to Eos")

# =============================================================================
# DOTPLOT 1: Eosinophil to Macrophage
# =============================================================================

cat("\nCreating Eos to Mac dotplot...\n")

# Get unique interactions for this dataset
n_interactions_eos <- length(unique(eos_to_mac_plot$Interaction_Label))

p1 <- ggplot(eos_to_mac_plot, aes(x = Condition_Display, y = Interaction_Label)) +
  geom_point(aes(size = Strength, color = Strength)) +
  scale_size_continuous(range = c(1, 8), name = "Strength") +
  scale_color_gradient(low = "lightblue", high = "darkred", name = "Strength") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Interaction (Source -> Target)",
    title = "CCL Signaling: Eosinophil to Macrophage\n(Top 20 Interactions)"
  )

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot.pdf"),
  p1,
  width = 10,
  height = max(6, n_interactions_eos * 0.4),
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot.png"),
  p1,
  width = 10,
  height = max(6, n_interactions_eos * 0.4),
  dpi = 300
)

cat("Saved Eos to Mac dotplot\n")

# =============================================================================
# DOTPLOT 2: Macrophage to Eosinophil
# =============================================================================

cat("\nCreating Mac to Eos dotplot...\n")

n_interactions_mac <- length(unique(mac_to_eos_plot$Interaction_Label))

p2 <- ggplot(mac_to_eos_plot, aes(x = Condition_Display, y = Interaction_Label)) +
  geom_point(aes(size = Strength, color = Strength)) +
  scale_size_continuous(range = c(1, 8), name = "Strength") +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "Strength") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Interaction (Source -> Target)",
    title = "CCL Signaling: Macrophage to Eosinophil\n(Top 20 Interactions)"
  )

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot.pdf"),
  p2,
  width = 10,
  height = max(6, n_interactions_mac * 0.4),
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot.png"),
  p2,
  width = 10,
  height = max(6, n_interactions_mac * 0.4),
  dpi = 300
)

cat("Saved Mac to Eos dotplot\n")

# =============================================================================
# ALTERNATIVE: AGGREGATED BY PATHWAY (Mean across cell type pairs)
# =============================================================================

cat("\nCreating aggregated dotplots by pathway...\n")

# Aggregate Eos to Mac by Pathway and Condition
eos_to_mac_agg <- eos_to_mac %>%
  group_by(Pathway, Condition) %>%
  summarise(
    Mean_Strength = mean(Strength),
    Max_Strength = max(Strength),
    N_Interactions = n(),
    .groups = 'drop'
  )

eos_to_mac_agg$Condition <- factor(eos_to_mac_agg$Condition, levels = condition_order)
eos_to_mac_agg$Condition_Display <- condition_labels[as.character(eos_to_mac_agg$Condition)]
eos_to_mac_agg$Condition_Display <- factor(eos_to_mac_agg$Condition_Display,
                                            levels = condition_labels[condition_order])

# Order pathways by mean strength
pathway_order_eos <- eos_to_mac_agg %>%
  group_by(Pathway) %>%
  summarise(total = sum(Mean_Strength)) %>%
  arrange(desc(total)) %>%
  pull(Pathway)

eos_to_mac_agg$Pathway <- factor(eos_to_mac_agg$Pathway, levels = rev(pathway_order_eos))

p3 <- ggplot(eos_to_mac_agg, aes(x = Condition_Display, y = Pathway)) +
  geom_point(aes(size = Mean_Strength, color = Mean_Strength)) +
  scale_size_continuous(range = c(2, 10), name = "Mean\nStrength") +
  scale_color_gradient(low = "lightyellow", high = "darkred", name = "Mean\nStrength") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "CCL Pathway",
    title = "CCL Signaling: Eosinophil to Macrophage\n(Aggregated by Pathway)"
  )

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_aggregated.pdf"),
  p3,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_aggregated.png"),
  p3,
  width = 9,
  height = 5,
  dpi = 300
)

# Aggregate Mac to Eos by Pathway and Condition
mac_to_eos_agg <- mac_to_eos %>%
  group_by(Pathway, Condition) %>%
  summarise(
    Mean_Strength = mean(Strength),
    Max_Strength = max(Strength),
    N_Interactions = n(),
    .groups = 'drop'
  )

mac_to_eos_agg$Condition <- factor(mac_to_eos_agg$Condition, levels = condition_order)
mac_to_eos_agg$Condition_Display <- condition_labels[as.character(mac_to_eos_agg$Condition)]
mac_to_eos_agg$Condition_Display <- factor(mac_to_eos_agg$Condition_Display,
                                            levels = condition_labels[condition_order])

pathway_order_mac <- mac_to_eos_agg %>%
  group_by(Pathway) %>%
  summarise(total = sum(Mean_Strength)) %>%
  arrange(desc(total)) %>%
  pull(Pathway)

mac_to_eos_agg$Pathway <- factor(mac_to_eos_agg$Pathway, levels = rev(pathway_order_mac))

p4 <- ggplot(mac_to_eos_agg, aes(x = Condition_Display, y = Pathway)) +
  geom_point(aes(size = Mean_Strength, color = Mean_Strength)) +
  scale_size_continuous(range = c(2, 10), name = "Mean\nStrength") +
  scale_color_gradient(low = "lightcyan", high = "darkblue", name = "Mean\nStrength") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "CCL Pathway",
    title = "CCL Signaling: Macrophage to Eosinophil\n(Aggregated by Pathway)"
  )

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_aggregated.pdf"),
  p4,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_aggregated.png"),
  p4,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Saved aggregated dotplots\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output directory:", output_dir, "\n\n")
cat("Generated files:\n")
cat("  1. CCL_Eos_to_Mac_dotplot.pdf/png - Individual interactions\n")
cat("  2. CCL_Mac_to_Eos_dotplot.pdf/png - Individual interactions\n")
cat("  3. CCL_Eos_to_Mac_dotplot_aggregated.pdf/png - Aggregated by pathway\n")
cat("  4. CCL_Mac_to_Eos_dotplot_aggregated.pdf/png - Aggregated by pathway\n")
cat("\nX-axis order:\n")
cat("  ", paste(condition_labels[condition_order], collapse = ", "), "\n")
