# CCL Interaction Dotplots - TOP 10 Unique Interactions (v3)
# Purpose: Extract TOP 10 unique interactions and create dotplots with SHARED SCALE
# Updated: 2024-12-13
# Requirements:
#   - Exclude "Other_Macrophages"
#   - TOP 10 unique interactions by strength
#   - Y-axis: 10 unique interactions
#     - Mac-Eos: ordered by Mac subtype (M1, M2 Surface, M2 Functional)
#     - Eos-Mac: ordered by Eos subcluster (C0, C1, C2, C3, C4, C5, C6)
#   - X-axis: 6 conditions (WT Female, WT Male, KO Male, KO Female, KO Female Inf, WT Female Inf)
#   - Color: Viridis (yellow=strongest, purple=weakest)
#   - Size: bigger=stronger, smaller=weaker
#   - SHARED SCALE across both plots for color and size

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# =============================================================================
# SETUP
# =============================================================================

# Input files - FULL raw data (local paths)
input_dir <- '/Users/xinruqiu/Desktop/302025/302025_git/scripts/R20_20251208'
eos_to_mac_file <- file.path(input_dir, "CCL_Eosinophil_to_Macrophage_FULL_raw.csv")
mac_to_eos_file <- file.path(input_dir, "CCL_Macrophage_to_Eosinophil_FULL_raw.csv")

# Output directory
output_dir <- '/Users/xinruqiu/Desktop/302025/302025_git/scripts/R20_20251208/CCL_dotplots_TOP10_v3'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define condition order as specified:
# WT Female, WT Male, KO Male, KO Female, KO Female Inf, WT Female Inf
condition_order <- c(
  "WT_Female_HFD",
  "WT_Male_HFD",
  "KO_Male_HFD",
  "KO_Female_HFD",
  "KO_Female_HFDInfected",
  "WT_Female_HFDInfected"
)

# Clean condition labels for display
condition_labels <- c(
  "WT_Female_HFD" = "WT Female",
  "WT_Male_HFD" = "WT Male",
  "KO_Male_HFD" = "KO Male",
  "KO_Female_HFD" = "KO Female",
  "KO_Female_HFDInfected" = "KO Female Inf",
  "WT_Female_HFDInfected" = "WT Female Inf"
)

# Mac subtype display labels (shorter names)
mac_labels <- c(
  "M1_Macrophages" = "M1",
  "M2_Macrophages_Surface_Only" = "M2 Surface",
  "M2_Macrophages_Functional" = "M2 Functional"
)

# Define Mac subtype order for Y-axis (Mac to Eos plot)
mac_subtype_order <- c("M1_Macrophages", "M2_Macrophages_Surface_Only", "M2_Macrophages_Functional")

# Define Eos subcluster order for Y-axis (Eos to Mac plot)
eos_subcluster_order <- c("Eos_C0", "Eos_C1", "Eos_C2", "Eos_C3", "Eos_C4", "Eos_C5", "Eos_C6")

# Figure dimensions - with legend
fig_height_cm <- 29
fig_width_cm <- 20
fig_height_in <- fig_height_cm / 2.54
fig_width_in <- fig_width_cm / 2.54

# Figure dimensions - no legend (smaller)
fig_height_nolegend_cm <- 7.1
fig_width_nolegend_cm <- 3.71
fig_height_nolegend_in <- fig_height_nolegend_cm / 2.54
fig_width_nolegend_in <- fig_width_nolegend_cm / 2.54

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

eos_to_mac_raw <- read.csv(eos_to_mac_file)
mac_to_eos_raw <- read.csv(mac_to_eos_file)

cat("Eos to Mac raw rows:", nrow(eos_to_mac_raw), "\n")
cat("Mac to Eos raw rows:", nrow(mac_to_eos_raw), "\n")

# =============================================================================
# FILTER OUT "Other_Macrophages"
# =============================================================================

cat("\nFiltering out Other_Macrophages...\n")

# For Eos to Mac: Target is Mac, so filter Target
eos_to_mac_filtered <- eos_to_mac_raw %>%
  filter(Target != "Other_Macrophages")

# For Mac to Eos: Source is Mac, so filter Source
mac_to_eos_filtered <- mac_to_eos_raw %>%
  filter(Source != "Other_Macrophages")

cat("After filtering:\n")
cat("Eos to Mac rows:", nrow(eos_to_mac_filtered), "\n")
cat("Mac to Eos rows:", nrow(mac_to_eos_filtered), "\n")

# =============================================================================
# EXTRACT TOP 20 UNIQUE INTERACTIONS BY STRENGTH
# =============================================================================

get_top_unique_interactions <- function(df, n = 20) {
  # Create unique interaction identifier
  df$Unique_Interaction <- paste(df$Source, df$Target, df$Pathway, sep = ".")

  # Find TOP N unique interactions by their MAX strength across all conditions
  top_interactions <- df %>%
    group_by(Unique_Interaction, Source, Target, Pathway) %>%
    summarise(Max_Strength = max(Strength), .groups = 'drop') %>%
    arrange(desc(Max_Strength)) %>%
    head(n) %>%
    pull(Unique_Interaction)

  cat("Top", n, "unique interactions:\n")
  for (i in 1:length(top_interactions)) {
    cat("  ", i, ". ", top_interactions[i], "\n", sep = "")
  }

  # Filter original data to keep ALL rows for these top N unique interactions
  result <- df %>%
    filter(Unique_Interaction %in% top_interactions)

  return(result)
}

cat("\n=== Extracting TOP 10 unique Eos to Mac interactions ===\n")
eos_to_mac <- get_top_unique_interactions(eos_to_mac_filtered, 10)
cat("\nResult:", nrow(eos_to_mac), "rows,",
    length(unique(eos_to_mac$Unique_Interaction)), "unique interactions\n")

cat("\n=== Extracting TOP 10 unique Mac to Eos interactions ===\n")
mac_to_eos <- get_top_unique_interactions(mac_to_eos_filtered, 10)
cat("\nResult:", nrow(mac_to_eos), "rows,",
    length(unique(mac_to_eos$Unique_Interaction)), "unique interactions\n")

# Export the TOP 10 data
write.csv(eos_to_mac, file.path(output_dir, "CCL_Eosinophil_to_Macrophage_TOP10.csv"), row.names = FALSE)
write.csv(mac_to_eos, file.path(output_dir, "CCL_Macrophage_to_Eosinophil_TOP10.csv"), row.names = FALSE)

# =============================================================================
# PREPARE DATA FOR EOS TO MAC (Y-axis ordered by Eos subcluster: C0, C1, C2...)
# =============================================================================

prepare_eos_to_mac_data <- function(df) {

  # Apply Mac labels to Target
  df$Target_Display <- ifelse(df$Target %in% names(mac_labels),
                               mac_labels[df$Target],
                               df$Target)

  # Create interaction label with short Mac names
  df$Interaction_Label <- paste0(df$Source, " -> ", df$Target_Display, ": ", df$Pathway)

  # Extract Eos subcluster from Source for ordering
  df$Eos_Cluster <- df$Source

  # Create ordering key based on Eos subcluster order (C0, C1, C2, C3, C4, C5, C6)
  df$Eos_Order <- match(df$Eos_Cluster, eos_subcluster_order)
  df$Eos_Order[is.na(df$Eos_Order)] <- 999

  # Get max strength for secondary ordering within cluster
  df <- df %>%
    group_by(Interaction_Label) %>%
    mutate(Max_Strength_Interaction = max(Strength)) %>%
    ungroup()

  # Order interactions: first by Eos cluster, then by strength within cluster (descending)
  interaction_order <- df %>%
    select(Interaction_Label, Eos_Order, Max_Strength_Interaction) %>%
    distinct() %>%
    arrange(Eos_Order, desc(Max_Strength_Interaction)) %>%
    pull(Interaction_Label)

  # Set factor levels (reverse for bottom-to-top display, so C0 appears at top)
  df$Interaction_Label <- factor(df$Interaction_Label, levels = rev(interaction_order))
  df$Condition <- factor(df$Condition, levels = condition_order)

  # Add display labels for conditions
  df$Condition_Display <- condition_labels[as.character(df$Condition)]
  df$Condition_Display <- factor(df$Condition_Display,
                                  levels = condition_labels[condition_order])

  return(df)
}

# =============================================================================
# PREPARE DATA FOR MAC TO EOS (Y-axis ordered by Mac subtype: M1, M2 Surface, M2 Functional)
# =============================================================================

prepare_mac_to_eos_data <- function(df) {

  # Apply Mac labels to Source
  df$Source_Display <- ifelse(df$Source %in% names(mac_labels),
                               mac_labels[df$Source],
                               df$Source)

  # Create interaction label with short Mac names
  df$Interaction_Label <- paste0(df$Source_Display, " -> ", df$Target, ": ", df$Pathway)

  # Extract Mac subtype from Source for ordering
  df$Mac_Subtype <- df$Source

  # Create ordering key based on Mac subtype order (M1, M2 Surface, M2 Functional)
  df$Mac_Order <- match(df$Mac_Subtype, mac_subtype_order)
  df$Mac_Order[is.na(df$Mac_Order)] <- 999

  # Get max strength for secondary ordering within subtype
  df <- df %>%
    group_by(Interaction_Label) %>%
    mutate(Max_Strength_Interaction = max(Strength)) %>%
    ungroup()

  # Order interactions: first by Mac subtype, then by strength within subtype (descending)
  interaction_order <- df %>%
    select(Interaction_Label, Mac_Order, Max_Strength_Interaction) %>%
    distinct() %>%
    arrange(Mac_Order, desc(Max_Strength_Interaction)) %>%
    pull(Interaction_Label)

  # Set factor levels (reverse for bottom-to-top display, so M1 appears at top)
  df$Interaction_Label <- factor(df$Interaction_Label, levels = rev(interaction_order))
  df$Condition <- factor(df$Condition, levels = condition_order)

  # Add display labels for conditions
  df$Condition_Display <- condition_labels[as.character(df$Condition)]
  df$Condition_Display <- factor(df$Condition_Display,
                                  levels = condition_labels[condition_order])

  return(df)
}

eos_to_mac_plot <- prepare_eos_to_mac_data(eos_to_mac)
mac_to_eos_plot <- prepare_mac_to_eos_data(mac_to_eos)

# =============================================================================
# CALCULATE SEPARATE SCALE LIMITS FOR EACH PLOT
# =============================================================================

cat("\nCalculating separate scale limits for each plot...\n")

# Eos to Mac scale
eos_to_mac_min <- min(eos_to_mac_plot$Strength, na.rm = TRUE)
eos_to_mac_max <- max(eos_to_mac_plot$Strength, na.rm = TRUE)
cat("Eos to Mac scale range:", eos_to_mac_min, "to", eos_to_mac_max, "\n")

# Mac to Eos scale
mac_to_eos_min <- min(mac_to_eos_plot$Strength, na.rm = TRUE)
mac_to_eos_max <- max(mac_to_eos_plot$Strength, na.rm = TRUE)
cat("Mac to Eos scale range:", mac_to_eos_min, "to", mac_to_eos_max, "\n")

# =============================================================================
# DOTPLOT 1: Eosinophil to Macrophage (Y-axis by Eos subcluster: C0, C1, C2...)
# =============================================================================

cat("\nCreating Eos to Mac dotplot...\n")

n_interactions_eos <- length(unique(eos_to_mac_plot$Interaction_Label))
cat("Number of unique interactions (Eos to Mac):", n_interactions_eos, "\n")

p1 <- ggplot(eos_to_mac_plot, aes(x = Condition_Display, y = Interaction_Label)) +
  geom_point(aes(size = Strength, color = Strength)) +
  scale_size_continuous(range = c(0.5, 6), limits = c(eos_to_mac_min, eos_to_mac_max), name = "Strength") +
  scale_color_viridis(option = "viridis", direction = 1, limits = c(eos_to_mac_min, eos_to_mac_max), name = "Strength") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8, face = "bold"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "mm")
  ) +
  labs(
    x = "Condition",
    title = "Eos -> Mac (Top 10)"
  )

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_TOP10.pdf"),
  p1,
  width = fig_width_in,
  height = fig_height_in,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_TOP10.png"),
  p1,
  width = fig_width_in,
  height = fig_height_in,
  dpi = 300
)

cat("Saved Eos to Mac dotplot\n")

# =============================================================================
# DOTPLOT 2: Macrophage to Eosinophil (Y-axis by Mac subtype: M1, M2 Surface, M2 Functional)
# =============================================================================

cat("\nCreating Mac to Eos dotplot...\n")

n_interactions_mac <- length(unique(mac_to_eos_plot$Interaction_Label))
cat("Number of unique interactions (Mac to Eos):", n_interactions_mac, "\n")

p2 <- ggplot(mac_to_eos_plot, aes(x = Condition_Display, y = Interaction_Label)) +
  geom_point(aes(size = Strength, color = Strength)) +
  scale_size_continuous(range = c(0.5, 6), limits = c(mac_to_eos_min, mac_to_eos_max), name = "Strength") +
  scale_color_viridis(option = "viridis", direction = 1, limits = c(mac_to_eos_min, mac_to_eos_max), name = "Strength") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8, face = "bold"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "mm")
  ) +
  labs(
    x = "Condition",
    title = "Mac -> Eos (Top 10)"
  )

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_TOP10.pdf"),
  p2,
  width = fig_width_in,
  height = fig_height_in,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_TOP10.png"),
  p2,
  width = fig_width_in,
  height = fig_height_in,
  dpi = 300
)

cat("Saved Mac to Eos dotplot\n")

# =============================================================================
# DOTPLOTS WITHOUT LEGENDS
# =============================================================================

cat("\nCreating dotplots without legends...\n")

# Eos to Mac - no legend, no axis labels
p1_nolegend <- p1 + theme(
  legend.position = "none",
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks = element_blank(),
  plot.title = element_blank()
)

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_TOP10_nolegend.pdf"),
  p1_nolegend,
  width = fig_width_nolegend_in,
  height = fig_height_nolegend_in,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Eos_to_Mac_dotplot_TOP10_nolegend.png"),
  p1_nolegend,
  width = fig_width_nolegend_in,
  height = fig_height_nolegend_in,
  dpi = 300
)

cat("Saved Eos to Mac dotplot (no legend)\n")

# Mac to Eos - no legend, no axis labels
p2_nolegend <- p2 + theme(
  legend.position = "none",
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks = element_blank(),
  plot.title = element_blank()
)

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_TOP10_nolegend.pdf"),
  p2_nolegend,
  width = fig_width_nolegend_in,
  height = fig_height_nolegend_in,
  dpi = 300
)

ggsave(
  file.path(output_dir, "CCL_Mac_to_Eos_dotplot_TOP10_nolegend.png"),
  p2_nolegend,
  width = fig_width_nolegend_in,
  height = fig_height_nolegend_in,
  dpi = 300
)

cat("Saved Mac to Eos dotplot (no legend)\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output directory:", output_dir, "\n\n")
cat("Figure dimensions:", fig_width_cm, "cm (W) x", fig_height_cm, "cm (H)\n\n")
cat("Generated files:\n")
cat("  1. CCL_Eos_to_Mac_dotplot_TOP10.pdf/png\n")
cat("  2. CCL_Mac_to_Eos_dotplot_TOP10.pdf/png\n")
cat("  3. CCL_Eos_to_Mac_dotplot_TOP10_nolegend.pdf/png\n")
cat("  4. CCL_Mac_to_Eos_dotplot_TOP10_nolegend.pdf/png\n")
cat("  5. CCL_Eosinophil_to_Macrophage_TOP10.csv\n")
cat("  6. CCL_Macrophage_to_Eosinophil_TOP10.csv\n")
cat("\nY-axis order:\n")
cat("  Eos->Mac: by Eos subcluster (C0, C1, C2, C3, C4, C5, C6)\n")
cat("  Mac->Eos: by Mac subtype (M1, M2 Surface, M2 Functional)\n")
cat("\nX-axis order:\n")
cat("  ", paste(condition_labels[condition_order], collapse = " -> "), "\n")
cat("\nExcluded: Other_Macrophages\n")
cat("SEPARATE SCALES for each plot:\n")
cat("  Eos to Mac range:", eos_to_mac_min, "to", eos_to_mac_max, "\n")
cat("  Mac to Eos range:", mac_to_eos_min, "to", mac_to_eos_max, "\n")
cat("Color scale: Viridis (yellow=strongest, purple=weakest)\n")
cat("Size: bigger=stronger, smaller=weaker\n")
