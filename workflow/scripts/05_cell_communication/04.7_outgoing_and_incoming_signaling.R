library(CellChat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)

# Set working directory
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat")

# Define color scheme for cell types
group_colors <- c(
  "Eosinophils" = "#777B7E",
  "M1_Mac" = "#cb416b",
  "M2_Functional" = "#8FC50E",
  "M2_Surface" = "#5633AA"
)

# Define pathways with biological annotations
# Categorize pathways by their significance
pathway_categories <- list(
  "Macrophage_specific_unchanged" = c("APP"),
  "Eosinophil_specific_unchanged" = c("THBS"),
  "Significantly_changed" = c("FN1", "CCL")
)

all_pathways <- c("APP", "THBS", "FN1", "CCL")

# Function to extract signaling data
extract_signaling_data <- function(cellchat_path, condition_name) {

  cat("\n=== Processing:", condition_name, "===\n")

  cellchat <- readRDS(cellchat_path)

  pathway_names <- cellchat@netP$pathways
  cell_groups <- levels(cellchat@idents)

  cat("Found", length(pathway_names), "pathways\n")

  if (is.null(cellchat@netP$prob)) {
    cat("ERROR: netP$prob is NULL\n")
    return(NULL)
  }

  prob_pathways <- cellchat@netP$prob

  # Initialize result lists
  outgoing_list <- list()
  incoming_list <- list()

  # Calculate for each pathway
  for (pathway in pathway_names) {
    if (pathway %in% dimnames(prob_pathways)[[3]]) {
      pathway_prob <- prob_pathways[,,pathway]

      outgoing <- rowSums(pathway_prob)
      incoming <- colSums(pathway_prob)

      outgoing_list[[pathway]] <- data.frame(
        Pathway = pathway,
        CellType = names(outgoing),
        Strength = as.numeric(outgoing),
        Direction = "Outgoing",
        Condition = condition_name,
        stringsAsFactors = FALSE
      )

      incoming_list[[pathway]] <- data.frame(
        Pathway = pathway,
        CellType = names(incoming),
        Strength = as.numeric(incoming),
        Direction = "Incoming",
        Condition = condition_name,
        stringsAsFactors = FALSE
      )
    }
  }

  combined_data <- rbind(
    do.call(rbind, outgoing_list),
    do.call(rbind, incoming_list)
  )

  cat("✓ Extracted data from", condition_name, "\n")

  return(combined_data)
}

# Process all conditions
cellchat_files <- list(
  "WT Female HFD" = "WT_Female_HFD_cellchat.rds",
  "WT Male HFD" = "WT_Male_HFD_cellchat.rds",
  "WT Female Infected" = "WT_Female_HFDInfected_cellchat.rds",
  "KO Female HFD" = "KO_Female_HFD_cellchat.rds",
  "KO Male HFD" = "KO_Male_HFD_cellchat.rds",
  "KO Female Infected" = "KO_Female_HFDInfected_cellchat.rds"
)

cat("\n========================================\n")
cat("EXTRACTING DATA FROM ALL CONDITIONS\n")
cat("========================================\n")

all_data_list <- list()

for (name in names(cellchat_files)) {
  file_path <- cellchat_files[[name]]
  if (file.exists(file_path)) {
    data <- extract_signaling_data(file_path, name)
    if (!is.null(data)) {
      all_data_list[[name]] <- data
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Combine all data
all_data <- do.call(rbind, all_data_list)

# Filter for our pathways of interest
plot_data <- all_data %>%
  filter(Pathway %in% all_pathways)

# Add pathway category
plot_data$Category <- NA
plot_data$Category[plot_data$Pathway == "APP"] <- "Mac-specific\n(unchanged)"
plot_data$Category[plot_data$Pathway == "THBS"] <- "Eos-specific\n(unchanged)"
plot_data$Category[plot_data$Pathway %in% c("FN1", "CCL")] <- "Changed between\ngroups"

# Reorder factors for better visualization
plot_data$Pathway <- factor(plot_data$Pathway, levels = c("APP", "THBS", "FN1", "CCL"))
plot_data$Category <- factor(plot_data$Category,
                             levels = c("Mac-specific\n(unchanged)",
                                      "Eos-specific\n(unchanged)",
                                      "Changed between\ngroups"))
plot_data$CellType <- factor(plot_data$CellType,
                             levels = c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface"))

# NEW: Reorder Condition factor from least protected to most protected
plot_data$Condition <- factor(plot_data$Condition,
                              levels = c("WT Male HFD",
                                       "KO Male HFD",
                                       "KO Female HFD",
                                       "KO Female Infected",
                                       "WT Female HFD",
                                       "WT Female Infected"))

cat("\n========================================\n")
cat("CREATING IMPROVED VISUALIZATIONS\n")
cat("========================================\n\n")

# Save the data
write.csv(plot_data, "signaling_patterns_all_conditions.csv", row.names = FALSE)
cat("✓ Saved: signaling_patterns_all_conditions.csv\n\n")

# ============================================================================
# Figure 1: Comprehensive heatmap with all conditions
# ============================================================================

cat("Creating comprehensive heatmap...\n")

# Create separate plots for outgoing and incoming
outgoing_data <- plot_data %>% filter(Direction == "Outgoing")
incoming_data <- plot_data %>% filter(Direction == "Incoming")

# Outgoing signaling heatmap
p_outgoing <- ggplot(outgoing_data, aes(x = CellType, y = Condition, fill = Strength)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Strength)), size = 3, fontface = "bold") +
  facet_grid(Category ~ Pathway, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "white",
    mid = "yellow",
    high = "darkred",
    midpoint = median(outgoing_data$Strength, na.rm = TRUE),
    name = "Signaling\nStrength"
  ) +
  labs(
    title = "Outgoing Signaling Patterns",
    subtitle = "Cell types as SENDERS of signals",
    x = "Cell Type (Source)",
    y = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

# Incoming signaling heatmap
p_incoming <- ggplot(incoming_data, aes(x = CellType, y = Condition, fill = Strength)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Strength)), size = 3, fontface = "bold") +
  facet_grid(Category ~ Pathway, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "white",
    mid = "yellow",
    high = "darkred",
    midpoint = median(incoming_data$Strength, na.rm = TRUE),
    name = "Signaling\nStrength"
  ) +
  labs(
    title = "Incoming Signaling Patterns",
    subtitle = "Cell types as RECEIVERS of signals",
    x = "Cell Type (Target)",
    y = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

# Combine plots
pdf("signaling_patterns_comprehensive_heatmap.pdf", width = 20, height = 14)
print(p_outgoing / p_incoming)
dev.off()

cat("✓ Saved: signaling_patterns_comprehensive_heatmap.pdf\n\n")

# ============================================================================
# Figure 2: Bar plots showing differences more clearly
# ============================================================================

cat("Creating bar plot comparison...\n")

# Calculate mean and SD across conditions for each pathway/cell type
summary_data <- plot_data %>%
  group_by(Pathway, CellType, Direction, Category) %>%
  summarise(
    Mean = mean(Strength, na.rm = TRUE),
    SD = sd(Strength, na.rm = TRUE),
    .groups = "drop"
  )

p_bars <- ggplot(summary_data, aes(x = CellType, y = Mean, fill = CellType)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, linewidth = 0.8) +
  facet_grid(Direction ~ Pathway) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Signaling Strength Across Pathways",
    subtitle = "Mean ± SD across all conditions | Pathways categorized by biological significance",
    x = "Cell Type",
    y = "Signaling Strength",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom"
  )

pdf("signaling_patterns_barplot_summary.pdf", width = 14, height = 10)
print(p_bars)
dev.off()

cat("✓ Saved: signaling_patterns_barplot_summary.pdf\n\n")

# ============================================================================
# Figure 3: Focused comparison of changed vs unchanged pathways
# ============================================================================

cat("Creating focused comparison plot...\n")

# Create a simplified version highlighting key differences
focused_data <- plot_data %>%
  mutate(
    PathwayType = ifelse(Pathway %in% c("FN1", "CCL"),
                        "Changed (FN1, CCL)",
                        "Unchanged (APP, THBS)")
  )

p_focused <- ggplot(focused_data, aes(x = Condition, y = Strength, color = CellType, group = CellType)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 2) +
  facet_grid(Direction ~ PathwayType + Pathway, scales = "free") +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Signaling Variability Across Conditions",
    subtitle = "Unchanged pathways (APP, THBS) vs Changed pathways (FN1, CCL) | Ordered from least to most protected",
    x = "Condition",
    y = "Signaling Strength",
    color = "Cell Type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold")
  )

pdf("signaling_patterns_focused_comparison.pdf", width = 18, height = 10)
print(p_focused)
dev.off()

cat("✓ Saved: signaling_patterns_focused_comparison.pdf\n\n")

# ============================================================================
# Summary statistics
# ============================================================================

cat("Calculating summary statistics...\n")

# Calculate coefficient of variation for each pathway across conditions
cv_summary <- plot_data %>%
  group_by(Pathway, CellType, Direction) %>%
  summarise(
    Mean = mean(Strength, na.rm = TRUE),
    SD = sd(Strength, na.rm = TRUE),
    CV = (SD / Mean) * 100,
    Category = first(Category),
    .groups = "drop"
  ) %>%
  arrange(Category, Pathway, Direction, CellType)

write.csv(cv_summary, "signaling_patterns_variability_summary.csv", row.names = FALSE)
cat("✓ Saved: signaling_patterns_variability_summary.csv\n\n")

# Print key findings
cat("========================================\n")
cat("KEY FINDINGS SUMMARY\n")
cat("========================================\n\n")

cat("UNCHANGED PATHWAYS (low variability):\n")
unchanged_cv <- cv_summary %>% filter(Pathway %in% c("APP", "THBS"))
print(unchanged_cv %>% select(Pathway, CellType, Direction, Mean, CV))

cat("\n\nCHANGED PATHWAYS (higher variability):\n")
changed_cv <- cv_summary %>% filter(Pathway %in% c("FN1", "CCL"))
print(changed_cv %>% select(Pathway, CellType, Direction, Mean, CV))

cat("\n========================================\n")
cat("✓✓✓ ALL VISUALIZATIONS COMPLETE ✓✓✓\n")
cat("========================================\n\n")

cat("Generated files:\n")
cat("  1. signaling_patterns_comprehensive_heatmap.pdf - Full heatmap with all conditions\n")
cat("  2. signaling_patterns_barplot_summary.pdf - Bar plots with error bars\n")
cat("  3. signaling_patterns_focused_comparison.pdf - Line plots showing variability\n")
cat("  4. signaling_patterns_all_conditions.csv - Raw data\n")
cat("  5. signaling_patterns_variability_summary.csv - Statistical summary\n\n")

cat("INTERPRETATION GUIDE:\n")
cat("  • APP (Macrophage-specific): Consistent across conditions\n")
cat("  • THBS (Eosinophil-specific): Consistent across conditions\n")
cat("  • FN1 & CCL: Show significant variation between groups\n")
cat("  • Condition order: From least protected to most protected\n")
cat("========================================\n")
