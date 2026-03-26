# ==============================================================================
# 00_theme.R
# Shared color palettes and ggplot2 theme for publication-quality figures
# Source this file at the top of each analysis script:
#   source(file.path(dirname(sys.frame(1)$ofile), "00_theme.R"))
#   — or —
#   source("/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/scripts/00_theme.R")
# ==============================================================================

library(ggplot2)

# --- Eosinophil cluster palette (Fig 2A) -------------------------------------
eos_colors <- c(

  "C0" = "#6CA6CD", "C1" = "#E84A5F", "C2" = "#FF6F1F",
  "C3" = "#CBAACB", "C4" = "#B03A48", "C5" = "#5BA4A4", "C6" = "#2F3C7E"
)

# Prefixed version for CellChat (Eos_C0, etc.)
eos_colors_prefixed <- setNames(eos_colors, paste0("Eos_", names(eos_colors)))

# --- Macrophage subtype palette -----------------------------------------------
mac_colors <- c(
  "Mac_M1"      = "#cb416b",
  "Mac_M2_Func" = "#8FC50E",
  "Mac_M2_Surf" = "#5633AA"
)

# --- Combined eos + mac palette -----------------------------------------------
celltype_colors <- c(eos_colors_prefixed, mac_colors)

# --- CellChat group colors (Fig 3A chord) ------------------------------------
group_colors <- c(
  "Eosinophils"    = "#777B7E",
  "M1_Mac"         = "#cb416b",
  "M2_Functional"  = "#8FC50E",
  "M2_Surface"     = "#5633AA"
)

# --- Condition colors (Fig 1B) -----------------------------------------------
cond_colors <- c(
  "WT_Female_HFD"          = "#B2182B",
  "WT_Male_HFD"            = "#2166AC",
  "KO_Female_HFD"          = "#E08080",
  "KO_Male_HFD"            = "#6BAED6",
  "WT_Female_HFDInfected"  = "#D6604D",
  "KO_Female_HFDInfected"  = "#FC8D59"
)

# --- Condition shapes (WT=circle, KO=triangle, Infected=square) ---------------
cond_shapes <- c(
  "WT_Female_HFD" = 16, "WT_Male_HFD" = 16,
  "KO_Female_HFD" = 17, "KO_Male_HFD" = 17,
  "WT_Female_HFDInfected" = 15, "KO_Female_HFDInfected" = 15
)

# --- Condition display labels (with sex symbols) ------------------------------
cond_labels <- c(
  "WT_Female_HFD"          = "WT \u2640",
  "WT_Male_HFD"            = "WT \u2642",
  "KO_Female_HFD"          = "KO \u2640",
  "KO_Male_HFD"            = "KO \u2642",
  "WT_Female_HFDInfected"  = "WT \u2640 Inf",
  "KO_Female_HFDInfected"  = "KO \u2640 Inf"
)

# --- All-cell-type broad palette (Fig 1C, matching manuscript) ----------------
# Names match label-transfer output (plural forms from R02 object)
broad_celltype_colors <- c(
  "Eosinophils"      = "#E07B5F",   # salmon/coral
  "Macrophages"      = "#D4A0C8",   # pink/lavender
  "Dendritic cells"  = "#4D8D8D",   # teal/green
  "B cells"          = "#8CB84B",   # olive/yellow-green
  "NK cells"         = "#3D7D7D",   # darker teal
  "Neutrophils"      = "#2F3C7E",   # dark blue/navy
  "Mast cells"       = "#CBAACB"    # light purple/mauve
)

# --- Volcano plot colors (Fig 1E-H) ------------------------------------------
vol_colors <- c("Up" = "#E84A5F", "Down" = "#2166AC", "NS" = "grey70")

# --- Diverging heatmap palette (Fig 3B difference) ----------------------------
diverging_palette <- colorRampPalette(
  c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
)(100)

# --- Viridis-like palette (Fig 3B baseline) -----------------------------------
viridis_palette <- colorRampPalette(
  c("#440154", "#482878", "#3e4a89", "#31688e",
    "#26838e", "#1f9e89", "#35b779", "#6ece58",
    "#b5de2b", "#fde725")
)(100)

# --- Publication theme (matches manuscript) -----------------------------------
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "") %+replace%
    theme(
      # Title
      plot.title    = element_text(size = rel(1.2), face = "bold", hjust = 0,
                                   margin = margin(b = 8)),
      plot.subtitle = element_text(size = rel(0.9), color = "grey40",
                                   margin = margin(b = 8)),
      # Axes
      axis.title    = element_text(size = rel(1.0), face = "bold"),
      axis.text     = element_text(size = rel(0.85), color = "black"),
      axis.line     = element_line(color = "black", linewidth = 0.5),
      axis.ticks    = element_line(color = "black", linewidth = 0.4),
      # Legend
      legend.title  = element_text(size = rel(0.9), face = "bold"),
      legend.text   = element_text(size = rel(0.8)),
      legend.key.size = unit(0.4, "cm"),
      # Panel
      panel.grid    = element_blank(),
      strip.text    = element_text(size = rel(0.9), face = "bold"),
      strip.background = element_blank(),
      # Margins
      plot.margin   = margin(10, 10, 10, 10)
    )
}

# --- UMAP theme (no axes, clean) ----------------------------------------------
theme_umap <- function(base_size = 12) {
  theme_pub(base_size = base_size) %+replace%
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank()
    )
}

# --- Helper: save publication PDF ---------------------------------------------
save_pdf <- function(plot, path, width, height, ...) {
  cairo_pdf(path, width = width, height = height, ...)
  print(plot)
  dev.off()
  cat("  Saved:", path, "\n")
}

# --- Sample-to-condition mapping (16 samples) ---------------------------------
sample_to_condition <- c(
  "1"  = "WT_Male_HFD_1",       "2"  = "WT_Female_HFD_1",
  "3"  = "KO_Male_HFD_1",       "4"  = "KO_Female_HFD_1",
  "5"  = "WT_Female_HFDInfected_1", "6"  = "WT_Female_HFDInfected_2",
  "7"  = "KO_Female_HFDInfected_1", "8"  = "KO_Female_HFDInfected_2",
  "9"  = "WT_Male_HFD_2",       "10" = "WT_Male_HFD_3",
  "11" = "WT_Female_HFD_2",     "12" = "WT_Female_HFD_3",
  "13" = "KO_Male_HFD_2",       "14" = "KO_Male_HFD_3",
  "15" = "KO_Female_HFD_2",     "16" = "KO_Female_HFD_3"
)

# GEM reaction mapping
gem_reaction_map <- c(
  "1" = "GEM1", "2" = "GEM1", "3" = "GEM1", "4" = "GEM1",
  "5" = "GEM2", "6" = "GEM2", "7" = "GEM2", "8" = "GEM2",
  "9" = "GEM3", "10" = "GEM3", "11" = "GEM3", "12" = "GEM3",
  "13" = "GEM4", "14" = "GEM4", "15" = "GEM4", "16" = "GEM4"
)

cat("  Loaded 00_theme.R: palettes, themes, and helpers\n")
