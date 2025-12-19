library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(tibble)

# Load the seurat object
seurat_object <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/seurat_object_for_cellchat.rds")

out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'
dir.create(out_dir, recursive = TRUE)

# ===============================================
# PSEUDOBULK FUNCTION FOR PCA (averaged log-normalized expression)
# ===============================================

create_pseudobulk_for_pca <- function(seurat_obj, group_by = "orig.ident", assay = "RNA", slot = "data") {
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  meta <- seurat_obj@meta.data
  meta$group <- meta[[group_by]]
  
  pseudobulk_list <- list()
  for (group in unique(meta$group)) {
    cells <- rownames(meta)[meta$group == group]
    if (length(cells) > 1) {
      pseudobulk_list[[group]] <- Matrix::rowMeans(expr[, cells])
    } else {
      pseudobulk_list[[group]] <- expr[, cells]
    }
  }
  
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  return(pseudobulk_matrix)
}

# Helper function to create pseudobulk and metadata
create_pseudobulk_and_metadata <- function(seurat_obj, output_prefix) {
  cat("Creating pseudobulk matrix for", output_prefix, "...\n")
  
  # Create pseudobulk matrix
  pseudobulk_matrix <- create_pseudobulk_for_pca(seurat_obj, group_by = "orig.ident")
  
  # Filter genes with zero variance
  gene_vars <- apply(pseudobulk_matrix, 1, var)
  pseudobulk_filtered <- pseudobulk_matrix[gene_vars > 0, ]
  cat(output_prefix, ": Removed", sum(gene_vars == 0), "zero-variance genes\n")
  
  # Create sample metadata
  sample_metadata <- seurat_obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      genotype = first(genotype),
      Sex = first(Sex),
      condition = first(condition),
      infection_status = first(infection_status),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    column_to_rownames("orig.ident")
  
  # Align metadata with pseudobulk columns
  sample_metadata <- sample_metadata[colnames(pseudobulk_filtered), ]
  
  # Save files
  saveRDS(pseudobulk_filtered, file.path(out_dir, paste0("pseudobulk_", output_prefix, ".rds")))
  saveRDS(sample_metadata, file.path(out_dir, paste0("sample_metadata_", output_prefix, ".rds")))
  
  cat(output_prefix, "pseudobulk matrix dimensions:", dim(pseudobulk_filtered), "\n")
  cat(output_prefix, "samples:", ncol(pseudobulk_filtered), "\n\n")
  
  return(list(pseudobulk = pseudobulk_filtered, metadata = sample_metadata))
}

# ===============================================
# 1. ALL CELLS - THREE CONDITIONS
# ===============================================

# 1a. All cells - All samples (16 samples)
all_cells_all <- create_pseudobulk_and_metadata(seurat_object, "all_cells_all_samples")

# 1b. All cells - Not infected samples only
seurat_not_infected <- subset(seurat_object, subset = infection_status == "Not_Infected")
all_cells_not_infected <- create_pseudobulk_and_metadata(seurat_not_infected, "all_cells_not_infected")

# 1c. All cells - Infected samples only
seurat_infected <- subset(seurat_object, subset = infection_status == "Infected")
all_cells_infected <- create_pseudobulk_and_metadata(seurat_infected, "all_cells_infected")

# ===============================================
# 2. EOSINOPHILS - THREE CONDITIONS
# ===============================================

# Subset to eosinophils
eos_seurat <- subset(seurat_object, subset = cell_type == "Eosinophils")

# 2a. Eosinophils - All samples
eos_all <- create_pseudobulk_and_metadata(eos_seurat, "eosinophils_all_samples")

# 2b. Eosinophils - Not infected samples only
eos_seurat_not_infected <- subset(eos_seurat, subset = infection_status == "Not_Infected")
eos_not_infected <- create_pseudobulk_and_metadata(eos_seurat_not_infected, "eosinophils_not_infected")

# 2c. Eosinophils - Infected samples only
eos_seurat_infected <- subset(eos_seurat, subset = infection_status == "Infected")
eos_infected <- create_pseudobulk_and_metadata(eos_seurat_infected, "eosinophils_infected")

# ===============================================
# 3. MACROPHAGES - THREE CONDITIONS
# ===============================================

# Subset to macrophages
mac_seurat <- subset(seurat_object, subset = cell_type == "Macrophages")

# 3a. Macrophages - All samples
mac_all <- create_pseudobulk_and_metadata(mac_seurat, "macrophages_all_samples")

# 3b. Macrophages - Not infected samples only
mac_seurat_not_infected <- subset(mac_seurat, subset = infection_status == "Not_Infected")
mac_not_infected <- create_pseudobulk_and_metadata(mac_seurat_not_infected, "macrophages_not_infected")

# 3c. Macrophages - Infected samples only
mac_seurat_infected <- subset(mac_seurat, subset = infection_status == "Infected")
mac_infected <- create_pseudobulk_and_metadata(mac_seurat_infected, "macrophages_infected")

# ===============================================
# 4. COMBINED EOSINOPHILS AND MACROPHAGES - THREE CONDITIONS
# ===============================================

# Subset to both eosinophils and macrophages
eos_mac_seurat <- subset(seurat_object, subset = cell_type %in% c("Eosinophils", "Macrophages"))

# 4a. Eosinophils + Macrophages - All samples
eos_mac_all <- create_pseudobulk_and_metadata(eos_mac_seurat, "eos_mac_all_samples")

# 4b. Eosinophils + Macrophages - Not infected samples only
eos_mac_seurat_not_infected <- subset(eos_mac_seurat, subset = infection_status == "Not_Infected")
eos_mac_not_infected <- create_pseudobulk_and_metadata(eos_mac_seurat_not_infected, "eos_mac_not_infected")

# 4c. Eosinophils + Macrophages - Infected samples only
eos_mac_seurat_infected <- subset(eos_mac_seurat, subset = infection_status == "Infected")
eos_mac_infected <- create_pseudobulk_and_metadata(eos_mac_seurat_infected, "eos_mac_infected")

# ===============================================
# SUMMARY
# ===============================================

cat("\n=== COMPREHENSIVE PSEUDOBULK SUMMARY ===\n")
cat("Output directory:", out_dir, "\n\n")

cat("1. ALL CELLS:\n")
cat("   - All samples (16):", dim(all_cells_all$pseudobulk), "\n")
cat("   - Not infected:", dim(all_cells_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(all_cells_infected$pseudobulk), "\n\n")

cat("2. EOSINOPHILS:\n")
cat("   - All samples:", dim(eos_all$pseudobulk), "\n")
cat("   - Not infected:", dim(eos_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(eos_infected$pseudobulk), "\n\n")

cat("3. MACROPHAGES:\n")
cat("   - All samples:", dim(mac_all$pseudobulk), "\n")
cat("   - Not infected:", dim(mac_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(mac_infected$pseudobulk), "\n\n")

cat("4. EOSINOPHILS + MACROPHAGES COMBINED:\n")
cat("   - All samples:", dim(eos_mac_all$pseudobulk), "\n")
cat("   - Not infected:", dim(eos_mac_not_infected$pseudobulk), "\n")
cat("   - Infected:", dim(eos_mac_infected$pseudobulk), "\n\n")

cat("FILES CREATED:\n")
files_created <- c(
  # All cells
  "pseudobulk_all_cells_all_samples.rds",
  "pseudobulk_all_cells_not_infected.rds", 
  "pseudobulk_all_cells_infected.rds",
  "sample_metadata_all_cells_all_samples.rds",
  "sample_metadata_all_cells_not_infected.rds",
  "sample_metadata_all_cells_infected.rds",
  
  # Eosinophils only
  "pseudobulk_eosinophils_all_samples.rds",
  "pseudobulk_eosinophils_not_infected.rds",
  "pseudobulk_eosinophils_infected.rds",
  "sample_metadata_eosinophils_all_samples.rds",
  "sample_metadata_eosinophils_not_infected.rds",
  "sample_metadata_eosinophils_infected.rds",
  
  # Macrophages only
  "pseudobulk_macrophages_all_samples.rds",
  "pseudobulk_macrophages_not_infected.rds",
  "pseudobulk_macrophages_infected.rds",
  "sample_metadata_macrophages_all_samples.rds",
  "sample_metadata_macrophages_not_infected.rds",
  "sample_metadata_macrophages_infected.rds",
  
  # Eosinophils + Macrophages combined
  "pseudobulk_eos_mac_all_samples.rds",
  "pseudobulk_eos_mac_not_infected.rds",
  "pseudobulk_eos_mac_infected.rds",
  "sample_metadata_eos_mac_all_samples.rds",
  "sample_metadata_eos_mac_not_infected.rds",
  "sample_metadata_eos_mac_infected.rds"
)

for (file in files_created) {
  cat("   -", file, "\n")
}

cat("\nTotal files created:", length(files_created), "\n")

# ===============================================
# CELL TYPE COMPOSITION SUMMARY
# ===============================================

cat("\n=== CELL TYPE COMPOSITION SUMMARY ===\n")

# Check cell counts in combined datasets
cat("Cell counts in combined eosinophils + macrophages datasets:\n")

cell_type_summary <- eos_mac_seurat@meta.data %>%
  group_by(orig.ident, cell_type, infection_status) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  pivot_wider(names_from = cell_type, values_from = n_cells, values_fill = 0)

print(cell_type_summary)

cat("\nTotal cells per sample in combined dataset:\n")
total_cells_summary <- eos_mac_seurat@meta.data %>%
  group_by(orig.ident, infection_status) %>%
  summarise(
    total_cells = n(),
    eosinophils = sum(cell_type == "Eosinophils"),
    macrophages = sum(cell_type == "Macrophages"),
    eos_percentage = round(100 * sum(cell_type == "Eosinophils") / n(), 1),
    mac_percentage = round(100 * sum(cell_type == "Macrophages") / n(), 1),
    .groups = 'drop'
  )

print(total_cells_summary)


# PCA and Clustering Analysis for Combined Eosinophils and Macrophages
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(dplyr)
library(tibble)
library(dbscan)
library(mclust)

# Set directories
out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'
pca_plots_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/pca_plots_eos_mac'
dir.create(pca_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# FUNCTION TO PERFORM PCA ANALYSIS
# ===============================================

perform_pca_analysis <- function(pseudobulk_file, metadata_file, output_prefix) {
  
  cat("\n===============================================\n")
  cat("Processing PCA for:", output_prefix, "\n")
  cat("===============================================\n")
  
  # Load data
  print("Loading pseudobulk data...")
  pseudobulk_all <- readRDS(file.path(out_dir, pseudobulk_file))
  sample_metadata <- readRDS(file.path(out_dir, metadata_file))
  
  print(paste("Loaded data:", nrow(pseudobulk_all), "genes x", ncol(pseudobulk_all), "samples"))
  print("Sample metadata:")
  print(sample_metadata)
  
  # Check if we have enough samples for analysis
  if (ncol(pseudobulk_all) < 3) {
    cat("Warning: Too few samples (", ncol(pseudobulk_all), ") for meaningful PCA analysis. Skipping...\n")
    return(NULL)
  }
  
  # ===============================================
  # PREPARE DATA FOR PCA
  # ===============================================
  
  # Transpose for PCA (samples as rows, genes as columns)
  pca_data <- t(pseudobulk_all)
  
  # Log2 transform with pseudocount
  pca_data_log <- log2(pca_data + 1)
  
  # Filter genes - keep most variable genes for PCA
  gene_vars <- apply(pca_data_log, 2, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(2000, length(gene_vars))]
  pca_data_filtered <- pca_data_log[, top_genes]
  
  print(paste("Using top", ncol(pca_data_filtered), "most variable genes for PCA"))
  
  # ===============================================
  # PERFORM PCA USING PRCOMP
  # ===============================================
  
  print("Performing PCA...")
  pca_result <- prcomp(pca_data_filtered, scale. = TRUE, center = TRUE)
  
  # Get variance explained
  n_pcs <- min(10, ncol(pca_result$x))
  variance_explained <- summary(pca_result)$importance[2, 1:n_pcs] * 100
  
  print("Variance explained by first PCs:")
  print(round(variance_explained, 2))
  
  # ===============================================
  # PREPARE PLOT DATA
  # ===============================================
  
  # Extract PC scores
  pc_scores <- data.frame(pca_result$x[, 1:n_pcs])
  pc_scores$sample <- rownames(pc_scores)
  
  # Merge with metadata
  plot_data <- pc_scores %>%
    rownames_to_column("sample_id") %>%
    left_join(sample_metadata %>% rownames_to_column("sample_id"), by = "sample_id")
  
  # Create plotting variables
  plot_data$infection_fill <- ifelse(plot_data$infection_status == "Infected", "filled", "clear")
  plot_data$sex_color <- ifelse(plot_data$Sex == "Male", "#6CA6CD", "#D82020")
  plot_data$genotype_shape <- ifelse(plot_data$genotype == "WT", "square", "triangle")
  
  # Convert to factors for ggplot
  plot_data$infection_status <- factor(plot_data$infection_status, levels = c("Not_Infected", "Infected"))
  plot_data$Sex <- factor(plot_data$Sex)
  plot_data$genotype <- factor(plot_data$genotype, levels = c("WT", "KO"))
  
  print("Plot data prepared:")
  print(table(plot_data$infection_status))
  print(table(plot_data$Sex))
  print(table(plot_data$genotype))
  
  # ===============================================
  # HIERARCHICAL CLUSTERING WITH DIFFERENT K VALUES
  # ===============================================
  
  # Try multiple cluster numbers (2 to min(5, n_samples-1))
  max_k <- min(5, nrow(plot_data) - 1)
  for (k in 2:max_k) {
    tryCatch({
      # Perform hierarchical clustering
      dist_matrix <- dist(plot_data[, c("PC1", "PC2")])
      hc <- hclust(dist_matrix, method = "ward.D2")
      cluster_assignments <- cutree(hc, k = k)
      
      # Add Cluster column
      plot_data$Cluster <- as.factor(cluster_assignments)
      
      # Generate plot
      pca_plot_k <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
        stat_ellipse(aes(group = Cluster, color = Cluster),
                     type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
        geom_point(aes(
          color = Sex,
          shape = genotype,
          fill = infection_status
        ), size = 4, stroke = 1.5) +
        scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
        scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
        scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
        guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "grey")),
          color = guide_legend(override.aes = list(shape = 16))
        ) +
        labs(
          x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
          y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
          title = paste("PCA:", output_prefix, "(k =", k, "clusters)"),
          subtitle = "Dashed ellipses represent similarity-based clusters",
          color = "Sex / Cluster",
          shape = "Genotype", 
          fill = "Infection Status"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right"
        )
      
      # Save the plot
      ggsave(
        filename = file.path(pca_plots_dir, paste0(output_prefix, "_pca_k", k, "_clusters.pdf")),
        plot = pca_plot_k,
        width = 10, height = 8
      )
      
      cat(paste("Saved PCA plot with k =", k, "clusters\n"))
    }, error = function(e) {
      cat("Error with k =", k, ":", e$message, "\n")
    })
  }
  
  # ===============================================
  # DIFFERENT HIERARCHICAL CLUSTERING METHODS
  # ===============================================
  
  methods <- c("ward.D2", "average", "complete", "single")
  
  for (m in methods) {
    tryCatch({
      hc <- hclust(dist(plot_data[, c("PC1", "PC2")]), method = m)
      cluster_assignments <- cutree(hc, k = 3)
      plot_data$Cluster <- as.factor(cluster_assignments)
      
      # Create plot
      p <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
        stat_ellipse(aes(group = Cluster, color = Cluster),
                     type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
        geom_point(aes(color = Sex, shape = genotype, fill = infection_status),
                   size = 4, stroke = 1.5) +
        scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
        scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
        scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
        guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "grey")),
          color = guide_legend(override.aes = list(shape = 16))
        ) +
        labs(
          x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
          y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
          title = paste("PCA:", output_prefix, "with", m, "clustering"),
          subtitle = "Dashed ellipses show similarity-based clusters",
          color = "Sex / Cluster",
          shape = "Genotype", 
          fill = "Infection Status"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right"
        )
      
      ggsave(file.path(pca_plots_dir, paste0(output_prefix, "_pca_cluster_", m, ".pdf")), 
             p, width = 10, height = 8)
      cat("Saved:", m, "clustering\n")
    }, error = function(e) {
      cat("Error with", m, "method:", e$message, "\n")
    })
  }
  
  # ===============================================
  # KMEANS CLUSTERING
  # ===============================================
  
  tryCatch({
    # Run k-means clustering
    set.seed(123)
    k_centers <- min(3, nrow(plot_data) - 1)
    kmeans_result <- kmeans(plot_data[, c("PC1", "PC2")], centers = k_centers, nstart = 25)
    plot_data$Cluster <- as.factor(kmeans_result$cluster)
    
    # Create PCA plot with cluster ellipses and k-means coloring
    pca_plot_kmeans <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
      stat_ellipse(aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
        title = paste("PCA:", output_prefix, "with K-means Clustering"),
        subtitle = paste("Dashed ellipses = k-means clusters (k=", k_centers, ")"),
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    # Save the plot
    ggsave(file.path(pca_plots_dir, paste0(output_prefix, "_pca_kmeans_clusters.pdf")),
           plot = pca_plot_kmeans, width = 10, height = 8)
    
    cat("Saved PCA plot with k-means clustering\n")
  }, error = function(e) {
    cat("Error with K-means:", e$message, "\n")
  })
  
  # ===============================================
  # DBSCAN CLUSTERING
  # ===============================================
  
  tryCatch({
    # Run DBSCAN clustering
    db <- dbscan::dbscan(plot_data[, c("PC1", "PC2")], eps = 3, minPts = 3)
    plot_data$Cluster <- as.factor(db$cluster)
    
    # Optional: Label noise points (cluster 0) as "Noise"
    plot_data$Cluster <- factor(plot_data$Cluster, 
                                levels = sort(unique(db$cluster)), 
                                labels = ifelse(sort(unique(db$cluster)) == 0, "Noise", paste0("Cluster_", sort(unique(db$cluster)))))
    
    # PCA plot with DBSCAN clusters
    pca_plot_dbscan <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
      stat_ellipse(data = subset(plot_data, Cluster != "Noise"), 
                   aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
        title = paste("PCA:", output_prefix, "with DBSCAN Clustering"),
        subtitle = paste("DBSCAN: eps = 3, minPts = 3 | 'Noise' = unclustered samples"),
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    # Save the PDF
    ggsave(file.path(pca_plots_dir, paste0(output_prefix, "_pca_dbscan_clusters.pdf")),
           plot = pca_plot_dbscan, width = 10, height = 8)
    
    cat("Saved PCA plot with DBSCAN clustering\n")
  }, error = function(e) {
    cat("Error with DBSCAN:", e$message, "\n")
  })
  
  # ===============================================
  # GMM CLUSTERING
  # ===============================================
  
  tryCatch({
    # Load library and run GMM clustering with forced fewer clusters
    max_clusters <- min(4, nrow(plot_data) - 1)
    gmm <- Mclust(plot_data[, c("PC1", "PC2")], G = 2:max_clusters)
    plot_data$Cluster <- as.factor(gmm$classification)
    
    # Optional: relabel clusters for clarity
    levels(plot_data$Cluster) <- paste0("Cluster_", levels(plot_data$Cluster))
    
    # Create PCA plot
    pca_plot_gmm <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
      stat_ellipse(aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
        title = paste("PCA:", output_prefix, "with GMM Clustering (mclust)"),
        subtitle = paste("Model:", gmm$modelName, "| Selected clusters:", gmm$G),
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    # Save the plot
    ggsave(file.path(pca_plots_dir, paste0(output_prefix, "_pca_gmm_clusters.pdf")),
           plot = pca_plot_gmm, width = 10, height = 8)
    
    cat("Saved PCA plot with GMM clustering\n")
  }, error = function(e) {
    cat("Error with GMM:", e$message, "\n")
  })
  
  cat("Completed PCA analysis for:", output_prefix, "\n")
}

# ===============================================
# DEFINE COMBINED EOS+MAC DATASETS TO ANALYZE
# ===============================================

datasets <- list(
  # Combined Eosinophils + Macrophages
  list(pseudobulk = "pseudobulk_eos_mac_all_samples.rds", 
       metadata = "sample_metadata_eos_mac_all_samples.rds", 
       prefix = "eos_mac_all_samples"),
  list(pseudobulk = "pseudobulk_eos_mac_not_infected.rds", 
       metadata = "sample_metadata_eos_mac_not_infected.rds", 
       prefix = "eos_mac_not_infected"),
  list(pseudobulk = "pseudobulk_eos_mac_infected.rds", 
       metadata = "sample_metadata_eos_mac_infected.rds", 
       prefix = "eos_mac_infected")
)

# ===============================================
# RUN ANALYSIS FOR COMBINED DATASETS
# ===============================================

cat("Starting PCA analysis for combined eosinophils and macrophages datasets...\n")
cat("Output directory:", pca_plots_dir, "\n")

for (dataset in datasets) {
  perform_pca_analysis(dataset$pseudobulk, dataset$metadata, dataset$prefix)
}

cat("\n===============================================\n")
cat("PCA ANALYSIS COMPLETE!\n")
cat("===============================================\n")
cat("All PCA plots saved to:", pca_plots_dir, "\n")
cat("Total datasets analyzed:", length(datasets), "\n")
cat("Each dataset includes multiple clustering methods:\n")
cat("- Hierarchical clustering (k=2 to 5)\n")
cat("- Different linkage methods (ward.D2, average, complete, single)\n")
cat("- K-means clustering\n")
cat("- DBSCAN clustering\n")
cat("- GMM clustering (mclust) with forced fewer clusters\n")


# MDS and Clustering Analysis for Combined Eosinophils and Macrophages
library(ggplot2)
library(dplyr)
library(tibble)
library(dbscan)
library(mclust)

# Set directories
out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'
mds_plots_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/mds_plots_eos_mac'
dir.create(mds_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# FUNCTION TO PERFORM COMPREHENSIVE MDS ANALYSIS
# ===============================================

perform_mds_analysis <- function(pseudobulk_file, metadata_file, output_prefix) {
  
  cat("\n===============================================\n")
  cat("Processing MDS for:", output_prefix, "\n")
  cat("===============================================\n")
  
  # Load data
  pseudobulk_data <- readRDS(file.path(out_dir, pseudobulk_file))
  sample_metadata <- readRDS(file.path(out_dir, metadata_file))
  
  cat("Loaded data:", nrow(pseudobulk_data), "genes x", ncol(pseudobulk_data), "samples\n")
  
  # Check if we have enough samples for analysis
  if (ncol(pseudobulk_data) < 3) {
    cat("Warning: Too few samples (", ncol(pseudobulk_data), ") for meaningful MDS analysis. Skipping...\n")
    return(NULL)
  }
  
  # Transpose for MDS (samples as rows, genes as columns)
  mds_data <- t(pseudobulk_data)
  
  # Log2 transform with pseudocount
  mds_data_log <- log2(mds_data + 1)
  
  # Calculate distance matrix
  dist_matrix <- dist(mds_data_log)
  
  # Perform MDS
  mds_result <- cmdscale(dist_matrix, k = min(5, nrow(mds_data_log) - 1), eig = TRUE)
  
  # Create MDS dataframe
  n_dims <- ncol(mds_result$points)
  mds_df <- data.frame(mds_result$points)
  colnames(mds_df) <- paste0("MDS", 1:n_dims)
  mds_df$Replicate <- rownames(mds_data_log)
  
  # Merge with metadata
  mds_df <- mds_df %>%
    left_join(sample_metadata %>% rownames_to_column("Replicate"), by = "Replicate")
  
  # Create condition variable
  mds_df$Condition <- paste(mds_df$genotype, mds_df$Sex, "HFD", sep = "_")
  mds_df$Condition[mds_df$infection_status == "Infected"] <- paste0(mds_df$Condition[mds_df$infection_status == "Infected"], "_Infected")
  
  # Ensure we have required columns
  if (!"Sex" %in% colnames(mds_df)) mds_df$Sex <- "Unknown"
  if (!"genotype" %in% colnames(mds_df)) mds_df$genotype <- "Unknown"
  if (!"infection_status" %in% colnames(mds_df)) mds_df$infection_status <- "Unknown"
  
  # Define colors
  replicate_colors <- rainbow(nrow(mds_df))
  names(replicate_colors) <- mds_df$Replicate
  
  condition_colors <- c("KO_Female_HFD" = "#3182BD", "KO_Female_HFD_Infected" = "#08306B",
                       "KO_Male_HFD" = "#238B45", "WT_Female_HFD" = "#DE2D26", 
                       "WT_Female_HFD_Infected" = "#67000D", "WT_Male_HFD" = "#E34A33")
  
  # ===============================================
  # BASIC MDS PLOTS (5 STANDARD PLOTS)
  # ===============================================
  
  # 1. MDS plot with each replicate as a point, colored by replicate
  p1 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Replicate)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = Replicate), hjust = 0.5, vjust = -0.7, size = 3, fontface = "bold") +
    scale_color_manual(values = replicate_colors) +
    labs(title = paste("MDS of Biological Replicates:", output_prefix),
         x = "First MDS dimension",
         y = "Second MDS dimension") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "none") +
    expand_limits(x = c(min(mds_df$MDS1) * 1.2, max(mds_df$MDS1) * 1.2),
                  y = c(min(mds_df$MDS2) * 1.2, max(mds_df$MDS2) * 1.2))
  
  # 2. MDS colored by condition
  p2 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Condition)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = gsub(".*_([123])$", "\\1", Replicate)), hjust = 0.5, vjust = -0.7, 
              size = 3, fontface = "bold", color = "black") +
    scale_color_manual(values = condition_colors) +
    labs(title = paste("MDS by Condition:", output_prefix),
         x = "First MDS dimension",
         y = "Second MDS dimension",
         color = "Condition") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right") +
    expand_limits(x = c(min(mds_df$MDS1) * 1.2, max(mds_df$MDS1) * 1.2),
                  y = c(min(mds_df$MDS2) * 1.2, max(mds_df$MDS2) * 1.2))
  
  # 3. MDS colored by genotype with sex as shape
  p3 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = genotype, shape = Sex)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = gsub(".*_([123])$", "\\1", Replicate)), hjust = 0.5, vjust = -0.7, 
              size = 3, fontface = "bold", color = "black") +
    scale_color_manual(values = c("KO" = "#2166AC", "WT" = "#D6604D")) +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
    labs(title = paste("MDS by Genotype and Sex:", output_prefix),
         x = "First MDS dimension",
         y = "Second MDS dimension") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right") +
    expand_limits(x = c(min(mds_df$MDS1) * 1.2, max(mds_df$MDS1) * 1.2),
                  y = c(min(mds_df$MDS2) * 1.2, max(mds_df$MDS2) * 1.2))
  
  # 4. MDS colored by infection status
  p4 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = infection_status, shape = genotype)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = gsub(".*_([123])$", "\\1", Replicate)), hjust = 0.5, vjust = -0.7, 
              size = 3, fontface = "bold", color = "black") +
    scale_color_manual(values = c("Not_Infected" = "#4DAF4A", "Infected" = "#E41A1C")) +
    scale_shape_manual(values = c("KO" = 16, "WT" = 17)) +
    labs(title = paste("MDS by Infection Status:", output_prefix),
         x = "First MDS dimension",
         y = "Second MDS dimension") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right") +
    expand_limits(x = c(min(mds_df$MDS1) * 1.2, max(mds_df$MDS1) * 1.2),
                  y = c(min(mds_df$MDS2) * 1.2, max(mds_df$MDS2) * 1.2))
  
  # 5. MDS1 vs MDS3 plot (if MDS3 exists)
  if (n_dims >= 3) {
    p5 <- ggplot(mds_df, aes(x = MDS1, y = MDS3, color = Condition)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_text(aes(label = gsub(".*_([123])$", "\\1", Replicate)), hjust = 0.5, vjust = -0.7, 
                size = 3, fontface = "bold", color = "black") +
      scale_color_manual(values = condition_colors) +
      labs(title = paste("MDS (MDS1 vs MDS3):", output_prefix),
           x = "First MDS dimension",
           y = "Third MDS dimension",
           color = "Condition") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right") +
      expand_limits(x = c(min(mds_df$MDS1) * 1.2, max(mds_df$MDS1) * 1.2),
                    y = c(min(mds_df$MDS3) * 1.2, max(mds_df$MDS3) * 1.2))
  }
  
  # Save basic plots
  ggsave(filename = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_BiologicalReplicates_as_points.pdf")),
         plot = p1, width = 12, height = 8, dpi = 300)
  
  ggsave(filename = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_BiologicalReplicates_by_condition.pdf")),
         plot = p2, width = 12, height = 8, dpi = 300)
  
  ggsave(filename = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_BiologicalReplicates_by_genotype_sex.pdf")),
         plot = p3, width = 12, height = 8, dpi = 300)
  
  ggsave(filename = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_BiologicalReplicates_by_infection.pdf")),
         plot = p4, width = 12, height = 8, dpi = 300)
  
  if (n_dims >= 3) {
    ggsave(filename = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_BiologicalReplicates_MDS1_MDS3.pdf")),
           plot = p5, width = 12, height = 8, dpi = 300)
  }
  
  # ===============================================
  # CLUSTERING ANALYSIS ON MDS COORDINATES
  # ===============================================
  
  # Prepare data for clustering (using MDS1 and MDS2)
  cluster_data <- mds_df[, c("MDS1", "MDS2")]
  
  # ===============================================
  # HIERARCHICAL CLUSTERING WITH DIFFERENT K VALUES
  # ===============================================
  
  max_k <- min(5, nrow(mds_df) - 1)
  for (k in 2:max_k) {
    tryCatch({
      dist_matrix_mds <- dist(cluster_data)
      hc <- hclust(dist_matrix_mds, method = "ward.D2")
      cluster_assignments <- cutree(hc, k = k)
      
      mds_df$Cluster <- as.factor(cluster_assignments)
      
      mds_plot_k <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
        stat_ellipse(aes(group = Cluster, color = Cluster),
                     type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
        geom_point(aes(
          color = Sex,
          shape = genotype,
          fill = infection_status
        ), size = 4, stroke = 1.5) +
        scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
        scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
        scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
        guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "grey")),
          color = guide_legend(override.aes = list(shape = 16))
        ) +
        labs(
          x = "First MDS dimension",
          y = "Second MDS dimension",
          title = paste("MDS:", output_prefix, "(k =", k, "clusters)"),
          subtitle = "Dashed ellipses represent similarity-based clusters",
          color = "Sex / Cluster",
          shape = "Genotype", 
          fill = "Infection Status"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right"
        )
      
      ggsave(
        filename = file.path(mds_plots_dir, paste0(output_prefix, "_mds_k", k, "_clusters.pdf")),
        plot = mds_plot_k,
        width = 10, height = 8
      )
      
      cat("Saved MDS plot with k =", k, "clusters\n")
    }, error = function(e) {
      cat("Error with k =", k, ":", e$message, "\n")
    })
  }
  
  # ===============================================
  # DIFFERENT HIERARCHICAL CLUSTERING METHODS
  # ===============================================
  
  methods <- c("ward.D2", "average", "complete", "single")
  
  for (m in methods) {
    tryCatch({
      hc <- hclust(dist(cluster_data), method = m)
      cluster_assignments <- cutree(hc, k = 3)
      mds_df$Cluster <- as.factor(cluster_assignments)
      
      p <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
        stat_ellipse(aes(group = Cluster, color = Cluster),
                     type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
        geom_point(aes(color = Sex, shape = genotype, fill = infection_status),
                   size = 4, stroke = 1.5) +
        scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
        scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
        scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
        guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "grey")),
          color = guide_legend(override.aes = list(shape = 16))
        ) +
        labs(
          x = "First MDS dimension",
          y = "Second MDS dimension",
          title = paste("MDS:", output_prefix, "with", m, "clustering"),
          subtitle = "Dashed ellipses show similarity-based clusters",
          color = "Sex / Cluster",
          shape = "Genotype", 
          fill = "Infection Status"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "right"
        )
      
      ggsave(file.path(mds_plots_dir, paste0(output_prefix, "_mds_cluster_", m, ".pdf")), 
             p, width = 10, height = 8)
      cat("Saved:", m, "clustering\n")
    }, error = function(e) {
      cat("Error with", m, "method:", e$message, "\n")
    })
  }
  
  # ===============================================
  # K-MEANS CLUSTERING
  # ===============================================
  
  tryCatch({
    set.seed(123)
    k_centers <- min(3, nrow(mds_df) - 1)
    kmeans_result <- kmeans(cluster_data, centers = k_centers, nstart = 25)
    mds_df$Cluster <- as.factor(kmeans_result$cluster)
    
    mds_plot_kmeans <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
      stat_ellipse(aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = "First MDS dimension",
        y = "Second MDS dimension",
        title = paste("MDS:", output_prefix, "with K-means Clustering"),
        subtitle = paste("Dashed ellipses = k-means clusters (k=", k_centers, ")"),
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(mds_plots_dir, paste0(output_prefix, "_mds_kmeans_clusters.pdf")),
           plot = mds_plot_kmeans, width = 10, height = 8)
    cat("Saved K-means clustering\n")
  }, error = function(e) {
    cat("Error with K-means:", e$message, "\n")
  })
  
  # ===============================================
  # DBSCAN CLUSTERING
  # ===============================================
  
  tryCatch({
    db <- dbscan::dbscan(cluster_data, eps = 0.5, minPts = 2)
    mds_df$Cluster <- as.factor(db$cluster)
    
    mds_df$Cluster <- factor(mds_df$Cluster, 
                            levels = sort(unique(db$cluster)), 
                            labels = ifelse(sort(unique(db$cluster)) == 0, "Noise", 
                                           paste0("Cluster_", sort(unique(db$cluster)))))
    
    mds_plot_dbscan <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
      stat_ellipse(data = subset(mds_df, Cluster != "Noise"), 
                   aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = "First MDS dimension",
        y = "Second MDS dimension",
        title = paste("MDS:", output_prefix, "with DBSCAN Clustering"),
        subtitle = "DBSCAN: eps = 0.5, minPts = 2 | 'Noise' = unclustered samples",
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(mds_plots_dir, paste0(output_prefix, "_mds_dbscan_clusters.pdf")),
           plot = mds_plot_dbscan, width = 10, height = 8)
    cat("Saved DBSCAN clustering\n")
  }, error = function(e) {
    cat("Error with DBSCAN:", e$message, "\n")
  })
  
  # ===============================================
  # GMM CLUSTERING
  # ===============================================
  
  tryCatch({
    # Force fewer clusters to avoid over-clustering
    max_clusters <- min(4, nrow(mds_df) - 1)
    gmm <- Mclust(cluster_data, G = 2:max_clusters)
    mds_df$Cluster <- as.factor(gmm$classification)
    levels(mds_df$Cluster) <- paste0("Cluster_", levels(mds_df$Cluster))
    
    mds_plot_gmm <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
      stat_ellipse(aes(group = Cluster, color = Cluster),
                   type = "norm", level = 0.68, linetype = "dashed", size = 0.7, alpha = 0.7) +
      geom_point(aes(
        color = Sex,
        shape = genotype,
        fill = infection_status
      ), size = 4, stroke = 1.5) +
      scale_color_manual(values = c("Male" = "#6CA6CD", "Female" = "#D82020")) +
      scale_shape_manual(values = c("WT" = 22, "KO" = 24)) +
      scale_fill_manual(values = c("Not_Infected" = "white", "Infected" = "black")) +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = guide_legend(override.aes = list(fill = "grey")),
        color = guide_legend(override.aes = list(shape = 16))
      ) +
      labs(
        x = "First MDS dimension",
        y = "Second MDS dimension",
        title = paste("MDS:", output_prefix, "with GMM Clustering (mclust)"),
        subtitle = paste("Model:", gmm$modelName, "| Selected clusters:", gmm$G),
        color = "Sex / Cluster",
        shape = "Genotype",
        fill = "Infection Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    ggsave(file.path(mds_plots_dir, paste0(output_prefix, "_mds_gmm_clusters.pdf")),
           plot = mds_plot_gmm, width = 10, height = 8)
    cat("Saved GMM clustering\n")
  }, error = function(e) {
    cat("Error with GMM:", e$message, "\n")
  })
  
  # ===============================================
  # SAVE COORDINATES AND ANALYSIS DATA
  # ===============================================
  
  # Calculate distance analysis
  same_condition_pairs <- data.frame(
    Replicate1 = character(),
    Replicate2 = character(),
    Distance = numeric(),
    Same_Condition = logical()
  )
  
  dist_df <- as.matrix(dist_matrix)
  for (i in 1:(nrow(mds_df)-1)) {
    for (j in (i+1):nrow(mds_df)) {
      rep1 <- mds_df$Replicate[i]
      rep2 <- mds_df$Replicate[j]
      same_cond <- mds_df$Condition[i] == mds_df$Condition[j]
      
      same_condition_pairs <- rbind(same_condition_pairs, data.frame(
        Replicate1 = rep1,
        Replicate2 = rep2,
        Distance = dist_df[rep1, rep2],
        Same_Condition = same_cond
      ))
    }
  }
  
  # Save coordinates and distance analysis
  write.csv(mds_df, file = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_biological_replicate_coordinates.csv")), row.names = FALSE)
  write.csv(same_condition_pairs, file = file.path(mds_plots_dir, paste0(output_prefix, "_MDS_replicate_distances.csv")), row.names = FALSE)
  
  cat("Completed MDS analysis for:", output_prefix, "\n")
}

# ===============================================
# DEFINE COMBINED EOS+MAC DATASETS TO ANALYZE
# ===============================================

datasets <- list(
  # Combined Eosinophils + Macrophages
  list(pseudobulk = "pseudobulk_eos_mac_all_samples.rds", 
       metadata = "sample_metadata_eos_mac_all_samples.rds", 
       prefix = "eos_mac_all_samples"),
  list(pseudobulk = "pseudobulk_eos_mac_not_infected.rds", 
       metadata = "sample_metadata_eos_mac_not_infected.rds", 
       prefix = "eos_mac_not_infected"),
  list(pseudobulk = "pseudobulk_eos_mac_infected.rds", 
       metadata = "sample_metadata_eos_mac_infected.rds", 
       prefix = "eos_mac_infected")
)

# ===============================================
# RUN ANALYSIS FOR COMBINED DATASETS
# ===============================================

cat("Starting MDS analysis for combined eosinophils and macrophages datasets...\n")
cat("Output directory:", mds_plots_dir, "\n")

for (dataset in datasets) {
  perform_mds_analysis(dataset$pseudobulk, dataset$metadata, dataset$prefix)
}

cat("\n===============================================\n")
cat("MDS ANALYSIS COMPLETE!\n")
cat("===============================================\n")
cat("All MDS plots saved to:", mds_plots_dir, "\n")
cat("Total datasets analyzed:", length(datasets), "\n")
cat("Each dataset includes:\n")
cat("- 5 basic MDS visualization plots\n")
cat("- Hierarchical clustering (k=2 to 5)\n")
cat("- Different linkage methods (ward.D2, average, complete, single)\n")
cat("- K-means clustering\n")
cat("- DBSCAN clustering\n")
cat("- GMM clustering (mclust) with forced fewer clusters\n")
cat("- Distance analysis and coordinates saved as CSV files\n")

# ===============================================
# SUMMARY OF FILES GENERATED PER DATASET
# ===============================================

cat("\nFile naming convention per dataset:\n")
cat("Basic MDS plots:\n")
cat("  - {prefix}_MDS_BiologicalReplicates_as_points.pdf\n")
cat("  - {prefix}_MDS_BiologicalReplicates_by_condition.pdf\n")
cat("  - {prefix}_MDS_BiologicalReplicates_by_genotype_sex.pdf\n")
cat("  - {prefix}_MDS_BiologicalReplicates_by_infection.pdf\n")
cat("  - {prefix}_MDS_BiologicalReplicates_MDS1_MDS3.pdf\n")

cat("\nClustering plots:\n")
cat("  - {prefix}_mds_k2_clusters.pdf through {prefix}_mds_k5_clusters.pdf\n")
cat("  - {prefix}_mds_cluster_ward.D2.pdf\n")
cat("  - {prefix}_mds_cluster_average.pdf\n")
cat("  - {prefix}_mds_cluster_complete.pdf\n")
cat("  - {prefix}_mds_cluster_single.pdf\n")
cat("  - {prefix}_mds_kmeans_clusters.pdf\n")
cat("  - {prefix}_mds_dbscan_clusters.pdf\n")
cat("  - {prefix}_mds_gmm_clusters.pdf\n")

cat("\nData files:\n")
cat("  - {prefix}_MDS_biological_replicate_coordinates.csv\n")
cat("  - {prefix}_MDS_replicate_distances.csv\n")

cat("\nEstimated total files for combined datasets: ~", length(datasets) * 15, "plots + ", length(datasets) * 2, "CSV files\n")

cat("\n===============================================\n")
cat("EXPECTED BENEFITS OF COMBINED ANALYSIS:\n")
cat("===============================================\n")
cat("1. INCREASED STATISTICAL POWER:\n")
cat("   - More cells per sample = more robust expression profiles\n")
cat("   - Better clustering due to reduced noise\n")
cat("   - Clearer ellipses in all clustering plots\n\n")

cat("2. IMMUNE CELL FOCUSED ANALYSIS:\n")
cat("   - Specific focus on eosinophils and macrophages\n")
cat("   - Cross-cell-type response patterns\n")
cat("   - Coordinated immune cell behavior\n\n")

cat("3. ENHANCED PATTERN DETECTION:\n")
cat("   - Better separation of experimental conditions\n")
cat("   - More reliable genotype and infection effects\n")
cat("   - Reduced individual cell type variability\n\n")

cat("OUTPUT DIRECTORIES:\n")
cat("- PCA plots: /bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/pca_plots_eos_mac/\n")
cat("- MDS plots: /bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/mds_plots_eos_mac/\n")
cat("- Original data: /bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/\n")