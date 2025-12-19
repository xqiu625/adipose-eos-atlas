# PCA and Correlation Circle Analysis for Pseudobulk Data
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(dplyr)
library(tibble)

# Set directories
out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'

# ===============================================
# LOAD DATA
# ===============================================

print("Loading pseudobulk data...")
pseudobulk_all <- readRDS(file.path(out_dir, "pseudobulk_all_celltypes.rds"))
sample_metadata <- readRDS(file.path(out_dir, "sample_metadata_all.rds"))

print(paste("Loaded data:", nrow(pseudobulk_all), "genes x", ncol(pseudobulk_all), "samples"))
print("Sample metadata:")
print(sample_metadata)

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
variance_explained <- summary(pca_result)$importance[2, 1:10] * 100

print("Variance explained by first 10 PCs:")
print(round(variance_explained, 2))

# ===============================================
# PREPARE PLOT DATA
# ===============================================

# Extract PC scores
pc_scores <- data.frame(pca_result$x[, 1:10])
pc_scores$sample <- rownames(pc_scores)

# Merge with metadata
plot_data <- pc_scores %>%
  rownames_to_column("sample_id") %>%
  left_join(sample_metadata %>% rownames_to_column("sample_id"), by = "sample_id")

# Create plotting variables
plot_data$infection_fill <- ifelse(plot_data$infection_status == "Infected", "filled", "clear")
plot_data$sex_color <- ifelse(plot_data$Sex == "Male",  "#6CA6CD", "#D82020")
plot_data$genotype_shape <- ifelse(plot_data$genotype == "WT", "square", "triangle")

# Convert to factors for ggplot
plot_data$infection_status <- factor(plot_data$infection_status, levels = c("Not_Infected", "Infected"))
plot_data$Sex <- factor(plot_data$Sex)
plot_data$genotype <- factor(plot_data$genotype, levels = c("WT", "KO"))

print("Plot data prepared:")
print(table(plot_data$infection_status))
print(table(plot_data$Sex))
print(table(plot_data$genotype))


# Try multiple cluster numbers (e.g., 2 to 5)
for (k in 2:5) {
  
  # Perform hierarchical clustering again for consistency
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
      title = paste("PCA of Pseudobulk Samples (k =", k, "clusters)"),
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
    filename = file.path(out_dir, paste0("pca_plot_k", k, "_clusters.pdf")),
    plot = pca_plot_k,
    width = 10, height = 8
  )
  
  cat(paste("Saved PCA plot with k =", k, "clusters\n"))
}


methods <- c("ward.D2", "average", "complete", "single")

for (m in methods) {
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
    labs(title = paste("PCA with", m, "clustering"),
         subtitle = "Dashed ellipses show similarity-based clusters") +
    theme_minimal()
  
  ggsave(file.path(out_dir, paste0("pca_cluster_", m, ".pdf")), p, width = 10, height = 8)
  cat("Saved:", m, "\n")
}

# ===============================================
# KMEANS CLUSTERING
# ===============================================

# Run k-means clustering (already done)
set.seed(123)
kmeans_result <- kmeans(plot_data[, c("PC1", "PC2")], centers = 3, nstart = 25)
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
    title = "PCA of Pseudobulk Samples with K-means Clustering",
    subtitle = "Dashed ellipses = k-means clusters (k=3)",
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
ggsave(file.path(out_dir, "pca_plot_kmeans_clusters.pdf"),
       plot = pca_plot_kmeans, width = 10, height = 8)

cat("Saved PCA plot with k-means clustering to:", file.path(out_dir, "pca_plot_kmeans_clusters.pdf"), "\n")


# ===============================================
# DBSCAN CLUSTERING
# ===============================================

# Run DBSCAN clustering
library(dbscan)
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
    title = "PCA of Pseudobulk Samples with DBSCAN Clustering",
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
ggsave(file.path(out_dir, "pca_plot_dbscan_clusters.pdf"),
       plot = pca_plot_dbscan, width = 10, height = 8)

cat("Saved PCA plot with DBSCAN clustering to:", file.path(out_dir, "pca_plot_dbscan_clusters.pdf"), "\n")


# ===============================================
# GMM CLUSTERING
# ===============================================

# Load library and run GMM clustering
library(mclust)
gmm <- Mclust(plot_data[, c("PC1", "PC2")])
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
    title = "PCA of Pseudobulk Samples with GMM Clustering (mclust)",
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
ggsave(file.path(out_dir, "pca_plot_gmm_clusters.pdf"),
       plot = pca_plot_gmm, width = 10, height = 8)

cat("Saved PCA plot with GMM (mclust) clustering to:", file.path(out_dir, "pca_plot_gmm_clusters.pdf"), "\n")