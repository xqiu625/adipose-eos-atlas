# MDS and Clustering Analysis for All Pseudobulk Datasets
library(ggplot2)
library(dplyr)
library(tibble)
library(dbscan)
library(mclust)

# Set directories
out_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster'
mds_plots_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/05_pseudobulk_cluster/mds_plots'
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
    gmm <- Mclust(cluster_data)
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
# DEFINE ALL DATASETS TO ANALYZE
# ===============================================

datasets <- list(
  # All cells
  list(pseudobulk = "pseudobulk_all_cells_all_samples.rds", 
       metadata = "sample_metadata_all_cells_all_samples.rds", 
       prefix = "all_cells_all_samples"),
  list(pseudobulk = "pseudobulk_all_cells_not_infected.rds", 
       metadata = "sample_metadata_all_cells_not_infected.rds", 
       prefix = "all_cells_not_infected"),
  list(pseudobulk = "pseudobulk_all_cells_infected.rds", 
       metadata = "sample_metadata_all_cells_infected.rds", 
       prefix = "all_cells_infected"),
  
  # Eosinophils
  list(pseudobulk = "pseudobulk_eosinophils_all_samples.rds", 
       metadata = "sample_metadata_eosinophils_all_samples.rds", 
       prefix = "eosinophils_all_samples"),
  list(pseudobulk = "pseudobulk_eosinophils_not_infected.rds", 
       metadata = "sample_metadata_eosinophils_not_infected.rds", 
       prefix = "eosinophils_not_infected"),
  list(pseudobulk = "pseudobulk_eosinophils_infected.rds", 
       metadata = "sample_metadata_eosinophils_infected.rds", 
       prefix = "eosinophils_infected"),
  
  # Macrophages
  list(pseudobulk = "pseudobulk_macrophages_all_samples.rds", 
       metadata = "sample_metadata_macrophages_all_samples.rds", 
       prefix = "macrophages_all_samples"),
  list(pseudobulk = "pseudobulk_macrophages_not_infected.rds", 
       metadata = "sample_metadata_macrophages_not_infected.rds", 
       prefix = "macrophages_not_infected"),
  list(pseudobulk = "pseudobulk_macrophages_infected.rds", 
       metadata = "sample_metadata_macrophages_infected.rds", 
       prefix = "macrophages_infected")
)

# ===============================================
# RUN ANALYSIS FOR ALL DATASETS
# ===============================================

cat("Starting comprehensive MDS analysis for all datasets...\n")
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
cat("- GMM clustering (mclust)\n")
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

cat("\nEstimated total files: ~", length(datasets) * 15, "plots + ", length(datasets) * 2, "CSV files\n")