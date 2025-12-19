library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# Load directories
DATA_DIR <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed"
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R03_draft_v2/Figure2/CellChat'

# Load the Eosinophils and Macrophages subset
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_object_with_mac_subtypes.rds"))

eos_seurat <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")

# Assuming seurat_object and eos_seurat are both loaded
# First, extract the eos_cluster information from eos_seurat
eos_clusters <- eos_seurat@meta.data[, c("eos_cluster")]

# Create a named vector for easier mapping
eos_clusters_named <- setNames(eos_clusters, rownames(eos_seurat@meta.data))

# Add eos_cluster to the main seurat_object
# This will add the column only for cells that exist in both objects
# For cells in seurat_object that don't exist in eos_seurat, the value will be NA
seurat_object@meta.data$eos_cluster <- eos_clusters_named[rownames(seurat_object@meta.data)]

# If you only want to add this information to cells identified as Eosinophils
# You can do:
seurat_object@meta.data$eos_cluster <- NA  # Initialize with NA
eosinophil_cells <- rownames(seurat_object@meta.data)[seurat_object@meta.data$cell_type == "Eosinophils"]
matching_cells <- intersect(eosinophil_cells, rownames(eos_seurat@meta.data))
seurat_object@meta.data[matching_cells, "eos_cluster"] <- eos_seurat@meta.data[matching_cells, "eos_cluster"]


eos_seurat <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")
table(eos_seurat@meta.data$eos_cluster)

# First, create the new cell type label column that includes eosinophil clusters
seurat_object$cell_type_with_eos_clusters <- seurat_object$cell_type

# For cells that are eosinophils, append their eos_cluster information to create new labels
eos_mask <- seurat_object$cell_type == "Eosinophils"
seurat_object$cell_type_with_eos_clusters[eos_mask] <- paste0("Eosinophil_C", seurat_object$eos_cluster[eos_mask])

# Verify the new column has been created correctly
table(seurat_object$eos_cluster)

seurat_object$eos_cluster <- seurat_object$eos_cluster - 1

# Save the updated seurat object
saveRDS(seurat_object, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/seurat_object_with_eos_clusters.rds")





seurat_object <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/seurat_object_with_eos_clusters.rds")
table(seurat_object$eos_cluster)
head(seurat_object@meta.data)

# Function to analyze interactions between a specific eosinophil cluster and macrophages
analyze_eos_cluster_mac_interactions <- function(seurat_object, eos_cluster) {
    # Filter to include only:
    # 1. The specific eosinophil cluster
    # 2. Macrophage cells
    
    eos_cells <- rownames(seurat_object@meta.data)[
        seurat_object$cell_type == "Eosinophils" & 
        seurat_object$eos_cluster == eos_cluster
    ]
    
    mac_cells <- rownames(seurat_object@meta.data)[
        !is.na(seurat_object$macrophage_type) & 
        seurat_object$macrophage_type %in% c("M1_Macrophages", "M2_Macrophages_Functional", "M2_Macrophages_Surface_Only")
    ]
    
    # Only proceed if we have enough eosinophil cells in this cluster
    if(length(eos_cells) < 10) {
        cat("Not enough cells in Eosinophil Cluster", eos_cluster, "(", length(eos_cells), "cells)\n")
        return(NULL)
    }
    
    # Combine cells for analysis
    cells_to_use <- c(eos_cells, mac_cells)
    subset_for_analysis <- subset(seurat_object, cells = cells_to_use)
    
    # Create metadata for CellChat
    meta <- data.frame(
        group = case_when(
            subset_for_analysis$cell_type == "Eosinophils" ~ paste0("Eosinophil_C", eos_cluster),
            subset_for_analysis$macrophage_type == "M1_Macrophages" ~ "M1_Mac",
            subset_for_analysis$macrophage_type == "M2_Macrophages_Functional" ~ "M2_Functional",
            subset_for_analysis$macrophage_type == "M2_Macrophages_Surface_Only" ~ "M2_Surface"
        ),
        row.names = colnames(subset_for_analysis)
    )
    
    # Check if we have all macrophage types represented
    mac_types <- c("M1_Mac", "M2_Functional", "M2_Surface")
    present_types <- intersect(mac_types, unique(meta$group))
    if(length(present_types) < 3) {
        cat("Warning: Not all macrophage types are present in the analysis for Eosinophil Cluster", eos_cluster, "\n")
        cat("Present types:", paste(present_types, collapse=", "), "\n")
    }
    
    # Create CellChat object and process
    data.input <- GetAssayData(subset_for_analysis, assay = "RNA", slot = "data")
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
    cellchat@DB <- CellChatDB.mouse
    
    # Process CellChat object
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    return(cellchat)
}

# Create directory for saving results
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R03_draft_v2/Figure2/CellChat'
dir.create(CellChat_dir, showWarnings = FALSE, recursive = TRUE)

# Process each eosinophil cluster
eos_clusters <- 0:6  # Clusters 0 through 6

for(eos_clust in eos_clusters) {
    cat("Analyzing Eosinophil Cluster", eos_clust, "...\n")
    
    cellchat_obj <- analyze_eos_cluster_mac_interactions(seurat_object, eos_clust)
    
    if(!is.null(cellchat_obj)) {
        # Save the CellChat object
        save_path <- file.path(CellChat_dir, paste0("EosC", eos_clust, "_mac_cellchat.rds"))
        saveRDS(cellchat_obj, save_path)
        cat("  Saved to:", save_path, "\n")
        
        # Create and save a basic circle plot to visualize interactions
        pdf(file.path(CellChat_dir, paste0("EosC", eos_clust, "_mac_circle.pdf")), width = 8, height = 8)
        # Use cellchat_obj instead of cellchat here
        netVisual_circle(cellchat_obj@net$weight, 
                        vertex.weight = rowSums(cellchat_obj@net$weight), 
                        weight.scale = TRUE, 
                        edge.weight.max = max(cellchat_obj@net$weight), 
                        title.name = paste0("Eosinophil Cluster ", eos_clust, " Interactions with Macrophages"), 
                        arrow.size = 0.3,
                        arrow.width = 1,
                        label.edge = TRUE)
        dev.off()
    }
}

cat("All CellChat analyses complete!\n")