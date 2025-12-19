library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)

# Load directories
DATA_DIR <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed"
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat'
# Create CellChat_dir if it does not exist
if (!dir.exists(CellChat_dir)) {
  dir.create(CellChat_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the Eosinophils and Macrophages subset
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_object_with_mac_subtypes_251021.rds"))
eos_object <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/seurat_object_for_cellchat.rds")
eos_object_metadata <- eos_object@meta.data
eos_object_eos_sub <- eos_object_metadata %>% select(eos_cluster)

seurat_object@meta.data$eos_cluster <- eos_object_eos_sub[rownames(seurat_object@meta.data), "eos_cluster"]


# Function 1: Analyze eosinophil subcluster interactions with macrophage SUBTYPES
analyze_eos_mac_subtypes <- function(seurat_subset, eos_subcluster_name) {
    # Create CellChat object and process
    data.input <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
    
    meta <- data.frame(
        group = case_when(
            !is.na(seurat_subset$eos_cluster) & seurat_subset$eos_cluster == eos_subcluster_name ~ paste0("Eos_", eos_subcluster_name),
            seurat_subset$macrophage_type == "M1_Macrophages" ~ "M1_Mac",
            seurat_subset$macrophage_type == "M2_Macrophages_Functional" ~ "M2_Functional",
            seurat_subset$macrophage_type == "M2_Macrophages_Surface_Only" ~ "M2_Surface",
            seurat_subset$macrophage_type == "Other_Macrophages" ~ "Other_Mac",
            TRUE ~ NA_character_
        ),
        row.names = colnames(seurat_subset)
    )
    
    # Remove NA values
    meta <- meta[!is.na(meta$group), , drop = FALSE]
    
    # Drop unused factor levels
    meta$group <- droplevels(as.factor(meta$group))
    
    # Print cell type composition
    cat("Cell type composition (with Mac Subtypes):\n")
    print(table(meta$group))
    cat("\n")
    
    # Check if we have enough cells in each group
    group_counts <- table(meta$group)
    if(any(group_counts < 10)) {
        cat("WARNING: Some groups have fewer than 10 cells:\n")
        print(group_counts[group_counts < 10])
        cat("These groups may not have reliable communication inference.\n\n")
    }
    
    # Check if we have at least 2 cell groups
    if(length(unique(meta$group)) < 2) {
        cat("ERROR: Need at least 2 cell groups for CellChat analysis. Skipping", eos_subcluster_name, "\n")
        return(NULL)
    }
    
    # Filter the data.input to match the filtered meta
    data.input <- data.input[, rownames(meta)]
    
    tryCatch({
        cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
        cellchat@DB <- CellChatDB.mouse
        
        # Process CellChat object
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)
        
        return(cellchat)
    }, error = function(e) {
        cat("ERROR in processing", eos_subcluster_name, ":\n")
        cat(as.character(e), "\n")
        return(NULL)
    })
}


# Function 2: Analyze eosinophil subcluster interactions with ALL macrophages (combined)
analyze_eos_all_mac <- function(seurat_subset, eos_subcluster_name) {
    # Create CellChat object and process
    data.input <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
    
    meta <- data.frame(
        group = case_when(
            !is.na(seurat_subset$eos_cluster) & seurat_subset$eos_cluster == eos_subcluster_name ~ paste0("Eos_", eos_subcluster_name),
            seurat_subset$cell_type == "Macrophages" ~ "All_Macrophages",
            TRUE ~ NA_character_
        ),
        row.names = colnames(seurat_subset)
    )
    
    # Remove NA values
    meta <- meta[!is.na(meta$group), , drop = FALSE]
    
    # Drop unused factor levels
    meta$group <- droplevels(as.factor(meta$group))
    
    # Print cell type composition
    cat("Cell type composition (with All Macrophages):\n")
    print(table(meta$group))
    cat("\n")
    
    # Check if we have enough cells in each group
    group_counts <- table(meta$group)
    if(any(group_counts < 10)) {
        cat("WARNING: Some groups have fewer than 10 cells:\n")
        print(group_counts[group_counts < 10])
        cat("These groups may not have reliable communication inference.\n\n")
    }
    
    # Check if we have at least 2 cell groups
    if(length(unique(meta$group)) < 2) {
        cat("ERROR: Need at least 2 cell groups for CellChat analysis. Skipping", eos_subcluster_name, "\n")
        return(NULL)
    }
    
    # Filter the data.input to match the filtered meta
    data.input <- data.input[, rownames(meta)]
    
    tryCatch({
        cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
        cellchat@DB <- CellChatDB.mouse
        
        # Process CellChat object
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)
        
        return(cellchat)
    }, error = function(e) {
        cat("ERROR in processing", eos_subcluster_name, ":\n")
        cat(as.character(e), "\n")
        return(NULL)
    })
}


# Get unique eosinophil subclusters (excluding NA)
eos_subclusters <- unique(seurat_object@meta.data$eos_cluster)
eos_subclusters <- eos_subclusters[!is.na(eos_subclusters)]  # Remove NA
eos_subclusters <- droplevels(eos_subclusters)  # Drop unused levels

cat("Eosinophil subclusters to analyze:\n")
print(eos_subclusters)
cat("\n")

# Create subsets for each eosinophil subcluster + ALL macrophages
eos_subclusters_subsets <- list()
for(eos_subcluster in eos_subclusters) {
    # Get cells that are either:
    # 1. Eosinophils from this specific subcluster, OR
    # 2. Any macrophage
    cells_to_keep <- colnames(seurat_object)[
        (seurat_object$eos_cluster == eos_subcluster & !is.na(seurat_object$eos_cluster)) |
        (seurat_object$cell_type == "Macrophages")
    ]
    
    eos_subclusters_subsets[[eos_subcluster]] <- subset(seurat_object, cells = cells_to_keep)
    
    cat("Subsetting", eos_subcluster, "- Total cells:", ncol(eos_subclusters_subsets[[eos_subcluster]]), "\n")
    cat("  Eosinophils (", eos_subcluster, "):", 
        sum(eos_subclusters_subsets[[eos_subcluster]]$eos_cluster == eos_subcluster, na.rm = TRUE), "\n")
    cat("  Macrophages:", 
        sum(eos_subclusters_subsets[[eos_subcluster]]$cell_type == "Macrophages", na.rm = TRUE), "\n\n")
}


# Process each eosinophil subcluster - BOTH with Mac subtypes AND with all Mac combined
successful_subtypes <- c()
failed_subtypes <- c()
successful_all_mac <- c()
failed_all_mac <- c()

for(eos_subcluster in eos_subclusters) {
    cat("\n########################################\n")
    cat("## Processing", eos_subcluster, "\n")
    cat("########################################\n\n")
    
    # 1. Analyze with Macrophage SUBTYPES
    cat("========================================\n")
    cat("Analysis 1: ", eos_subcluster, " with Mac SUBTYPES\n")
    cat("========================================\n")
    
    cellchat_subtypes <- analyze_eos_mac_subtypes(eos_subclusters_subsets[[eos_subcluster]], eos_subcluster)
    
    if(!is.null(cellchat_subtypes)) {
        saveRDS(cellchat_subtypes, file.path(CellChat_dir, paste0("eos_", eos_subcluster, "_vs_MacSubtypes_cellchat.rds")))
        cat("✓ Successfully saved", eos_subcluster, "vs Mac Subtypes CellChat object\n\n")
        successful_subtypes <- c(successful_subtypes, eos_subcluster)
    } else {
        cat("✗ Failed to process", eos_subcluster, "vs Mac Subtypes\n\n")
        failed_subtypes <- c(failed_subtypes, eos_subcluster)
    }
    
    # 2. Analyze with ALL Macrophages (combined)
    cat("========================================\n")
    cat("Analysis 2: ", eos_subcluster, " with ALL Macrophages\n")
    cat("========================================\n")
    
    cellchat_all_mac <- analyze_eos_all_mac(eos_subclusters_subsets[[eos_subcluster]], eos_subcluster)
    
    if(!is.null(cellchat_all_mac)) {
        saveRDS(cellchat_all_mac, file.path(CellChat_dir, paste0("eos_", eos_subcluster, "_vs_AllMac_cellchat.rds")))
        cat("✓ Successfully saved", eos_subcluster, "vs All Mac CellChat object\n\n")
        successful_all_mac <- c(successful_all_mac, eos_subcluster)
    } else {
        cat("✗ Failed to process", eos_subcluster, "vs All Mac\n\n")
        failed_all_mac <- c(failed_all_mac, eos_subcluster)
    }
}


# Summary
cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n\n")

cat("Analysis 1: Eosinophil Subclusters vs Macrophage SUBTYPES\n")
cat("----------------------------------------------------------\n")
cat("Successfully processed:", length(successful_subtypes), "subclusters\n")
if(length(successful_subtypes) > 0) {
    cat("  -", paste(successful_subtypes, collapse="\n  - "), "\n")
}
cat("Failed to process:", length(failed_subtypes), "subclusters\n")
if(length(failed_subtypes) > 0) {
    cat("  -", paste(failed_subtypes, collapse="\n  - "), "\n")
}

cat("\n")

cat("Analysis 2: Eosinophil Subclusters vs ALL Macrophages (combined)\n")
cat("----------------------------------------------------------\n")
cat("Successfully processed:", length(successful_all_mac), "subclusters\n")
if(length(successful_all_mac) > 0) {
    cat("  -", paste(successful_all_mac, collapse="\n  - "), "\n")
}
cat("Failed to process:", length(failed_all_mac), "subclusters\n")
if(length(failed_all_mac) > 0) {
    cat("  -", paste(failed_all_mac, collapse="\n  - "), "\n")
}