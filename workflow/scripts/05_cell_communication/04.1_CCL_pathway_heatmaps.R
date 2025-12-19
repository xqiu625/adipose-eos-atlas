library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set directories
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat'
output_dir <- file.path(CellChat_dir, 'CCL_signaling_role_plots')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# DIAGNOSTIC: Check individual objects first
# ============================================================================

cat("========================================\n")
cat("DIAGNOSTIC MODE\n")
cat("========================================\n\n")

cat("Loading CellChat objects...\n")

eos_clusters <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6")

eos_mac_subtypes <- list()
for(cluster in eos_clusters) {
    file_subtypes <- file.path(CellChat_dir, paste0("eos_", cluster, "_vs_MacSubtypes_cellchat.rds"))
    if(file.exists(file_subtypes)) {
        eos_mac_subtypes[[cluster]] <- readRDS(file_subtypes)
        cat("  Loaded: Eos", cluster, "\n")
    }
}

cat("\nTotal loaded:", length(eos_mac_subtypes), "\n\n")

# Check first object structure
cat("Checking first object (C0) structure:\n")
obj <- eos_mac_subtypes[[1]]
cat("  Class:", class(obj), "\n")
cat("  Slots:", slotNames(obj), "\n")
cat("  Cell types:", levels(obj@idents), "\n")
cat("  Number of pathways:", length(obj@netP$pathways), "\n")
cat("  Has centrality?:", !is.null(obj@netP$centr), "\n")
if(!is.null(obj@netP$centr)) {
    cat("  Centrality pathways:", names(obj@netP$centr), "\n")
}
cat("\n")

# Check for CCL
ccl_pathways <- c()
for(cluster in names(eos_mac_subtypes)) {
    pathways <- eos_mac_subtypes[[cluster]]@netP$pathways
    ccl_matches <- pathways[grepl("^CCL", pathways, ignore.case = TRUE)]
    ccl_pathways <- union(ccl_pathways, ccl_matches)
}

cat("CCL pathways found:", paste(ccl_pathways, collapse=", "), "\n")
if(length(ccl_pathways) == 0) {
    stop("No CCL pathways found")
}
cat("\n")

# ============================================================================
# OPTION 1: Try merging WITHOUT cell.prefix (might preserve centrality better)
# ============================================================================

cat("========================================\n")
cat("OPTION 1: Merging WITHOUT cell.prefix\n")
cat("========================================\n\n")

cellchat_merged_v1 <- NULL

tryCatch({
    # Don't use cell.prefix this time
    cellchat_merged_v1 <- mergeCellChat(
        eos_mac_subtypes, 
        add.names = paste0("Eos_", names(eos_mac_subtypes)),
        cell.prefix = FALSE  # Changed to FALSE
    )
    cat("✓ Merged successfully without cell.prefix\n")
    
    # Check pathways
    if(is.null(cellchat_merged_v1@netP$pathways)) {
        cat("⚠ Pathways are NULL - reconstructing...\n")
        all_pathways <- c()
        for(obj in eos_mac_subtypes) {
            all_pathways <- union(all_pathways, obj@netP$pathways)
        }
        cellchat_merged_v1@netP$pathways <- all_pathways
    }
    
    cat("  Pathways:", length(cellchat_merged_v1@netP$pathways), "\n")
    cat("  Has centrality?:", !is.null(cellchat_merged_v1@netP$centr), "\n")
    
    # Try to compute centrality
    cat("\nComputing centrality...\n")
    cellchat_merged_v1 <- netAnalysis_computeCentrality(
        cellchat_merged_v1, 
        slot.name = "netP"
    )
    cat("✓ Centrality computed!\n")
    
    # Check centrality
    ccl_in_merged <- intersect(ccl_pathways, cellchat_merged_v1@netP$pathways)
    if(length(ccl_in_merged) > 0) {
        ccl_pathway <- ccl_in_merged[1]
        if(!is.null(cellchat_merged_v1@netP$centr[[ccl_pathway]])) {
            cat("✓ Centrality data exists for", ccl_pathway, "\n")
            
            # Try to plot
            cat("\nAttempting to plot...\n")
            pdf(file.path(output_dir, "CCL_OPTION1_no_prefix.pdf"), width = 16, height = 10)
            
            tryCatch({
                gg <- netAnalysis_signalingRole_network(
                    cellchat_merged_v1, 
                    signaling = ccl_pathway,
                    width = 14, 
                    height = 2.5, 
                    font.size = 10
                )
                print(gg)
                cat("✓✓✓ PLOT SUCCESSFUL!\n")
            }, error = function(e) {
                cat("✗ Plot failed:", as.character(e), "\n")
            })
            
            dev.off()
        }
    }
    
}, error = function(e) {
    cat("✗ Option 1 failed:", as.character(e), "\n")
})

cat("\n")

# ============================================================================
# OPTION 2: Don't merge - use individual objects and combine plots manually
# ============================================================================

cat("========================================\n")
cat("OPTION 2: Generate plots from individual objects\n")
cat("========================================\n\n")

cat("Since merging causes issues, let's work with individual objects\n")
cat("We'll create a plot for each Eos subcluster\n\n")

# For each individual object, generate the plot
for(cluster in names(eos_mac_subtypes)) {
    cat("Processing", cluster, "...\n")
    
    obj <- eos_mac_subtypes[[cluster]]
    
    # Check for CCL in this object
    ccl_in_obj <- intersect(ccl_pathways, obj@netP$pathways)
    
    if(length(ccl_in_obj) == 0) {
        cat("  No CCL pathway\n\n")
        next
    }
    
    ccl_pathway <- ccl_in_obj[1]
    
    # Check centrality
    if(is.null(obj@netP$centr) || is.null(obj@netP$centr[[ccl_pathway]])) {
        cat("  Computing centrality...\n")
        obj <- netAnalysis_computeCentrality(obj, slot.name = "netP")
    }
    
    # Generate plot
    pdf_file <- file.path(output_dir, paste0("CCL_individual_", cluster, "_detailed.pdf"))
    pdf(pdf_file, width = 12, height = 8)
    
    tryCatch({
        gg <- netAnalysis_signalingRole_network(
            obj, 
            signaling = ccl_pathway,
            width = 10, 
            height = 2, 
            font.size = 10
        )
        print(gg)
        cat("  ✓ Plot generated\n")
    }, error = function(e) {
        cat("  ✗ Error:", as.character(e), "\n")
    })
    
    dev.off()
    
    # Also save centrality data
    if(!is.null(obj@netP$centr[[ccl_pathway]])) {
        write.csv(obj@netP$centr[[ccl_pathway]], 
                 file.path(output_dir, paste0("CCL_centrality_", cluster, ".csv")),
                 row.names = TRUE)
        cat("  ✓ Centrality data saved\n")
    }
    
    cat("\n")
}

# ============================================================================
# OPTION 3: Manual combination of centrality data
# ============================================================================

cat("========================================\n")
cat("OPTION 3: Manually combine centrality data\n")
cat("========================================\n\n")

cat("Extracting centrality data from all individual objects...\n")

all_centrality <- list()

for(cluster in names(eos_mac_subtypes)) {
    obj <- eos_mac_subtypes[[cluster]]
    
    # Ensure centrality is computed
    if(is.null(obj@netP$centr)) {
        obj <- netAnalysis_computeCentrality(obj, slot.name = "netP")
    }
    
    # Get CCL pathway
    ccl_in_obj <- intersect(ccl_pathways, obj@netP$pathways)
    if(length(ccl_in_obj) == 0) next
    
    ccl_pathway <- ccl_in_obj[1]
    
    # Extract centrality
    if(!is.null(obj@netP$centr[[ccl_pathway]])) {
        centr_data <- obj@netP$centr[[ccl_pathway]]
        
        # Add cluster info to rownames
        rownames(centr_data) <- paste0(cluster, "_", rownames(centr_data))
        
        all_centrality[[cluster]] <- centr_data
        cat("  Extracted from", cluster, "-", nrow(centr_data), "cell types\n")
    }
}

# Combine all centrality data
if(length(all_centrality) > 0) {
    combined_centrality <- do.call(rbind, all_centrality)
    
    cat("\n✓ Combined centrality data:\n")
    cat("  Total cell types:", nrow(combined_centrality), "\n")
    cat("  Measures:", colnames(combined_centrality), "\n")
    
    # Save combined data
    write.csv(combined_centrality, 
             file.path(output_dir, "CCL_centrality_ALL_COMBINED.csv"),
             row.names = TRUE)
    cat("  ✓ Saved: CCL_centrality_ALL_COMBINED.csv\n\n")
    
    # Create manual heatmap-like visualization
    cat("Creating manual visualization...\n")
    
    pdf(file.path(output_dir, "CCL_MANUAL_HEATMAP.pdf"), width = 16, height = 10)
    
    # Normalize data for plotting (0-1 scale for each measure)
    normalized_data <- combined_centrality
    for(col in colnames(normalized_data)) {
        col_data <- normalized_data[, col]
        normalized_data[, col] <- (col_data - min(col_data)) / (max(col_data) - min(col_data) + 1e-10)
    }
    
    # Create layout
    par(mfrow = c(2, 2), mar = c(12, 4, 4, 2))
    
    # Plot each measure
    measures <- c("outdeg", "indeg", "btw", "influ")
    titles <- c("Sender (Out-degree)", "Receiver (In-degree)", "Mediator (Betweenness)", "Influencer")
    
    for(i in 1:4) {
        if(measures[i] %in% colnames(normalized_data)) {
            barplot(normalized_data[, measures[i]], 
                   las = 2, 
                   main = paste("CCL -", titles[i]),
                   ylab = "Normalized Score",
                   cex.names = 0.6,
                   col = colorRampPalette(c("white", "lightblue", "darkblue"))(100)[cut(normalized_data[, measures[i]], 100)])
        }
    }
    
    dev.off()
    
    cat("✓ Saved: CCL_MANUAL_HEATMAP.pdf\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("========================================\n\n")

cat("Files in output directory:\n")
all_files <- list.files(output_dir)
for(f in all_files) {
    file_path <- file.path(output_dir, f)
    file_size <- file.info(file_path)$size
    cat("  -", f)
    if(grepl("\\.pdf$", f)) {
        if(file_size > 10000) {
            cat(" ✓✓✓ (", file_size, "bytes - HAS CONTENT)\n")
        } else {
            cat(" (", file_size, "bytes - likely empty)\n")
        }
    } else {
        cat(" (", file_size, "bytes)\n")
    }
}

cat("\n========================================\n")
cat("RECOMMENDATIONS:\n")
cat("  1. Check CCL_MANUAL_HEATMAP.pdf for the combined view\n")
cat("  2. Check CCL_individual_XX_detailed.pdf for each subcluster\n")
cat("  3. Check CCL_OPTION1_no_prefix.pdf if merge worked\n")
cat("  4. Use CCL_centrality_ALL_COMBINED.csv for custom plotting\n")
cat("========================================\n")