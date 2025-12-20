library(CellChat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(pheatmap)

# Set directories
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat'
output_dir <- file.path(CellChat_dir, 'CCL_selected_ligands_heatmaps')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading CellChat objects...\n")

eos_clusters <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6")

eos_mac_subtypes <- list()
for(cluster in eos_clusters) {
    file_path <- file.path(CellChat_dir, paste0("eos_", cluster, "_vs_MacSubtypes_cellchat.rds"))
    if(file.exists(file_path)) {
        eos_mac_subtypes[[cluster]] <- readRDS(file_path)
        cat("  Loaded:", cluster, "\n")
    } else {
        cat("  Missing:", cluster, "\n")
    }
}

cat("\nTotal loaded:", length(eos_mac_subtypes), "objects\n\n")

cat("Extracting CCL ligand-receptor interactions...\n\n")

all_interactions <- data.frame()

for(cluster in names(eos_mac_subtypes)) {
    cat("Processing", cluster, "...\n")

    obj <- eos_mac_subtypes[[cluster]]
    lr_db <- obj@LR$LRsig

    if(!is.null(obj@net$prob)) {
        prob_matrix <- obj@net$prob
        pathways <- obj@netP$pathways
        ccl_pathways <- pathways[grepl("^CCL", pathways, ignore.case = TRUE)]

        if(length(ccl_pathways) > 0) {
            cat("  CCL pathways:", paste(ccl_pathways, collapse=", "), "\n")

            for(pathway in ccl_pathways) {
                lr_pairs <- lr_db[lr_db$pathway_name == pathway, ]

                if(nrow(lr_pairs) > 0) {
                    cat("    Pathway", pathway, ":", nrow(lr_pairs), "L-R pairs\n")

                    for(i in 1:nrow(lr_pairs)) {
                        interaction_name <- lr_pairs$interaction_name[i]
                        ligand <- lr_pairs$ligand[i]
                        receptor <- lr_pairs$receptor[i]

                        if(interaction_name %in% dimnames(prob_matrix)[[3]]) {
                            prob_mat <- prob_matrix[,,interaction_name]
                            cell_types <- rownames(prob_mat)

                            for(sender in cell_types) {
                                for(receiver in cell_types) {
                                    prob_value <- prob_mat[sender, receiver]

                                    if(!is.na(prob_value) && prob_value > 0) {
                                        all_interactions <- rbind(all_interactions, data.frame(
                                            eos_cluster = cluster,
                                            pathway = pathway,
                                            ligand = ligand,
                                            receptor = receptor,
                                            interaction = interaction_name,
                                            sender = sender,
                                            receiver = receiver,
                                            probability = prob_value,
                                            stringsAsFactors = FALSE
                                        ))
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            cat("  No CCL pathways found\n")
        }
    }
    cat("\n")
}

cat("Total interactions extracted:", nrow(all_interactions), "\n\n")

cat("APPLYING FILTERS\n")
cat("========================================\n\n")

# UPDATED: Removed Ccl8, reordered to CCL6, CCL9, CCL2, CCL7, CCL24
selected_ligands <- c("Ccl6", "Ccl9", "Ccl2", "Ccl7", "Ccl24")

cat("Selected ligands:", paste(selected_ligands, collapse=", "), "\n")
cat("Original total interactions:", nrow(all_interactions), "\n")

all_interactions_filtered <- all_interactions %>%
    filter(!grepl("other", sender, ignore.case = TRUE)) %>%
    filter(!grepl("other", receiver, ignore.case = TRUE))

cat("After removing 'other Mac':", nrow(all_interactions_filtered), "\n")

all_interactions_filtered <- all_interactions_filtered %>%
    filter(ligand %in% selected_ligands)

cat("After filtering for selected ligands:", nrow(all_interactions_filtered), "\n")

present_ligands <- unique(all_interactions_filtered$ligand)
cat("Ligands present in data:", paste(present_ligands, collapse=", "), "\n\n")

write.csv(all_interactions_filtered,
         file.path(output_dir, "CCL_selected_ligands_interactions.csv"),
         row.names = FALSE)
cat("Saved filtered interaction data\n\n")

cat("========================================\n")
cat("HEATMAP 1: Eos to Mac\n")
cat("========================================\n\n")

eos_to_mac <- all_interactions_filtered %>%
    filter(grepl("^Eos_", sender) & !grepl("^Eos_", receiver))

if(nrow(eos_to_mac) > 0) {
    cat("Found", nrow(eos_to_mac), "Eos to Mac interactions\n")
    cat("Ligands present:", paste(unique(eos_to_mac$ligand), collapse=", "), "\n")

    heatmap_data_1 <- eos_to_mac %>%
        group_by(ligand, sender, receiver) %>%
        summarise(avg_prob = mean(probability), .groups = 'drop')

    heatmap_data_1$pair <- paste(heatmap_data_1$sender, "->", heatmap_data_1$receiver)

    matrix_1 <- dcast(heatmap_data_1, ligand ~ pair, value.var = "avg_prob", fill = 0)
    rownames(matrix_1) <- matrix_1$ligand
    matrix_1 <- matrix_1[, -1]
    matrix_1 <- as.matrix(matrix_1)

    existing_ligands <- intersect(selected_ligands, rownames(matrix_1))
    if(length(existing_ligands) > 0) {
        matrix_1 <- matrix_1[existing_ligands, , drop=FALSE]
    }

    cat("Heatmap dimensions:", nrow(matrix_1), "ligands x", ncol(matrix_1), "cell pairs\n")

    pdf(file.path(output_dir, "Heatmap1_Eos_to_Mac_selected_CCL.pdf"),
        width = max(12, ncol(matrix_1) * 0.3),
        height = max(6, nrow(matrix_1) * 0.5))

    pheatmap(matrix_1,
             color = colorRampPalette(c("white", "yellow", "red", "darkred"))(100),
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             scale = "none",
             main = "Selected CCL Ligands: Eosinophil Subclusters to Macrophage Subtypes\n(Excluding other Mac)",
             fontsize = 10,
             fontsize_row = 10,
             fontsize_col = 8,
             angle_col = 45,
             cellwidth = 15,
             cellheight = 18)

    dev.off()

    cat("Saved: Heatmap1_Eos_to_Mac_selected_CCL.pdf\n")

    write.csv(matrix_1,
             file.path(output_dir, "Matrix1_Eos_to_Mac_selected_CCL.csv"),
             row.names = TRUE)
    cat("Saved: Matrix1_Eos_to_Mac_selected_CCL.csv\n\n")
} else {
    cat("WARNING: No Eos to Mac interactions found with selected filters\n\n")
}

cat("========================================\n")
cat("HEATMAP 2: Mac to Eos\n")
cat("========================================\n\n")

mac_to_eos <- all_interactions_filtered %>%
    filter(!grepl("^Eos_", sender) & grepl("^Eos_", receiver))

if(nrow(mac_to_eos) > 0) {
    cat("Found", nrow(mac_to_eos), "Mac to Eos interactions\n")
    cat("Ligands present:", paste(unique(mac_to_eos$ligand), collapse=", "), "\n")

    heatmap_data_2 <- mac_to_eos %>%
        group_by(ligand, sender, receiver) %>%
        summarise(avg_prob = mean(probability), .groups = 'drop')

    heatmap_data_2$pair <- paste(heatmap_data_2$sender, "->", heatmap_data_2$receiver)

    matrix_2 <- dcast(heatmap_data_2, ligand ~ pair, value.var = "avg_prob", fill = 0)
    rownames(matrix_2) <- matrix_2$ligand
    matrix_2 <- matrix_2[, -1]
    matrix_2 <- as.matrix(matrix_2)

    existing_ligands <- intersect(selected_ligands, rownames(matrix_2))
    if(length(existing_ligands) > 0) {
        matrix_2 <- matrix_2[existing_ligands, , drop=FALSE]
    }

    cat("Heatmap dimensions:", nrow(matrix_2), "ligands x", ncol(matrix_2), "cell pairs\n")

    pdf(file.path(output_dir, "Heatmap2_Mac_to_Eos_selected_CCL.pdf"),
        width = max(12, ncol(matrix_2) * 0.3),
        height = max(6, nrow(matrix_2) * 0.5))

    pheatmap(matrix_2,
             color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             scale = "none",
             main = "Selected CCL Ligands: Macrophage Subtypes to Eosinophil Subclusters\n(Excluding other Mac)",
             fontsize = 10,
             fontsize_row = 10,
             fontsize_col = 8,
             angle_col = 45,
             cellwidth = 15,
             cellheight = 18)

    dev.off()

    cat("Saved: Heatmap2_Mac_to_Eos_selected_CCL.pdf\n")

    write.csv(matrix_2,
             file.path(output_dir, "Matrix2_Mac_to_Eos_selected_CCL.csv"),
             row.names = TRUE)
    cat("Saved: Matrix2_Mac_to_Eos_selected_CCL.csv\n\n")
} else {
    cat("WARNING: No Mac to Eos interactions found with selected filters\n\n")
}

cat("========================================\n")
cat("HEATMAP 3: CCL Ligands Senders\n")
cat("========================================\n\n")

ligand_cell_data <- all_interactions_filtered %>%
    group_by(ligand, sender, receiver) %>%
    summarise(avg_prob = mean(probability), .groups = 'drop')

sender_matrix <- ligand_cell_data %>%
    group_by(ligand, sender) %>%
    summarise(total_send = sum(avg_prob), .groups = 'drop') %>%
    dcast(ligand ~ sender, value.var = "total_send", fill = 0)

rownames(sender_matrix) <- sender_matrix$ligand
sender_matrix <- sender_matrix[, -1]
sender_matrix <- as.matrix(sender_matrix)

existing_ligands <- intersect(selected_ligands, rownames(sender_matrix))
if(length(existing_ligands) > 0) {
    sender_matrix <- sender_matrix[existing_ligands, , drop=FALSE]
}

cat("Sender heatmap dimensions:", nrow(sender_matrix), "ligands x", ncol(sender_matrix), "cell types\n")

if(ncol(sender_matrix) > 0 && nrow(sender_matrix) > 0) {
    pdf(file.path(output_dir, "CCL_ligands_senders_heatmap_filtered.pdf"),
        width = max(10, ncol(sender_matrix) * 0.8),
        height = max(6, nrow(sender_matrix) * 0.5))

    pheatmap(sender_matrix,
             color = colorRampPalette(c("white", "pink", "red", "darkred"))(100),
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             scale = "none",
             main = "CCL Ligand Senders: Eos Subclusters and Mac Subtypes\n(Selected ligands, excluding other Mac)",
             fontsize = 10,
             fontsize_row = 10,
             fontsize_col = 9,
             angle_col = 45,
             cellwidth = 20,
             cellheight = 20)

    dev.off()
    cat("Saved: CCL_ligands_senders_heatmap_filtered.pdf\n")

    write.csv(sender_matrix,
             file.path(output_dir, "CCL_ligands_senders_matrix_filtered.csv"))
    cat("Saved: CCL_ligands_senders_matrix_filtered.csv\n\n")
}

cat("========================================\n")
cat("HEATMAP 4: CCL Ligands Receivers\n")
cat("========================================\n\n")

receiver_matrix <- ligand_cell_data %>%
    group_by(ligand, receiver) %>%
    summarise(total_receive = sum(avg_prob), .groups = 'drop') %>%
    dcast(ligand ~ receiver, value.var = "total_receive", fill = 0)

rownames(receiver_matrix) <- receiver_matrix$ligand
receiver_matrix <- receiver_matrix[, -1]
receiver_matrix <- as.matrix(receiver_matrix)

existing_ligands <- intersect(selected_ligands, rownames(receiver_matrix))
if(length(existing_ligands) > 0) {
    receiver_matrix <- receiver_matrix[existing_ligands, , drop=FALSE]
}

cat("Receiver heatmap dimensions:", nrow(receiver_matrix), "ligands x", ncol(receiver_matrix), "cell types\n")

if(ncol(receiver_matrix) > 0 && nrow(receiver_matrix) > 0) {
    pdf(file.path(output_dir, "CCL_ligands_receivers_heatmap_filtered.pdf"),
        width = max(10, ncol(receiver_matrix) * 0.8),
        height = max(6, nrow(receiver_matrix) * 0.5))

    pheatmap(receiver_matrix,
             color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             scale = "none",
             main = "CCL Ligand Receivers: Eos Subclusters and Mac Subtypes\n(Selected ligands, excluding other Mac)",
             fontsize = 10,
             fontsize_row = 10,
             fontsize_col = 9,
             angle_col = 45,
             cellwidth = 20,
             cellheight = 20)

    dev.off()
    cat("Saved: CCL_ligands_receivers_heatmap_filtered.pdf\n")

    write.csv(receiver_matrix,
             file.path(output_dir, "CCL_ligands_receivers_matrix_filtered.csv"))
    cat("Saved: CCL_ligands_receivers_matrix_filtered.csv\n\n")
}

cat("========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n\n")

ligand_summary <- all_interactions_filtered %>%
    group_by(ligand) %>%
    summarise(
        total_interactions = n(),
        total_probability = sum(probability),
        avg_probability = mean(probability),
        eos_clusters_involved = n_distinct(eos_cluster),
        unique_senders = n_distinct(sender),
        unique_receivers = n_distinct(receiver),
        .groups = 'drop'
    ) %>%
    arrange(match(ligand, selected_ligands))

write.csv(ligand_summary,
         file.path(output_dir, "Summary_selected_CCL_ligands.csv"),
         row.names = FALSE)

cat("Summary of selected CCL ligands:\n")
print(ligand_summary)
cat("\n")

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("Output directory:\n")
cat("  ", output_dir, "\n\n")

cat("Filters applied:\n")
cat("  - Excluded: other Mac cell types\n")
cat("  - Selected ligands: CCL6, CCL9, CCL2, CCL7, CCL24 (CCL8 removed)\n")
cat("  - Ligand order: CCL6 → CCL9 → CCL2 → CCL7 → CCL24\n\n")

cat("Generated files:\n")
cat("  HEATMAPS:\n")
cat("    1. Heatmap1_Eos_to_Mac_selected_CCL.pdf\n")
cat("    2. Heatmap2_Mac_to_Eos_selected_CCL.pdf\n")
cat("    3. CCL_ligands_senders_heatmap_filtered.pdf\n")
cat("    4. CCL_ligands_receivers_heatmap_filtered.pdf\n\n")

cat("  DATA MATRICES:\n")
cat("    - Matrix1_Eos_to_Mac_selected_CCL.csv\n")
cat("    - Matrix2_Mac_to_Eos_selected_CCL.csv\n")
cat("    - CCL_ligands_senders_matrix_filtered.csv\n")
cat("    - CCL_ligands_receivers_matrix_filtered.csv\n\n")

cat("  RAW DATA:\n")
cat("    - CCL_selected_ligands_interactions.csv\n")
cat("    - Summary_selected_CCL_ligands.csv\n\n")

cat("========================================\n")
cat("All 4 heatmaps show only the 5 selected CCL ligands\n")
cat("(CCL6, CCL9, CCL2, CCL7, CCL24) in the specified order\n")
cat("and exclude all other Mac interactions!\n")
cat("========================================\n")
