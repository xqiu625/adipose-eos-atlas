library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

eos_subset_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds"
eos <- readRDS(eos_subset_file)

eos@meta.data$Clusters <- paste0("C", eos@meta.data$RNA_snn_res.0.8)
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB")
write.csv(eos@meta.data, file = "eos_sub_metadata.csv", row.names = TRUE, quote = FALSE)


library(tibble)
library(dplyr)

# Load both CSV files
eos_sub_metadata <- read.csv("eos_sub_metadata.csv", row.names = 1)  # Assuming cell names are in the first column
eos_mac_metadata <- read.csv("eos_mac_metadata.csv", row.names = 1)  # Assuming cell names are in the first column

# Ensure the row names are named consistently
rownames(eos_sub_metadata) <- make.names(rownames(eos_sub_metadata), unique = TRUE)
rownames(eos_mac_metadata) <- make.names(rownames(eos_mac_metadata), unique = TRUE)

# Merge the dataframes, keeping all rows from eos_mac_metadata
merged_metadata <- left_join(
  eos_mac_metadata %>% rownames_to_column("Cell"),
  eos_sub_metadata %>% rownames_to_column("Cell") %>% select(Cell, Clusters),
  by = "Cell"
)

# Fill NA values in the Cluster column with 'Mac'
merged_metadata$Clusters <- ifelse(is.na(merged_metadata$Clusters), 'Mac', merged_metadata$Clusters)

# Set Cell column back as row names
merged_metadata <- merged_metadata %>% column_to_rownames("Cell")

# Save the result
write.csv(merged_metadata, "eos_mac_metadata_with_eos_subclusters.csv")

# Print some information about the result
cat("Total cells:", nrow(merged_metadata), "\n")
cat("Cells with original clusters:", sum(merged_metadata$Clusters != 'Mac'), "\n")
cat("Cells labeled as Mac:", sum(merged_metadata$Clusters == 'Mac'), "\n")
cat("\nCluster distribution:\n")
print(table(merged_metadata$Clusters))


#   C0    C1    C2    C3    C4    C5    C6   Mac
# 3092  1899  1831  1338  1143  1126   269 24489
