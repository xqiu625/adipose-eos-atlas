lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
library(Seurat)
library(SeuratDisk)

# subset eos and mac
seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(seurat_file)
eos_mac_cells <- subset(seurat_object, cell_type %in% c("Eosinophils", "Macrophages"))

# working with another specie such as mouse, we recommend you to convert the gene ids to their corresponding orthologous.

## Download the file manually curated ortholog list from the Jackson Laboratory:
download.file("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
              "mouse_human_orthologs.txt")
## Read the ortholog file
orthologs <- read.delim("mouse_human_orthologs.txt", sep="\t", header=TRUE)
## Check the column names
print(colnames(orthologs))
## Check the first few rows
print(head(orthologs))
## Check unique values in the species column (it might not be called 'Common.Organism.Name')
print(unique(orthologs$Common.Organism.Name))
## Filter for mouse and human genes
mouse_human <- orthologs[orthologs$Common.Organism.Name %in% c("mouse, laboratory", "human"), ]
## Create separate dataframes for mouse and human genes
mouse_genes_df <- mouse_human[mouse_human$Common.Organism.Name == "mouse, laboratory", c("Symbol", "EntrezGene.ID", "DB.Class.Key")]
human_genes_df <- mouse_human[mouse_human$Common.Organism.Name == "human", c("Symbol", "EntrezGene.ID", "DB.Class.Key" )]


## Rename columns
colnames(mouse_genes_df) <- c("MGI.symbol", "mEntrezGene.ID", "DB.Class.Key")
colnames(human_genes_df) <- c("HGNC.symbol", "hEntrezGene.ID", "DB.Class.Key")

# Merge based on symbol
gene_map <- merge(mouse_genes_df, human_genes_df, by = "DB.Class.Key", all = FALSE)

# Check for any duplicate MGI symbols
print("Number of duplicate MGI symbols:")
print(sum(duplicated(gene_map$MGI.symbol)))

# If there are duplicates, let's see what they are
if(sum(duplicated(gene_map$MGI.symbol)) > 0) {
  print("Duplicate MGI symbols:")
  print(gene_map$MGI.symbol[duplicated(gene_map$MGI.symbol)])
}  
  
library(dplyr)
library(tibble)

# Step 1: Extract the counts as a matrix
counts_matrix <- GetAssayData(eos_mac_cells, layer = "counts")

counts_tibble <- counts_matrix %>%
  as.matrix() %>%
  as_tibble(rownames = "MGI.symbol")
  

# Step 2: Merge with gene_map
merged_counts <- counts_tibble %>%
  left_join(gene_map, by = "MGI.symbol")
  
# Step 3: Replace MGI.symbol with HGNC.symbol
human_counts <- merged_counts %>%
  select(-DB.Class.Key, -MGI.symbol, -mEntrezGene.ID, -hEntrezGene.ID) %>%
  filter(!is.na(HGNC.symbol) & HGNC.symbol != "")
  
# Check for any duplicate HGNC symbols
duplicate_count <- sum(duplicated(human_counts$HGNC.symbol))
print(paste("Number of duplicate HGNC symbols:", duplicate_count))

  
# Step 4: Convert back to matrix with HGNC.symbol as rownames
# Sum the counts for duplicate HGNC symbols
human_counts_deduplicated <- human_counts %>%
  group_by(HGNC.symbol) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup()
  
# create the matrix 
human_counts_matrix <- human_counts_deduplicated %>%
  column_to_rownames("HGNC.symbol") %>%
  as.matrix()

# Check the dimensions of the new matrix
print("Dimensions of human_counts_matrix:")
print(dim(human_counts_matrix))

# Print the first few row names
print("First few row names (HGNC symbols):")
print(head(rownames(human_counts_matrix)))


# Extract the count matrix and save as CSV
human_counts_matrix <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/human_counts_matrix.rds")
write.csv(human_counts_matrix, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/hgnc_eos_mac_counts.csv")


# Save metadata
seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(seurat_file)
eos_mac_cells <- subset(seurat_object, cell_type %in% c("Eosinophils", "Macrophages"))
metadata <- eos_mac_cells@meta.data
write.csv(metadata, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/CellphoneDB/eos_mac_metadata.csv")

