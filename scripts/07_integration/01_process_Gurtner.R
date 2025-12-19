lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# Load necessary libraries
library(dplyr)
library(tidyr)
library(Seurat)

# Set the data directory
data_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/"

# List of file names
file_names <- c("GSM5515944_Colon.st.txt.gz", 
                "GSM5515944_SI.st.txt.gz", 
                "GSM5515944_Stomach.st.txt.gz", 
                "GSM5515944_Spleen.st.txt.gz", 
                "GSM5515945_blood.st.txt.gz", 
                "GSM5515945_bonemarrow.st.txt.gz",
                "GSM6675068_adipose.st.txt.gz",
                "GSM6675068_lung.st.txt.gz",
                "GSM6675068_uterus.st.txt.gz",
                "GSM5515945_stomachHP.st.txt.gz",
                "GSM5515946_bloodCR.st.txt.gz",
                "GSM5515946_bonemarrowCR.st.txt.gz",
                "GSM5515946_colonCR.st.txt.gz")

# Function to process a single file
process_file <- function(file_name) {
  file_path <- file.path(data_dir, file_name)
  
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }
  
  # Read the gzipped file
  data <- try({
    conn <- gzfile(file_path, "rt")
    d <- read.table(conn, header = TRUE, sep = "\t", row.names = NULL)
    close(conn)
    d
  })
  
  if (inherits(data, "try-error")) {
    warning(paste("Error reading file:", file_path))
    return(NULL)
  }
  
  # Reshape the data
  reshaped_data <- data %>%
    dplyr::select(Cell_Index, Gene, RSEC_Adjusted_Molecules) %>%
    tidyr::pivot_wider(names_from = Cell_Index, values_from = RSEC_Adjusted_Molecules, values_fill = list(RSEC_Adjusted_Molecules = 0))
  

  # Convert to matrix
  gene_names <- reshaped_data$Gene
  reshaped_data <- reshaped_data %>% dplyr::select(-Gene)
  reshaped_matrix <- as.matrix(reshaped_data)
  rownames(reshaped_matrix) <- gene_names
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = reshaped_matrix)
  
  # Add tissue information
  tissue <- gsub("GSM[0-9]+_|.st.txt.gz", "", file_name)
  seurat_obj$tissue <- tissue
  
  return(seurat_obj)
}

# Process all files
seurat_objects <- lapply(file_names, process_file)
names(seurat_objects) <- gsub(".st.txt.gz", "", file_names)

# Print summaries
for (name in names(seurat_objects)) {
  cat("\nSummary for", name, ":\n")
  print(seurat_objects[[name]])
}

# Merge all objects into one (if you want to analyze them together)
merged_seurat <- merge(seurat_objects[[1]], 
                       y = seurat_objects[2:length(seurat_objects)],
                       add.cell.ids = names(seurat_objects))

# Print summary of merged object
cat("\nSummary of merged Seurat object:\n")
print(merged_seurat)

# Function to add dataset-specific cell IDs, organ, and condition
add_metadata <- function(seurat_obj, dataset_name) {
  seurat_obj$dataset <- dataset_name
  seurat_obj$cell_id <- paste0(dataset_name, "_", colnames(seurat_obj))
  
  # Extract organ information
  organ <- case_when(
    grepl("Colon", dataset_name) ~ "Colon",
    grepl("SI", dataset_name) ~ "Small Intestine",
    grepl("Stomach", dataset_name) ~ "Stomach",
    grepl("Spleen", dataset_name) ~ "Spleen",
    grepl("blood", dataset_name) ~ "Blood",
    grepl("bonemarrow", dataset_name) ~ "Bone Marrow",
    grepl("adipose", dataset_name) ~ "Adipose",
    grepl("lung", dataset_name) ~ "Lung",
    grepl("uterus", dataset_name) ~ "Uterus",
    TRUE ~ "Unknown"
  )
  seurat_obj$organ <- organ
  
  # Extract condition information
  condition <- case_when(
    grepl("CR", dataset_name) ~ "C.Rodentium infection",
    grepl("HP", dataset_name) ~ "H.pylori infection",
    TRUE ~ "Steady-state"
  )
  seurat_obj$condition <- condition
  
  seurat_obj
}
# Add metadata to each Seurat object
seurat_objects_with_metadata <- lapply(names(seurat_objects), function(name) {
  add_metadata(seurat_objects[[name]], name)
})

# Merge all objects into one
merged_seurat <- merge(seurat_objects_with_metadata[[1]],
                       y = seurat_objects_with_metadata[2:length(seurat_objects_with_metadata)],
                       add.cell.ids = names(seurat_objects),
                       project = "Gurtner")

merged_seurat <- JoinLayers(merged_seurat)

# Save the merged Seurat object
saveRDS(merged_seurat, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_merged_tissues_seurat_object.rds")


# here are the rds offered from https://zenodo.org/records/10225349

colon_eosinophils <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/B6_colon_eosinophils.rds")
# 326 cells

eosinophil_challenge <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/eosinophil_challenge.rds")
# 1942 cells

eosinophils_steadystate <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/eosinophils_steadystate.rds")
# cannot load this file