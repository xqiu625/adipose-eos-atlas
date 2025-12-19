# Load required libraries
library(Seurat)
library(harmony)

# Load the datasets
eos_this_study <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")
eos_this_study@meta.data$cluster <- eos_this_study@meta.data$seurat_clusters
eos_gurtner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_eosinophil_subset.rds")
eos_posner <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Current/eos_posner_updated_geneSymbols.rds")

data_dir <- "/bigdata/godziklab/shared/Xinru/302025/dataset/Immunity_24May/GSM7919060"
GSM7919060_seurat_object <- Read10X(data.dir = data_dir)
eos_jorssen <- CreateSeuratObject(counts = GSM7919060_seurat_object, min.cells = 3, min.features = 200)

print(eos_this_study)
print(eos_gurtner) # BD Rhapsody Single-Cell Analysis System
print(eos_posner)
print(eos_jorssen)

print(head(rownames(eos_this_study[["RNA"]]), n = 10))
print(head(rownames(eos_gurtner[["RNA"]]), n = 10))
print(head(rownames(eos_posner[["RNA"]]), n = 10))
print(head(rownames(eos_jorssen[["RNA"]]), n = 10))


# List of Seurat objects
seurat_objects <- list(eos_this_study, eos_gurtner, eos_posner, eos_jorssen)
names(seurat_objects) <- c("this_study", "gurtner", "posner", "jorssen")

# Function to preprocess each dataset starting from counts
preprocess_seurat <- function(seurat_obj, project_name) {
  # Ensure we're working with raw counts
  DefaultAssay(seurat_obj) <- "RNA"
  if ("counts" %in% Layers(seurat_obj)) {
    seurat_obj <- SetAssayData(seurat_obj, slot = "data", new.data = GetAssayData(seurat_obj, slot = "counts"))
  }
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj$dataset <- project_name
  return(seurat_obj)
}

# Preprocess each dataset
seurat_objects <- lapply(names(seurat_objects), function(name) {
  preprocess_seurat(seurat_objects[[name]], name)
})

# Find common features
features <- Reduce(intersect, lapply(seurat_objects, function(x) rownames(x)))

# Subset to common features
seurat_objects <- lapply(seurat_objects, function(x) {
  subset(x, features = features)
})

# Merge datasets
merged_seurat <- merge(seurat_objects[[1]], 
                       y = seurat_objects[2:length(seurat_objects)], 
                       add.cell.ids = names(seurat_objects),
                       project = "integrated")

# Function to extract and combine count data from all layers
extract_combined_counts <- function(seurat_obj) {
  layers <- Layers(seurat_obj[["RNA"]])
  count_layers <- layers[grep("^counts", layers)]
  
  combined_counts <- lapply(count_layers, function(layer) {
    GetAssayData(seurat_obj, assay = "RNA", layer = layer)
  })
  
  # Combine the count matrices
  combined_counts <- do.call(cbind, combined_counts)
  return(combined_counts)
}

# Extract and combine count data
combined_counts <- extract_combined_counts(merged_seurat)

# Create a new assay with the combined counts
merged_seurat[["RNA_combined"]] <- CreateAssayObject(counts = combined_counts)

# Set the new assay as default and preprocess
DefaultAssay(merged_seurat) <- "RNA_combined"
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)
merged_seurat <- ScaleData(merged_seurat)

# Run PCA
merged_seurat <- RunPCA(merged_seurat, npcs = 30, verbose = FALSE)

# Run Harmony
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "orig.ident", assay.use = "RNA_combined")

# Run UMAP and clustering on Harmony-corrected embeddings
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

# Function to get cell barcodes from a Seurat object
get_cell_barcodes <- function(seurat_obj) {
  return(colnames(seurat_obj))
}

# Get cell barcodes from each original Seurat object
this_study_barcodes <- get_cell_barcodes(eos_this_study)
gurtner_barcodes <- get_cell_barcodes(eos_gurtner)
posner_barcodes <- get_cell_barcodes(eos_posner)
jorssen_barcodes <- get_cell_barcodes(eos_jorssen)

# Function to identify the source study for each cell
identify_cell_source <- function(cell_barcode) {
  if (cell_barcode %in% this_study_barcodes) return("This Study")
  if (cell_barcode %in% gurtner_barcodes) return("Gurtner")
  if (cell_barcode %in% posner_barcodes) return("Posner")
  if (cell_barcode %in% jorssen_barcodes) return("Jorssen")
  return("Unknown")
}

# Apply the function to all cells in merged_seurat
merged_seurat$source_study <- sapply(colnames(merged_seurat), identify_cell_source)

# Modify metadata
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  mutate(
    # Set organ and condition for Jorsse study
    organ = case_when(
      source_study == "Jorsse" ~ "Bone Marrow",
      source_study == "Posner" ~ "Dura Mater",
      source_study == "This Study" ~ "Adipose",
      TRUE ~ organ
    ),
    
    # Set condition for Jorsse and Posner studies, and change Steady-state to Normal
    condition = case_when(
      source_study %in% c("Jorsse", "Posner") ~ "Normal",
      condition == "Steady-state" ~ "Normal",
      TRUE ~ condition
    ),
    
    # Merge Sex and sex columns
    sex = coalesce(Sex, sex)
  ) %>%
  select(-Sex)  # Remove the lowercase 'sex' column


# Update organ information based on rownames and source_study
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  mutate(
    organ = case_when(
      grepl("stomach", rownames(merged_seurat@meta.data), ignore.case = TRUE) ~ "Stomach",
      grepl("colon", rownames(merged_seurat@meta.data), ignore.case = TRUE) ~ "Colon",
      source_study == "Jorssen" ~ "Bone Marrow",
      TRUE ~ organ
    )
  )

# Update condition information based on source_study
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  mutate(
    condition = case_when(
      source_study == "Jorssen" ~ "Normal",
      TRUE ~ condition
    )
  )

# Assign cluster value if it is NA
merged_seurat@meta.data <- merged_seurat@meta.data %>%
  mutate(
    cluster = ifelse(is.na(cluster), paste(source_study, condition, organ, sep = "_"), cluster)
  )

# Print the first few rows of the updated metadata to verify the changes
print(head(merged_seurat@meta.data))
print(table(merged_seurat$cluster))

# Print the first few rows of the updated metadata
print(head(merged_seurat@meta.data))

merged_seurat$sex <- ifelse(merged_seurat$sex == "female", "Female", merged_seurat$sex)
merged_seurat$sex <- ifelse(merged_seurat$sex == "male", "Male", merged_seurat$sex)
merged_seurat$sex <- ifelse(is.na(merged_seurat$sex), "Unknown", merged_seurat$sex)
merged_seurat$condition <- ifelse(is.na(merged_seurat$condition), "Unknown", merged_seurat$condition)
merged_seurat$organ <- ifelse(is.na(merged_seurat$organ), "Unknown", merged_seurat$organ)

# Print summary of cell sources
print(table(merged_seurat$source_study))
print(table(merged_seurat$condition))
print(table(merged_seurat$sex))
print(table(merged_seurat$organ))
print(table(merged_seurat$cluster))


# Define the save path
save_path <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/"
filename <- "harmony_integrated_eos_with_gurtner_posner_jorssen.rds"
full_path <- file.path(save_path, filename)
saveRDS(merged_seurat, file = full_path)

