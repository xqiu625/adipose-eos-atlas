library(Seurat)
library(ggplot2)

# Define paths to the filtered feature barcode matrices
data_paths <- c(
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/1/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/2/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/3/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/4/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0508/outs/per_sample_outs/5/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0508/outs/per_sample_outs/6/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0508/outs/per_sample_outs/7/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0508/outs/per_sample_outs/8/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0912/outs/per_sample_outs/9/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0912/outs/per_sample_outs/10/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0912/outs/per_sample_outs/11/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0912/outs/per_sample_outs/12/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR1316/outs/per_sample_outs/13/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR1316/outs/per_sample_outs/14/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR1316/outs/per_sample_outs/15/count/sample_filtered_feature_bc_matrix/",
    "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR1316/outs/per_sample_outs/16/count/sample_filtered_feature_bc_matrix/"
)

cell_identifiers <- c(
  "WT_Male_HFD_1",
  "WT_Female_HFD_1",
  "KO_Male_HFD_1",
  "KO_Female_HFD_1",
  "WT_Female_HFDInfected_1",
  "WT_Female_HFDInfected_2",
  "KO_Female_HFDInfected_1",
  "KO_Female_HFDInfected_2",
  "WT_Male_HFD_2",
  "WT_Male_HFD_3",
  "WT_Female_HFD_2",
  "WT_Female_HFD_3",
  "KO_Male_HFD_2",
  "KO_Male_HFD_3",
  "KO_Female_HFD_2",
  "KO_Female_HFD_3"
)


# Create a list of Seurat objects with cell identifiers
seurat_objects <- mapply(function(path, id) {
  # Read the data
  data <- Read10X(data.dir = path)
  # Create a Seurat object and add cell identifier
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj <- RenameCells(seurat_obj, new.names = paste(id, Cells(seurat_obj), sep = "_"))
  return(seurat_obj)
}, data_paths, cell_identifiers, SIMPLIFY = FALSE)

# Merge Seurat objects
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)
merged_seurat@meta.data$orig.ident <- sapply(strsplit(rownames(merged_seurat@meta.data), "_"), function(x) paste(x[1:4], collapse="_"))

# 1. Normalize the data
merged_seurat <- NormalizeData(merged_seurat)

# 2. Find variable features (optional, but recommended)
merged_seurat <- FindVariableFeatures(merged_seurat)

# 3. Scale the data
merged_seurat <- ScaleData(merged_seurat)

# 4. Run PCA
merged_seurat <- RunPCA(merged_seurat)

# 5. Run UMAP
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

# Plot UMAP and color by orig.ident
umap_plot <- DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident")

# Save UMAP plot
umap_output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/figures/UMAP_plot_by_orig_ident_merged.png"
ggsave(umap_output_file, plot = umap_plot, width = 10, height = 8)
                                             
# Define the output file path
output_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged.rds"

# Save the merged_seurat object to the specified file
saveRDS(merged_seurat, file = output_file)

# No batch effect found, so did not move batch effect.

