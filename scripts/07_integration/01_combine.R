
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

library(dplyr)
library(Seurat)
library(ggplot2)

# local Eos 
eos <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")

# GSM7919060 
data_dir <- "/bigdata/godziklab/shared/Xinru/302025/dataset/Immunity_24May/GSM7919060"
GSM7919060_seurat_object <- Read10X(data.dir = data_dir)
GSM7919060_seurat_object <- CreateSeuratObject(counts = GSM7919060_seurat_object)

# Create a list of Seurat objects
seurat_list <- list(eos, GSM7919060_seurat_object)

# Add a batch column to each object's metadata
for (i in 1:length(seurat_list)) {
    seurat_list[[i]]$batch <- paste0("Batch", i)
}

# Remove the existing dimensional reductions from each object
for (i in 1:length(seurat_list)) {
    for (dr in Reductions(seurat_list[[i]])) {
        seurat_list[[i]][[dr]] <- NULL
    }
}

# Merge the objects
combined <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], project = "integrated")
saveRDS(combined, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Immunity/eos_GSM7919060.rds")


