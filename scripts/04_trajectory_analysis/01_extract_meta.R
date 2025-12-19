
library(Seurat)
library(magrittr)


df <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/eos_subset_cluster_res0.8.rds")
head(Embeddings(df, reduction = "umap"))

umap <- Embeddings(df, reduction = "umap")
write.csv(umap, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_embeddings.csv")

write.csv(Cells(df), file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_cellID_obs.csv", row.names = FALSE)

clusters_df <- data.frame(CellID = row.names(df@meta.data), SeuratClusters = df@meta.data$seurat_clusters)
write.csv(clusters_df, file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/EOS_clusters.csv", row.names = FALSE)

