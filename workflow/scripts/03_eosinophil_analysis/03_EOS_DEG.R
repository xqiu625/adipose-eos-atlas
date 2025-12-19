
library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)

seurat_file <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Mouse_WT_KO_HFD_Infected_merged_clustered.rds"
seurat_object <- readRDS(seurat_file)

# 1. Subset eosinophils (if not already done)
eosinophils <- subset(seurat_object, cell_type == "Eosinophils")

# 2. Set up the variables for comparison
eosinophils$Sex <- factor(eosinophils$Sex)
eosinophils$condition <- factor(eosinophils$condition)
eosinophils$infection_status <- factor(eosinophils$infection_status)
eosinophils$genotype <- factor(eosinophils$genotype)

# 3. Perform differential expression analysis for each variable
# Sex
Idents(eosinophils) <- "Sex"
deg_sex <- FindMarkers(eosinophils, ident.1 = "Male", ident.2 = "Female", test.use = "MAST")

# Infection status
Idents(eosinophils) <- "infection_status"
deg_infected <- FindMarkers(eosinophils, ident.1 = "Infected", ident.2 = "Not_Infected", test.use = "MAST")

# Genotype
Idents(eosinophils) <- "genotype"
deg_genotype <- FindMarkers(eosinophils, ident.1 = "WT", ident.2 = levels(eosinophils$genotype)[levels(eosinophils$genotype) != "WT"][1], test.use = "MAST")


# Condition (compare KO_ Female _HFD vs. KO_ Male _HFD)
# Condition (compare WT_ Female _HFD vs. WT_ Male _HFD)
# Condition (compare KO_ Female _HFD vs. WT_ Female _HFD)
# Condition (compare KO_ Male _HFD vs. WT_ Male _HFD)
# Condition comparisons
# 1. KO Female HFD vs. KO Male HFD
Idents(eosinophils) <- "condition"
deg_condition_1 <- FindMarkers(eosinophils, 
                               ident.1 = "KO_Female_HFD", 
                               ident.2 = "KO_Male_HFD", 
                               test.use = "MAST")

# 2. WT Female HFD vs. WT Male HFD
deg_condition_2 <- FindMarkers(eosinophils, 
                               ident.1 = "WT_Female_HFD", 
                               ident.2 = "WT_Male_HFD", 
                               test.use = "MAST")

# 3. KO Female HFD vs. WT Female HFD
deg_condition_3 <- FindMarkers(eosinophils, 
                               ident.1 = "KO_Female_HFD", 
                               ident.2 = "WT_Female_HFD", 
                               test.use = "MAST")

# 4. KO Male HFD vs. WT Male HFD
deg_condition_4 <- FindMarkers(eosinophils, 
                               ident.1 = "KO_Male_HFD", 
                               ident.2 = "WT_Male_HFD", 
                               test.use = "MAST")



# 4. Add gene names and sort results
add_gene_names <- function(deg_results) {
  deg_results$gene <- rownames(deg_results)
  deg_results <- deg_results %>% arrange(p_val_adj)
  return(deg_results)
}

deg_sex <- add_gene_names(deg_sex)
deg_infected <- add_gene_names(deg_infected)
deg_genotype <- add_gene_names(deg_genotype)

# 5. Save results
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data")
write.csv(deg_sex, "DEG_eosinophils_sex.csv", row.names = FALSE)
write.csv(deg_infected, "DEG_eosinophils_infected.csv", row.names = FALSE)
write.csv(deg_genotype, "DEG_eosinophils_genotype.csv", row.names = FALSE)

# Add gene names and sort results
deg_condition_1 <- add_gene_names(deg_condition_1)
deg_condition_2 <- add_gene_names(deg_condition_2)
deg_condition_3 <- add_gene_names(deg_condition_3)
deg_condition_4 <- add_gene_names(deg_condition_4)

# Save results
setwd("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data")
write.csv(deg_condition_1, "DEG_eosinophils_KO_Female_vs_KO_Male_HFD.csv", row.names = FALSE)
write.csv(deg_condition_2, "DEG_eosinophils_WT_Female_vs_WT_Male_HFD.csv", row.names = FALSE)
write.csv(deg_condition_3, "DEG_eosinophils_KO_Female_vs_WT_Female_HFD.csv", row.names = FALSE)
write.csv(deg_condition_4, "DEG_eosinophils_KO_Male_vs_WT_Male_HFD.csv", row.names = FALSE)


