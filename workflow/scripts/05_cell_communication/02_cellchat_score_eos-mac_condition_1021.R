library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)

# Load directories
DATA_DIR <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed"
CellChat_dir <- '/bigdata/nairlab/shared/XQiu/EOS_project/R17_20251021/CellChat'
# Create CellChat_dir if it does not exist
if (!dir.exists(CellChat_dir)) {
  dir.create(CellChat_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the Eosinophils and Macrophages subset
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_object_with_mac_subtypes_251021.rds"))
eos_mac_object <- subset(seurat_object, cell_type %in% c("Eosinophils", "Macrophages"))
eos_mac_object@meta.data$groups <- ifelse(eos_mac_object@meta.data$cell_type == "Eosinophils", "Eosinophils",eos_mac_object@meta.data$macrophage_type)


analyze_mac_eos_interactions <- function(seurat_subset) {
    # Create CellChat object and process as before
    data.input <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
    
    meta <- data.frame(
        group = case_when(
            seurat_subset$cell_type == "Eosinophils" ~ "Eosinophils",
            seurat_subset$macrophage_type == "M1_Macrophages" ~ "M1_Mac",
            seurat_subset$macrophage_type == "M2_Macrophages_Functional" ~ "M2_Functional",
            seurat_subset$macrophage_type == "M2_Macrophages_Surface_Only" ~ "M2_Surface",
            TRUE ~ NA_character_  # Handle any cells that don't match
        ),
        row.names = colnames(seurat_subset)
    )
    
    # Remove NA values if any
    meta <- meta[!is.na(meta$group), , drop = FALSE]
    
    # Drop unused factor levels
    meta$group <- droplevels(as.factor(meta$group))
    
    # Print cell type composition for this condition
    cat("Cell type composition:\n")
    print(table(meta$group))
    cat("\n")
    
    # Filter the data.input to match the filtered meta
    data.input <- data.input[, rownames(meta)]
    
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
    cellchat@DB <- CellChatDB.mouse
    
    # Process CellChat object
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # Save full CellChat object
    return(cellchat)
}

# Create condition subsets
conditions <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Female_HFD", 
               "KO_Male_HFD", "WT_Female_HFDInfected", "KO_Female_HFDInfected")

condition_subsets <- list()
for(cond in conditions) {
    condition_subsets[[cond]] <- subset(seurat_object, subset = condition == cond)
    cat("Subsetting", cond, "- Cells:", ncol(condition_subsets[[cond]]), "\n")
}

# Process each condition and save CellChat object
for(cond in conditions) {
    cat("\n========================================\n")
    cat("Processing", cond, "...\n")
    cat("========================================\n")
    cellchat_obj <- analyze_mac_eos_interactions(condition_subsets[[cond]])
    saveRDS(cellchat_obj, file.path(CellChat_dir, paste0(cond, "_cellchat.rds")))
    cat("Saved", cond, "CellChat object\n")
}




cellchat_obj <- analyze_mac_eos_interactions(eos_mac_object)
saveRDS(cellchat_obj, file.path(CellChat_dir, "eos_mac_cellchat.rds"))


R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Rocky Linux 8.8 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.15.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Los_Angeles
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  utils     datasets  grDevices methods   base     

other attached packages:
 [1] patchwork_1.3.0     CellChat_2.1.2      bigmemory_4.6.4     Biobase_2.66.0      BiocGenerics_0.52.0 ggplot2_3.5.1      
 [7] igraph_2.1.1        dplyr_1.1.4         Seurat_5.1.0        SeuratObject_5.0.2  sp_2.1-4            setwidth_1.0-4     
[13] colorout_1.3-2     

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22       splines_4.4.2          later_1.3.2            tibble_3.2.1           polyclip_1.10-7       
  [6] ggnetwork_0.5.13       fastDummies_1.7.4      lifecycle_1.0.4        rstatix_0.7.2          doParallel_1.0.17     
 [11] globals_0.16.3         lattice_0.22-6         MASS_7.3-61            backports_1.5.0        magrittr_2.0.3        
 [16] plotly_4.10.4          sass_0.4.9             jquerylib_0.1.4        httpuv_1.6.15          NMF_0.28              
 [21] sctransform_0.4.1      spam_2.11-0            spatstat.sparse_3.1-0  reticulate_1.39.0      cowplot_1.1.3         
 [26] pbapply_1.7-2          RColorBrewer_1.1-3     abind_1.4-8            Rtsne_0.17             purrr_1.0.2           
 [31] circlize_0.4.16        IRanges_2.40.0         S4Vectors_0.44.0       ggrepel_0.9.6          irlba_2.3.5.1         
 [36] listenv_0.9.1          spatstat.utils_3.1-1   goftest_1.2-3          RSpectra_0.16-2        spatstat.random_3.3-2 
 [41] fitdistrplus_1.2-1     parallelly_1.38.0      svglite_2.1.3          leiden_0.4.3.1         codetools_0.2-20      
 [46] tidyselect_1.2.1       shape_1.4.6.1          farver_2.1.2           matrixStats_1.4.1      stats4_4.4.2          
 [51] spatstat.explore_3.3-3 jsonlite_1.8.9         GetoptLong_1.0.5       BiocNeighbors_2.0.0    progressr_0.15.0      
 [56] Formula_1.2-5          ggridges_0.5.6         ggalluvial_0.12.5      survival_3.7-0         iterators_1.0.14      
 [61] systemfonts_1.1.0      foreach_1.5.2          tools_4.4.2            sna_2.8                ica_1.0-3             
 [66] Rcpp_1.0.13-1          glue_1.8.0             gridExtra_2.3          withr_3.0.2            BiocManager_1.30.25   
 [71] fastmap_1.2.0          fansi_1.0.6            digest_0.6.37          R6_2.5.1               mime_0.12             
 [76] colorspace_2.1-1       scattermore_1.2        tensor_1.5             spatstat.data_3.1-2    utf8_1.2.4            
 [81] tidyr_1.3.1            generics_0.1.3         data.table_1.16.2      FNN_1.1.4.1            httr_1.4.7            
 [86] htmlwidgets_1.6.4      uwot_0.2.2             pkgconfig_2.0.3        gtable_0.3.6           registry_0.5-1        
 [91] ComplexHeatmap_2.22.0  lmtest_0.9-40          htmltools_0.5.8.1      carData_3.0-5          dotCall64_1.2         
 [96] clue_0.3-65            scales_1.3.0           png_0.1-8              spatstat.univar_3.1-1  bigmemory.sri_0.1.8   
[101] reshape2_1.4.4         rjson_0.2.23           uuid_1.2-1             coda_0.19-4.1          statnet.common_4.10.0 
[106] nlme_3.1-166           zoo_1.8-12             cachem_1.1.0           GlobalOptions_0.1.2    stringr_1.5.1         
[111] KernSmooth_2.23-24     parallel_4.4.2         miniUI_0.1.1.1         pillar_1.9.0           grid_4.4.2            
[116] vctrs_0.6.5            RANN_2.6.2             promises_1.3.0         ggpubr_0.6.0           car_3.1-3             
[121] xtable_1.8-4           cluster_2.1.6          cli_3.6.3              compiler_4.4.2         rlang_1.1.4           
[126] crayon_1.5.3           rngtools_1.5.2         future.apply_1.11.3    ggsignif_0.6.4         plyr_1.8.9            
[131] stringi_1.8.4          network_1.18.2         viridisLite_0.4.2      deldir_2.0-4           gridBase_0.4-7        
[136] munsell_0.5.1          lazyeval_0.2.2         spatstat.geom_3.3-3    Matrix_1.7-1           RcppHNSW_0.6.0        
[141] future_1.34.0          shiny_1.9.1            ROCR_1.0-11            broom_1.0.7            bslib_0.8.0           
> 