# =============================================================================
# Preprocessing Rules
# =============================================================================

rule seurat_preprocess:
    """
    Load 10x Genomics FLEX data and create Seurat object
    - Quality control filtering
    - Normalization and scaling
    - Variable feature selection
    """
    input:
        raw_data = config["data"]["raw_data_dir"]
    output:
        seurat = "results/01_preprocessing/seurat_processed.rds",
        qc_plots = "results/01_preprocessing/qc_metrics.pdf"
    params:
        min_cells = config["preprocessing"]["min_cells"],
        min_features = config["preprocessing"]["min_features"],
        max_features = config["preprocessing"]["max_features"],
        max_mt = config["preprocessing"]["max_mt_percent"],
        n_var_features = config["preprocessing"]["n_variable_features"]
    threads: config["resources"]["threads"]
    log:
        "logs/preprocessing/seurat_preprocess.log"
    script:
        "../scripts/01_preprocessing/01_local_data_seurat_preprocess.R"


rule run_umap:
    """
    Perform PCA and UMAP dimensionality reduction
    """
    input:
        seurat = "results/01_preprocessing/seurat_processed.rds"
    output:
        seurat = "results/01_preprocessing/seurat_umap.rds",
        umap_plot = "results/01_preprocessing/umap_all_cells.pdf"
    params:
        n_pcs = config["dimred"]["n_pcs"],
        umap_dims = config["dimred"]["umap_dims"]
    threads: config["resources"]["threads"]
    log:
        "logs/preprocessing/umap.log"
    script:
        "../scripts/01_preprocessing/03_UMAP.R"


rule find_markers:
    """
    Identify cluster marker genes using FindAllMarkers
    """
    input:
        seurat = "results/01_preprocessing/seurat_umap.rds"
    output:
        markers = "results/01_preprocessing/cluster_markers.csv",
        heatmap = "results/01_preprocessing/marker_heatmap.pdf"
    params:
        test = config["deg"]["test"],
        min_pct = config["deg"]["min_pct"],
        logfc_threshold = config["deg"]["logfc_threshold"]
    threads: config["resources"]["threads"]
    log:
        "logs/preprocessing/find_markers.log"
    script:
        "../scripts/01_preprocessing/04_FindAllMarkers_cluster.R"
