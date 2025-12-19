# =============================================================================
# Trajectory Analysis Rules
# =============================================================================

rule extract_metadata:
    """
    Extract metadata from Seurat for Python analysis
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        metadata = "results/04_trajectory/cell_metadata.csv",
        umap = "results/04_trajectory/umap_coordinates.csv"
    log:
        "logs/trajectory/extract_meta.log"
    script:
        "../scripts/04_trajectory_analysis/01_extract_meta.R"


rule scvelo_prep:
    """
    Prepare data for RNA velocity analysis
    - Convert Seurat to AnnData
    - Add spliced/unspliced counts from velocyto
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds",
        loom = config["data"]["raw_data_dir"] + "/velocyto/"
    output:
        adata = "results/04_trajectory/adata_velocity.h5ad"
    threads: config["resources"]["threads"]
    log:
        "logs/trajectory/scvelo_prep.log"
    script:
        "../scripts/04_trajectory_analysis/02_EOS_scVelo_preparation.py"


rule scvelo_dynamics:
    """
    RNA velocity analysis using dynamical model
    """
    input:
        adata = "results/04_trajectory/adata_velocity.h5ad"
    output:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad",
        stream_plot = "results/04_trajectory/velocity_stream.pdf",
        embedding_plot = "results/04_trajectory/velocity_embedding.pdf"
    params:
        mode = config["trajectory"]["scvelo_mode"],
        n_genes = config["trajectory"]["n_top_genes"]
    threads: config["resources"]["threads"]
    log:
        "logs/trajectory/scvelo_dynamics.log"
    script:
        "../scripts/04_trajectory_analysis/03_EOS_scVelo_dynamics.py"


rule scvelo_visualization:
    """
    Generate velocity stream plots with bold arrows
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        plot = "results/04_trajectory/velocity_bold_arrows.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/trajectory/scvelo_viz.log"
    script:
        "../scripts/04_trajectory_analysis/01_EOS_scvelo_bold_arrow_1024.py"


rule cellrank_analysis:
    """
    CellRank2 trajectory inference
    - Identify terminal states
    - Compute fate probabilities
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        terminal_states = "results/04_trajectory/cellrank_terminal_states.csv",
        fate_probs = "results/04_trajectory/cellrank_fate_probabilities.csv",
        plots = "results/04_trajectory/cellrank_analysis.pdf"
    threads: config["resources"]["threads"]
    log:
        "logs/trajectory/cellrank.log"
    script:
        "../scripts/04_trajectory_analysis/02_CellRank2.py"


rule pseudotime_heatmap:
    """
    Generate pseudotime heatmap of top variable genes
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        heatmap = "results/04_trajectory/pseudotime_heatmap.pdf",
        gene_order = "results/04_trajectory/pseudotime_gene_order.csv"
    params:
        n_genes = 1000,
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/trajectory/pseudotime_heatmap.log"
    script:
        "../scripts/04_trajectory_analysis/04_peudotime_heatmap_1024.py"


rule scvelo_by_condition:
    """
    RNA velocity analysis stratified by condition
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        plots = "results/04_trajectory/velocity_by_condition.pdf"
    params:
        conditions = config["samples"]["conditions"]
    log:
        "logs/trajectory/scvelo_condition.log"
    script:
        "../scripts/04_trajectory_analysis/04_EOS_scVelo_by_condition.py"
