# =============================================================================
# Publication Figure Generation Rules
# =============================================================================

rule figure2_eosinophil_atlas:
    """
    Figure 2: Eosinophil single-cell atlas
    - Panel B: Top 5 markers per cluster
    - Panel D: Cluster composition by condition
    - Panel G: RNA velocity UMAP
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds",
        velocity = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        figure = "results/08_figures/Figure2_eosinophil_atlas.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/figures/figure2.log"
    shell:
        """
        Rscript workflow/scripts/08_figure_generation/Fig8E_top_5markers_Eos.R && \
        Rscript workflow/scripts/08_figure_generation/Fig8D_Eos_condition_cluster_composition.R && \
        python workflow/scripts/08_figure_generation/Fig8C_UMAP_plot_by_clusters_scVelo_dynamics.py
        """


rule figure2b_markers:
    """
    Figure 2B: Top 5 marker genes per eosinophil cluster
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        plot = "results/08_figures/Fig2B_top5_markers.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/figures/fig2b.log"
    script:
        "../scripts/08_figure_generation/Fig8E_top_5markers_Eos.R"


rule figure2d_composition:
    """
    Figure 2D: Eosinophil cluster composition by condition
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        plot = "results/08_figures/Fig2D_cluster_composition.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/figures/fig2d.log"
    script:
        "../scripts/08_figure_generation/Fig8D_Eos_condition_cluster_composition.R"


rule figure2g_velocity:
    """
    Figure 2G: RNA velocity stream plot
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad"
    output:
        plot = "results/08_figures/Fig2G_velocity_umap.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/figures/fig2g.log"
    script:
        "../scripts/08_figure_generation/Fig8C_UMAP_plot_by_clusters_scVelo_dynamics.py"


rule figure3_trajectory:
    """
    Figure 3: Trajectory analysis
    - Pseudotime heatmap
    - CellRank fate probabilities
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad",
        cellrank = "results/04_trajectory/cellrank_terminal_states.csv"
    output:
        figure = "results/08_figures/Figure3_trajectory.pdf"
    log:
        "logs/figures/figure3.log"
    script:
        "../scripts/04_trajectory_analysis/04_peudotime_heatmap_1024.py"


rule figure4_cellchat:
    """
    Figure 4: Cell-cell communication
    - CellChat chord diagrams
    - CCL pathway heatmaps
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        figure = "results/08_figures/Figure4_cellchat.pdf"
    log:
        "logs/figures/figure4.log"
    script:
        "../scripts/05_cell_communication/04.1_CCL_pathway_heatmaps.R"


rule figure5_scenic:
    """
    Figure 5: Gene regulatory networks
    - Regulon activity heatmap
    - GRN visualization
    """
    input:
        regulons = "results/06_scenic/regulons.csv",
        grn_eos = "results/06_scenic/grn_eosinophil.pdf"
    output:
        figure = "results/08_figures/Figure5_scenic.pdf"
    log:
        "logs/figures/figure5.log"
    script:
        "../scripts/06_scenic_grn/SCENIC_Publication_Figures.py"


rule go_enrichment_figures:
    """
    GO enrichment analysis figures for eosinophils
    """
    input:
        deg = "results/03_eosinophil/deg_results/"
    output:
        plot = "results/08_figures/GO_enrichment_eosinophils.pdf"
    log:
        "logs/figures/go_eos.log"
    script:
        "../scripts/08_figure_generation/02_EOS_subcluster_GO_enrichment_no_gridlines.R"


rule go_enrichment_macrophages:
    """
    GO enrichment analysis figures for macrophages
    """
    input:
        deg = "results/03_eosinophil/macrophage_deg.csv"
    output:
        plot = "results/08_figures/GO_enrichment_macrophages.pdf"
    log:
        "logs/figures/go_mac.log"
    script:
        "../scripts/08_figure_generation/GO_Analysis_All_Terms_Macrophages.R"


rule deg_summary_stats:
    """
    Generate DEG summary statistics for manuscript
    """
    input:
        deg = "results/03_eosinophil/deg_results/"
    output:
        stats = "results/08_figures/DEG_summary_statistics.csv"
    log:
        "logs/figures/deg_stats.log"
    script:
        "../scripts/08_figure_generation/DEG_Summary_Statistics.R"


rule infection_deg_figure:
    """
    Infected vs non-infected DEG and GO dotplot
    """
    input:
        deg = "results/03_eosinophil/deg_results/DEG_infection.csv"
    output:
        plot = "results/08_figures/infection_DEG_GO_dotplot.pdf"
    log:
        "logs/figures/infection_deg.log"
    script:
        "../scripts/08_figure_generation/Fig9_inf_vs_noninf_DEG_GO_dotplot.R"
