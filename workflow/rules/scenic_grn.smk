# =============================================================================
# SCENIC Gene Regulatory Network Rules
# =============================================================================

rule subcluster_regulon_analysis:
    """
    Analyze regulon activity across eosinophil subclusters
    """
    input:
        adata = "results/04_trajectory/scvelo_dynamics.h5ad",
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv"
    output:
        regulon_activity = "results/06_scenic/subcluster_regulon_activity.csv",
        heatmap = "results/06_scenic/regulon_heatmap.pdf"
    threads: config["scenic"]["n_workers"]
    log:
        "logs/scenic/subcluster_regulon.log"
    script:
        "../scripts/06_scenic_grn/1201_Subcluster_Regulon_Analysis.py"


rule eos_vs_mac_regulons:
    """
    Identify eosinophil vs macrophage specific regulons
    """
    input:
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv"
    output:
        zscore_matrix = "results/06_scenic/eos_mac_zscore.csv",
        unique_regulons = "results/06_scenic/cell_type_specific_regulons.csv",
        heatmap = "results/06_scenic/eos_mac_regulon_heatmap.pdf"
    log:
        "logs/scenic/eos_mac_regulons.log"
    script:
        "../scripts/06_scenic_grn/1202_Eos_vs_Mac_Subcluster_Unique_Zscore.py"


rule eosinophil_grn:
    """
    Build eosinophil gene regulatory network
    """
    input:
        regulons = "results/06_scenic/cell_type_specific_regulons.csv",
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv"
    output:
        network = "results/06_scenic/grn_eosinophil.csv",
        plot = "results/06_scenic/grn_eosinophil.pdf"
    log:
        "logs/scenic/eos_grn.log"
    script:
        "../scripts/06_scenic_grn/1203_Eos_GRN_Network.py"


rule macrophage_grn:
    """
    Build macrophage gene regulatory network
    """
    input:
        regulons = "results/06_scenic/cell_type_specific_regulons.csv",
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv"
    output:
        network = "results/06_scenic/grn_macrophage.csv",
        plot = "results/06_scenic/grn_macrophage.pdf"
    log:
        "logs/scenic/mac_grn.log"
    script:
        "../scripts/06_scenic_grn/1203_Mac_GRN_Network.py"


rule scenic_condition_comparison:
    """
    Compare regulon activity across experimental conditions
    """
    input:
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv",
        metadata = "results/04_trajectory/cell_metadata.csv"
    output:
        zscore = "results/06_scenic/condition_regulon_zscore.csv",
        heatmap = "results/06_scenic/condition_comparison_heatmap.pdf"
    log:
        "logs/scenic/condition_comparison.log"
    script:
        "../scripts/06_scenic_grn/SCENIC_Condition_Comparison_Zscore.py"


rule scenic_publication_figures:
    """
    Generate publication-ready SCENIC figures
    """
    input:
        zscore = "results/06_scenic/condition_regulon_zscore.csv",
        eos_grn = "results/06_scenic/grn_eosinophil.csv",
        mac_grn = "results/06_scenic/grn_macrophage.csv"
    output:
        figures = directory("results/06_scenic/publication_figures/")
    log:
        "logs/scenic/publication_figures.log"
    script:
        "../scripts/06_scenic_grn/SCENIC_Publication_Figures.py"


rule relma_grn:
    """
    RELMα-dependent gene regulatory network visualization
    Firework-style network plot
    """
    input:
        regulons = "results/06_scenic/cell_type_specific_regulons.csv",
        auc_matrix = config["data"]["processed_dir"] + "/scenic/auc_mtx.csv"
    output:
        network = "results/06_scenic/relma_grn.csv",
        firework_plot = "results/06_scenic/relma_firework_grn.pdf"
    log:
        "logs/scenic/relma_grn.log"
    script:
        "../scripts/06_scenic_grn/SCENIC_RELMa_Firework_GRN.py"


rule combined_eos_mac_heatmap:
    """
    Combined eosinophil and macrophage regulon heatmap
    """
    input:
        eos_regulons = "results/06_scenic/grn_eosinophil.csv",
        mac_regulons = "results/06_scenic/grn_macrophage.csv"
    output:
        heatmap = "results/06_scenic/combined_eos_mac_heatmap.pdf"
    log:
        "logs/scenic/combined_heatmap.log"
    script:
        "../scripts/06_scenic_grn/1208_Combined_Eos_Mac_Heatmap.py"


rule export_regulons:
    """
    Export final regulon list for manuscript
    """
    input:
        regulon_activity = "results/06_scenic/subcluster_regulon_activity.csv",
        unique_regulons = "results/06_scenic/cell_type_specific_regulons.csv"
    output:
        regulons = "results/06_scenic/regulons.csv"
    run:
        import pandas as pd
        activity = pd.read_csv(input.regulon_activity)
        unique = pd.read_csv(input.unique_regulons)
        combined = pd.merge(activity, unique, on='regulon', how='outer')
        combined.to_csv(output.regulons, index=False)
