# =============================================================================
# Cell-Cell Communication Rules
# =============================================================================

# -----------------------------------------------------------------------------
# CellChat Analysis
# -----------------------------------------------------------------------------

rule prepare_cellchat:
    """
    Prepare data for CellChat analysis
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        cellchat_input = "results/05_cellchat/cellchat_input.rds"
    log:
        "logs/cellchat/prepare.log"
    script:
        "../scripts/05_cell_communication/01_prepare_cellchat_data_by_EosSubclusters.R"


rule cellchat_by_condition:
    """
    Run CellChat analysis comparing conditions
    Eosinophil-Macrophage interactions
    """
    input:
        cellchat_input = "results/05_cellchat/cellchat_input.rds"
    output:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds",
        comparison = "results/05_cellchat/condition_comparison.pdf"
    params:
        database = config["cellchat"]["database"],
        min_cells = config["cellchat"]["min_cells"]
    threads: config["resources"]["threads"]
    log:
        "logs/cellchat/by_condition.log"
    script:
        "../scripts/05_cell_communication/02_cellchat_score_eos-mac_condition_1021.R"


rule cellchat_eos_subtypes:
    """
    CellChat analysis with eosinophil subtype resolution
    """
    input:
        cellchat_input = "results/05_cellchat/cellchat_input.rds"
    output:
        cellchat = "results/05_cellchat/cellchat_eos_subtypes.rds",
        heatmap = "results/05_cellchat/eos_subtype_interactions.pdf"
    threads: config["resources"]["threads"]
    log:
        "logs/cellchat/eos_subtypes.log"
    script:
        "../scripts/05_cell_communication/03_cellchat_score_eos-mac_eos-sub_1021.R"


rule ccl_pathway_heatmaps:
    """
    CCL chemokine pathway visualization
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        heatmaps = "results/05_cellchat/ccl_pathway_analysis.pdf"
    log:
        "logs/cellchat/ccl_heatmaps.log"
    script:
        "../scripts/05_cell_communication/04.1_CCL_pathway_heatmaps.R"


rule ccl_ligand_analysis:
    """
    Detailed CCL ligand-receptor analysis
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        analysis = "results/05_cellchat/ccl_ligand_receptor.pdf"
    log:
        "logs/cellchat/ccl_ligands.log"
    script:
        "../scripts/05_cell_communication/04.4_Ccl_selected_ligands_clean_v2_1113.R"


rule signaling_direction:
    """
    Outgoing and incoming signaling analysis
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        outgoing = "results/05_cellchat/outgoing_signaling.pdf",
        incoming = "results/05_cellchat/incoming_signaling.pdf"
    log:
        "logs/cellchat/signaling_direction.log"
    script:
        "../scripts/05_cell_communication/04.7.1-outgoing_and_incoming_signaling_condition_improved_v2_1113.R"


rule cellchat_multipanel:
    """
    Multi-panel difference heatmap visualization
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        multipanel = "results/05_cellchat/multipanel_difference_heatmap.pdf"
    log:
        "logs/cellchat/multipanel.log"
    script:
        "../scripts/05_cell_communication/1128_CellChat_MultiPanel_Difference_Heatmap.R"


rule ccl_dotplots:
    """
    CCL interaction dotplots (top 20)
    """
    input:
        cellchat = "results/05_cellchat/cellchat_eos_mac.rds"
    output:
        dotplot = "results/05_cellchat/ccl_dotplot_top20.pdf"
    log:
        "logs/cellchat/ccl_dotplots.log"
    script:
        "../scripts/05_cell_communication/1211_CCL_Dotplots_TOP20_v2.R"


# -----------------------------------------------------------------------------
# CellPhoneDB Analysis
# -----------------------------------------------------------------------------

rule prepare_cellphonedb:
    """
    Prepare inputs for CellPhoneDB
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        counts = "results/05_cellphonedb/counts.h5ad",
        meta = "results/05_cellphonedb/meta.tsv"
    log:
        "logs/cellphonedb/prepare.log"
    script:
        "../scripts/05_cell_communication/03_save_h5ad.py"


rule run_cellphonedb:
    """
    Run CellPhoneDB statistical analysis
    """
    input:
        counts = "results/05_cellphonedb/counts.h5ad",
        meta = "results/05_cellphonedb/meta.tsv"
    output:
        means = "results/05_cellphonedb/means.txt",
        pvalues = "results/05_cellphonedb/pvalues.txt",
        significant = "results/05_cellphonedb/significant_means.txt"
    params:
        iterations = config["cellphonedb"]["iterations"],
        threshold = config["cellphonedb"]["threshold"]
    threads: config["resources"]["threads"]
    log:
        "logs/cellphonedb/run.log"
    script:
        "../scripts/05_cell_communication/04_run_cellphonedb_analysis.py"


rule visualize_cellphonedb:
    """
    Visualize CellPhoneDB results
    """
    input:
        means = "results/05_cellphonedb/means.txt",
        pvalues = "results/05_cellphonedb/pvalues.txt"
    output:
        plots = "results/05_cellphonedb/interaction_plots.pdf"
    log:
        "logs/cellphonedb/visualize.log"
    script:
        "../scripts/05_cell_communication/05_visualize_inter_result.py"


rule cellphonedb_by_condition:
    """
    CellPhoneDB analysis stratified by condition
    """
    input:
        counts = "results/05_cellphonedb/counts.h5ad",
        meta = "results/05_cellphonedb/meta.tsv"
    output:
        results = directory("results/05_cellphonedb/by_condition/")
    threads: config["resources"]["threads"]
    log:
        "logs/cellphonedb/by_condition.log"
    script:
        "../scripts/05_cell_communication/06_run_cellphonedb_analysis_byCondition.py"


rule compare_interaction_strength:
    """
    Compare interaction strength across conditions
    """
    input:
        results = "results/05_cellphonedb/by_condition/"
    output:
        comparison = "results/05_cellphonedb/interaction_strength_comparison.pdf"
    log:
        "logs/cellphonedb/compare_strength.log"
    script:
        "../scripts/05_cell_communication/11_compare_interaction_strength.py"
