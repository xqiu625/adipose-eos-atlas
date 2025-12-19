# =============================================================================
# Eosinophil Analysis Rules
# =============================================================================

rule eosinophil_subcluster:
    """
    Subcluster eosinophils to identify C0-C6 subtypes
    Resolution 0.8 yields 7 distinct clusters
    """
    input:
        seurat = "results/02_annotation/seurat_annotated.rds"
    output:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds",
        umap = "results/03_eosinophil/eos_subcluster_umap.pdf"
    params:
        resolution = config["eosinophil"]["resolution"],
        colors = config["colors"]["eosinophil_clusters"]
    threads: config["resources"]["threads"]
    log:
        "logs/eosinophil/subcluster.log"
    script:
        "../scripts/03_eosinophil_analysis/07_EOS_subcluster_0.8.R"


rule eosinophil_deg:
    """
    Differential expression analysis across conditions
    Comparisons: Sex, Genotype, Infection status
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        deg_dir = directory("results/03_eosinophil/deg_results/"),
        sex_deg = "results/03_eosinophil/deg_results/DEG_sex.csv",
        genotype_deg = "results/03_eosinophil/deg_results/DEG_genotype.csv",
        infection_deg = "results/03_eosinophil/deg_results/DEG_infection.csv"
    params:
        test = config["deg"]["test"],
        logfc = config["deg"]["logfc_threshold"],
        padj = config["deg"]["padj_threshold"]
    threads: config["resources"]["threads"]
    log:
        "logs/eosinophil/deg.log"
    script:
        "../scripts/03_eosinophil_analysis/03_EOS_DEG.R"


rule eosinophil_pathway:
    """
    GO/KEGG pathway enrichment using clusterProfiler
    """
    input:
        deg = "results/03_eosinophil/deg_results/"
    output:
        go_results = "results/03_eosinophil/pathway/GO_enrichment.csv",
        kegg_results = "results/03_eosinophil/pathway/KEGG_enrichment.csv",
        plots = "results/03_eosinophil/pathway/enrichment_plots.pdf"
    threads: 4
    log:
        "logs/eosinophil/pathway.log"
    script:
        "../scripts/03_eosinophil_analysis/04_EOS_ClusterProfiler.R"


rule gsva_analysis:
    """
    Gene Set Variation Analysis for pathway scoring
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        scores = "results/03_eosinophil/gsva_scores.csv",
        heatmap = "results/03_eosinophil/gsva_heatmap.pdf"
    threads: config["resources"]["threads"]
    log:
        "logs/eosinophil/gsva.log"
    script:
        "../scripts/03_eosinophil_analysis/08_GSVA.R"


rule define_m1_m2_signatures:
    """
    Define M1/M2 macrophage polarization signatures
    """
    input:
        seurat = "results/02_annotation/seurat_annotated.rds"
    output:
        signatures = "results/03_eosinophil/m1_m2_signatures.rds"
    params:
        m1_markers = config["macrophage"]["m1_markers"],
        m2_surface = config["macrophage"]["m2_surface_markers"],
        m2_functional = config["macrophage"]["m2_functional_markers"]
    log:
        "logs/eosinophil/m1_m2_signatures.log"
    script:
        "../scripts/03_eosinophil_analysis/01_define_M1_M2_signatures.R"


rule score_m1_m2:
    """
    Score macrophages for M1/M2 polarization
    """
    input:
        seurat = "results/02_annotation/seurat_annotated.rds",
        signatures = "results/03_eosinophil/m1_m2_signatures.rds"
    output:
        scores = "results/03_eosinophil/m1_m2_scores.csv",
        plots = "results/03_eosinophil/m1_m2_polarization.pdf"
    threads: 4
    log:
        "logs/eosinophil/score_m1_m2.log"
    script:
        "../scripts/03_eosinophil_analysis/02_score_cells_M1_M2.R"


rule macrophage_deg:
    """
    Macrophage cluster-specific differential expression
    """
    input:
        seurat = "results/02_annotation/seurat_annotated.rds"
    output:
        deg = "results/03_eosinophil/macrophage_deg.csv"
    params:
        test = config["deg"]["test"]
    threads: config["resources"]["threads"]
    log:
        "logs/eosinophil/macrophage_deg.log"
    script:
        "../scripts/03_eosinophil_analysis/01_Macrophage_Cluster_DEG.R"


rule pseudobulk_analysis:
    """
    Pseudobulk analysis by sample for robust statistics
    """
    input:
        seurat = "results/03_eosinophil/subclusters_res0.8.rds"
    output:
        pseudobulk = "results/03_eosinophil/pseudobulk_counts.csv",
        de_results = "results/03_eosinophil/pseudobulk_de.csv"
    threads: 4
    log:
        "logs/eosinophil/pseudobulk.log"
    script:
        "../scripts/03_eosinophil_analysis/01_pseudobulk_by_samplev3.R"
