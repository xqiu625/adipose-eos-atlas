# =============================================================================
# Data Integration Rules
# =============================================================================

rule process_gurtner:
    """
    Process Gurtner et al. (GSE182001) reference data
    Blood and bone marrow eosinophils
    """
    input:
        gurtner_data = config["data"]["raw_data_dir"] + "/gurtner/"
    output:
        seurat = "results/07_integration/gurtner_processed.rds"
    log:
        "logs/integration/process_gurtner.log"
    script:
        "../scripts/07_integration/01_process_Gurtner.R"


rule extract_published_eos:
    """
    Extract eosinophils from published datasets
    """
    input:
        gurtner = "results/07_integration/gurtner_processed.rds"
    output:
        eos_data = "results/07_integration/published_eosinophils.rds"
    log:
        "logs/integration/extract_eos.log"
    script:
        "../scripts/07_integration/02_extract_EOS_published_data.R"


rule harmony_integration:
    """
    Harmony batch correction for dataset integration
    """
    input:
        our_data = "results/03_eosinophil/subclusters_res0.8.rds",
        reference = "results/07_integration/published_eosinophils.rds"
    output:
        integrated = "results/07_integration/harmony_integrated.rds"
    params:
        method = config["integration"]["method"]
    threads: config["resources"]["threads"]
    log:
        "logs/integration/harmony.log"
    script:
        "../scripts/07_integration/04_harmony_integration.R"


rule combined_umap:
    """
    Generate UMAP of integrated datasets
    """
    input:
        integrated = "results/07_integration/harmony_integrated.rds"
    output:
        umap = "results/07_integration/combined_umap.pdf",
        split_umap = "results/07_integration/split_by_dataset_umap.pdf"
    params:
        colors = config["colors"]["eosinophil_clusters"]
    log:
        "logs/integration/combined_umap.log"
    script:
        "../scripts/07_integration/05_combined_UMAP.R"


rule combine_datasets:
    """
    Combine multiple eosinophil datasets
    """
    input:
        datasets = expand(
            "results/07_integration/{dataset}_processed.rds",
            dataset=["gurtner", "our_data"]
        )
    output:
        combined = "results/07_integration/all_datasets_combined.rds"
    log:
        "logs/integration/combine.log"
    script:
        "../scripts/07_integration/01_combine.R"
