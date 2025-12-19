# =============================================================================
# Cell Annotation Rules
# =============================================================================

rule uce_annotation:
    """
    Annotate cells using Universal Cell Embedding (UCE)
    Reference: Mouse kidney normalized embeddings
    """
    input:
        seurat = "results/01_preprocessing/seurat_umap.rds"
    output:
        embeddings = "results/02_annotation/uce_embeddings.csv",
        predictions = "results/02_annotation/uce_predictions.csv"
    params:
        reference = config["annotation"]["reference"],
        threshold = config["annotation"]["confidence_threshold"]
    threads: config["resources"]["threads"]
    log:
        "logs/annotation/uce_annotation.log"
    script:
        "../scripts/02_cell_annotation/05_UCE_annotate_by_mouse_kidney_norn.py"


rule label_cell_types:
    """
    Assign final cell type labels based on UCE predictions
    Identifies: Eosinophils, Macrophages, and other immune cells
    """
    input:
        seurat = "results/01_preprocessing/seurat_umap.rds",
        predictions = "results/02_annotation/uce_predictions.csv"
    output:
        seurat = "results/02_annotation/seurat_annotated.rds",
        assignments = "results/02_annotation/cell_type_assignments.csv",
        umap = "results/02_annotation/umap_annotated.pdf"
    threads: 4
    log:
        "logs/annotation/label_cell_types.log"
    script:
        "../scripts/02_cell_annotation/06_label_UCE_identity.R"
