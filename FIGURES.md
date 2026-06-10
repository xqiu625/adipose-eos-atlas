# Figure → Script Manifest

This file maps each published figure panel to the script(s) that produce it. Scripts named `NN_figX_*` are figure-named ground truth; other panels are mapped to the analysis script that generates the underlying plot. **All paths inside the scripts are HPCC-specific** (hard-coded to the UCR `/bigdata/nairlab/...` layout) and must be edited by hand to run elsewhere. Experimental (wet-lab / imaging) panels are noted as not produced by code in this repo.

> **Panel numbering note.** Script-internal panel labels and the result-PDF filenames in `adipose-eos-atlas-results/` follow an **older submission numbering** that is offset/swapped relative to the final (0508) manuscript. Mappings below were confirmed by rendering each candidate PDF and matching its content to the final figure legend; where the script's own label differs from the manuscript panel, that is called out inline.

## Figure 1 — Atlas overview and eosinophil DEGs

| Panel | Description | Script(s) |
|-------|-------------|-----------|
| 1A | Experimental schematic | Experimental / not computational |
| 1B | All-cells UMAP, 7 immune cell types | `02_all_cells_umap.R`, `03_annotate_celltypes.R` |
| 1C | Dot plot of canonical markers across immune cell types | `25_fig1_revision_v3.R` (`DotPlot`, 7 immune cell types × canonical markers, C6 removed → `fig1D_dotplot_markers.pdf`; the script's internal "Fig 1D" label is manuscript **1C**) |
| 1D | Pseudobulk DEG volcano — WT male vs WT female eos | `05_eos_pseudobulk.R` / `07c_deg_pseudobulk.R` (DESeq2), `07_replot_volcanoes.R`, `25_fig1_revision_v3.R` |
| 1E | Pseudobulk DEG volcano — WT vs KO male eos | `05_eos_pseudobulk.R` / `07c_deg_pseudobulk.R`, `07_replot_volcanoes.R`, `25_fig1_revision_v3.R` |
| 1F | Pseudobulk DEG volcano — WT vs KO female eos | `05_eos_pseudobulk.R` / `07c_deg_pseudobulk.R`, `07_replot_volcanoes.R`, `25_fig1_revision_v3.R` |
| 1G | Pseudobulk DEG volcano — naïve vs helminth-infected female eos (covariate-controlled) | `03c_deg_infection_split.R`, `07c_deg_pseudobulk.R`, `07_replot_volcanoes.R`, `25_fig1_revision_v3.R` |
| 1D–G bottom | GSEA GO BP per comparison (butterfly plots) | `25_fig1_revision_v3.R` (computes GSEA), `26_gsea_custom_butterfly.R` (plots) |

## Figure 2 — Eosinophil heterogeneity and trajectory

| Panel | Description | Script(s) |
|-------|-------------|-----------|
| 2A | Eosinophil UMAP (7,461 cells, six clusters C0–C5) | `04_eos_subcluster.R`, `04b_rename_clusters.R` |
| 2B | Dot plot of top-6 signature genes per cluster | `31_fig2B_dotplot.R` |
| 2C | Hierarchical GO BP enrichment heatmap of cluster markers | `08_pathway_eos_clusters.R` |
| 2D | Stacked bar of eos subcluster proportions by group | `20_compositional_analysis.R` (sccomp), `43_composition_biorep.R` (per-biorep supplement) |
| 2E | Spectral flow cytometry UMAP | Experimental / not computational |
| 2F | MHC2 MFI | Experimental / not computational |
| 2G | Monocle 3 trajectory of eos clusters integrated with Gurtner BM+blood, rooted at BM | `24_integrate_qiu_hu_gurtner.R` (integration), `11b_monocle3_integrated.R` (Monocle 3 trajectory) |
| 2H | Pseudotime distributions by group (boxplots) | `11b_monocle3_integrated.R` (main-figure boxplot, groups = Gurtner BM + Gurtner Blood + C0–C5, root = Gurtner BM → `monocle3_main_pseudotime_boxplot.pdf`). Not `24c_pseudotime_violin.py` — that renders violins, not the boxplot panel. |
| 2I | Cross-dataset transcriptional similarity dendrogram vs Hu depots + Gurtner | `30_cross_dataset_similarity.R` |
| 2J | Gene expression dynamics along pseudotime heatmap | `32_fig2J_pseudotime_heatmap.R` |

## Figure 3 — Eosinophil–macrophage communication (CellChat)

| Panel | Description | Script(s) |
|-------|-------------|-----------|
| 3A | Chord diagram of signaling between eos and M1 / M2-functional / M2-surface macrophages (4 cell-type sectors: gray / pink / green / purple) | **`44_fig3A_chord_revised.R` (canonical)** — see note below |
| 3B | Outgoing / incoming signaling pathway strength heatmaps across 6 conditions | `10b_cellchat_manuscript_figures.R` |
| 3C | CCL signaling network centrality heatmap (sender/receiver/mediator/influencer) for eos / C3 / mac clusters | `10b_cellchat_manuscript_figures.R` |
| 3D | CCL ligand sender + receiver heatmaps | `10c_cellchat_fig3d.R` |
| 3E | CCL signaling Eos → Mac heatmap | `33_fig3EF_cellchat_ccl_heatmap.R` |
| 3F | CCL signaling Mac → Eos heatmap | `33_fig3EF_cellchat_ccl_heatmap.R` |

**Fig 3A chord — canonical vs superseded (resolved).** Five scripts render a Fig 3A chord diagram. `44_fig3A_chord_revised.R` is canonical: it is the latest committed version (2026-04-25, newest of all five), it is cited for Fig 3A in the README pipeline table, and it produces an **all-pathways, cell-type-sector** chord — consistent with the final repo decision (commit `9fcb506` "Fig3A: all pathways aggregated, not just CCL") and with the manuscript legend, which describes four cell-type sectors colored gray/pink/green/purple. The rendered manuscript-style panel (4 sectors Eosinophils=gray, M1_Mac=pink, M2_Functional=green, M2_Surface=purple) was confirmed by inspection. The candidate that labels **individual L-R pairs as sectors** (`36`) does **not** match the panel and is superseded. Provenance of the alternates:

- `10b_cellchat_manuscript_figures.R` — original Fig 3A chord; `9fcb506` updated it from CCL-only to all-pathways aggregated (the CCL-only render survives in `adipose-eos-atlas-results/figure3_revision/fig3A_ccl_chord_diagram.pdf`, an older version of this panel).
- `35_fig3A_chord_revised.R` — all-pathways chord with enhanced arrow visibility (R1 comment 18).
- `36_fig3A_chord_LR.R` — chord with individual CCL ligand–receptor pairs as sectors. **Confirmed NOT the published panel** (the manuscript uses cell-type sectors, not per-L-R sectors).
- `40_fig3A_pathway_chord.R` — pathway-level chord (precursor to 44).
- `44b_fig3A_chord_cellchat_style.R` — CellChat/MDPI-style chord with cell-type arcs subdivided by pathway (alternate styling of 44).

## Figure 4 — CCL signaling dot plots and spatial validation

| Panel | Description | Script(s) |
|-------|-------------|-----------|
| 4A | Top-10 CCL signaling dot plot, eos → mac | `34_fig4AB_ccl_dotplot.R` — the "Eosinophil → Macrophage (CCL)" plot |
| 4B | Top-10 CCL signaling dot plot, mac → eos | `34_fig4AB_ccl_dotplot.R` — the "Macrophage → Eosinophil (CCL)" plot |
| 4C | RNAscope confocal | Experimental / not computational |
| 4D | Imaris 3D reconstruction | Experimental / not computational |
| 4E–F | CCR3/CCL6 positive-cell bar charts | Experimental imaging quantification / not computational |
| 4G–I | G-cross spatial analysis violin plots | Experimental imaging quantification / not computational |

> **Orientation (resolved).** Both panels come from `34_fig4AB_ccl_dotplot.R`. The script and the result-PDF filenames have A/B **swapped** relative to the final manuscript: the file named `fig4A_mac_to_eos_dotplot.pdf` renders "Macrophage → Eosinophil (CCL)", and `fig4B_eos_to_mac_dotplot.pdf` renders "Eosinophil → Macrophage (CCL)". The manuscript legend orders them **4A = eos→mac, 4B = mac→eos**. Map panels by plot title (content), not by filename.

## Figure 5 — Gene regulatory networks (SCENIC)

| Panel | Description | Script(s) |
|-------|-------------|-----------|
| 5A | Cluster-enriched regulon activity heatmap (z-scored AUCell), eos C0–C5 + mac clusters | `16b_scenic_eos_mac_heatmap.py` (canonical), `16c_per_cluster_heatmap.py` (alternate) |
| 5B | Eos GRN map (nodes = TF regulons) | `16d_grn_network.py` |
| 5C | Sex-differential regulon activity heatmap (5 TFs, WT/KO × M/F, infected excluded) | `16e_sex_diff_regulons.py` |
| 5D | Regulon-based cross-dataset similarity dendrogram vs Hu | `41b_hu_regulon_boxplot.py` (boxplots + dendrogram); `41_hu_regulon_comparison.py` (AUCell scoring of Hu) |
| 5E | Regulon activity dot plot across eos clusters + Hu depots | `42h_regulon_merged_dotplot.py` (canonical); `42g_regulon_dotplots.py`, `42f_regulon_venn_and_comparison.py` (related) |

**SCENIC upstream (feeds all Fig 5 panels).** Regulon inference and AUCell scoring are produced by, in order: `14_scenic_export.R`, `14b_export_mac_h5ad.R`, `14c_export_per_cluster_loom.R`, `14d_create_looms.py` (export) → `15_scenic_run.sh`, `15c_scenic_per_cluster.sh` (pySCENIC GRNBoost2/cisTarget/AUCell) → `16_scenic_visualize.py` (QC). Hu cross-dataset regulons use `42_scenic_hu.sh`, `42c_scenic_hu_combined.sh`, `42d_scenic_hu_bat.sh`, then `42b_regulon_compare_aucell.py` / `42e_regulon_unique_analysis.py`. Condition-resolved regulon AUC violins (`22_scenic_by_condition.py`, `38_scenic_auc_violin.py`) support supplementary SCENIC analyses.
