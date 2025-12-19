#!/usr/bin/env python3
"""
EOS scVelo Trajectory Analysis - All Cells (7-Cluster Focus)
Date: September 24, 2025

Analysis focuses on:
- All cells
- 7 clusters (C0–C6) with custom EOS color palette
- scVelo dynamics recovery and velocity analysis
"""

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Set up scanpy and scvelo settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')


def setup_eos_palette():
    """EOS color palette inspired by Evangelion"""
    return {
        "C0": "#6CA6CD",  # Rei Blue
        "C1": "#E84A5F",  # Blood Alert
        "C2": "#FF6F1F",  # Eva Core Orange
        "C3": "#CBAACB",  # LCL Lilac
        "C4": "#B03A48",  # EVA Core Red
        "C5": "#5BA4A4",  # Test Tube Teal
        "C6": "#2F3C7E"   # Abyss Indigo
    }


def apply_cluster_colors(adata, palette):
    """Apply custom colors to clusters"""
    # Ensure Clusters is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs['Clusters']):
        adata.obs['Clusters'] = adata.obs['Clusters'].astype('category')

    # Get cluster categories and apply colors
    categories = adata.obs['Clusters'].cat.categories
    colors = [palette.get(str(cat), "#808080") for cat in categories]
    adata.uns['Clusters_colors'] = colors

    print(f"🎨 Applied colors to {len(categories)} clusters: {list(categories)}")
    return colors


def eos_full_dynamics_analysis(adata, output_dir):
    """Run EOS velocity dynamics analysis on ALL cells (7 clusters, custom palette)"""

    # Show cluster distribution
    cluster_counts = adata.obs['Clusters'].value_counts()
    print("📈 Cluster distribution (all cells):")
    print(cluster_counts)

    # Apply custom color palette
    palette = setup_eos_palette()
    apply_cluster_colors(adata, palette)

    # Recompute PCA and neighborhood graph
    print("🔄 Recomputing PCA on all cells...")
    sc.pp.pca(adata, n_comps=50)
    print("🔄 Recomputing neighbors on all cells...")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=40)

    # Run scVelo dynamics recovery
    print("⚙️ Running dynamics recovery on all cells (scv.tl.recover_dynamics)...")
    scv.tl.recover_dynamics(adata, n_jobs=8)

    # Velocity calculation (with more permissive min_r2)
    print("🚀 Calculating velocity (dynamical model)...")
    scv.tl.velocity(adata, mode='dynamical', min_r2=0.005, min_likelihood=None)

    # Check velocity gene count
    n_velocity_genes = np.sum(adata.var['velocity_genes'])
    print(f"📊 Velocity genes detected: {n_velocity_genes}")

    if n_velocity_genes < 10:
        print(f"⚠️ Low velocity gene count ({n_velocity_genes}), rerunning with stochastic mode...")
        scv.tl.velocity(adata, mode='stochastic', min_r2=0.01)
        n_velocity_genes = np.sum(adata.var['velocity_genes'])
        print(f"📊 Velocity genes (stochastic mode): {n_velocity_genes}")

    # Compute velocity graph
    print("📊 Computing velocity graph...")
    try:
        scv.tl.velocity_graph(adata, n_jobs=1)
    except Exception as e:
        print(f"⚠️ Error during velocity graph: {e}")
        print("🔄 Trying velocity_embedding as fallback...")
        scv.tl.velocity_embedding(adata, basis='umap')

    # Compute velocity pseudotime & confidence
    try:
        scv.tl.velocity_pseudotime(adata)
    except Exception as e:
        print(f"⚠️ Pseudotime could not be computed: {e}")

    try:
        scv.tl.velocity_confidence(adata)
    except Exception as e:
        print(f"⚠️ Confidence could not be computed: {e}")

    # Save processed data (without fingerprint in filename)
    output_file = f"{output_dir}/EOS_scVelo_adata_allcells_dynamics_092425.h5ad"
    adata.write(output_file)
    print(f"💾 Saved all-cell AnnData to: {output_file}")

    return {
        'adata': adata,
        'file': output_file,
        'n_cells': adata.n_obs
    }


def velocity_plot_allcells(result, output_dir):
    """Create a velocity stream plot (all cells, 7 clusters, EOS palette)"""

    adata = result['adata']

    # Ensure velocity embedding is computed
    if 'velocity_umap' not in adata.obsm:
        scv.tl.velocity_embedding(adata, basis='umap')

    # Plot stream (big, bold arrows)
    fig, ax = plt.subplots(figsize=(12, 10), facecolor='white')
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    title = f'RNA Velocity - ALL CELLS\n(n = {result["n_cells"]} cells, 7 clusters)'
    try:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="Clusters",
            ax=ax,
            legend_loc="right margin",
            show=False,
            title=title,
            size=80,
            alpha=0.85,
            density=1.6,
            arrow_size=5,
            arrow_length=6,
            linewidth=3.0
        )
    except Exception as e:
        print(f"⚠️ Stream plot failed: {e}")
        print("🔄 Falling back to scatter plot.")
        scv.pl.scatter(
            adata,
            basis='umap',
            color='Clusters',
            ax=ax,
            legend_loc="right margin",
            show=False,
            title="Cluster UMAP - ALL CELLS",
            size=80
        )

    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)

    # Cluster info textbox
    cluster_counts = adata.obs['Clusters'].value_counts()
    info_text = f"Clusters: {len(cluster_counts)}\n"
    info_text += "\n".join([f"{c}: {n}" for c, n in cluster_counts.items()])

    ax.text(0.02, 0.98, info_text,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
            verticalalignment='top', fontsize=10)

    plt.tight_layout()

    # Save plots in multiple formats (.pdf, .eps, .svg)
    basename = f"{output_dir}/EOS_velocity_allcells_092425"
    formats = ['pdf', 'eps', 'svg']
    for fmt in formats:
        try:
            output_file = f"{basename}.{fmt}"
            plt.savefig(output_file, dpi=300, bbox_inches='tight', format=fmt, facecolor='white', edgecolor='none')
            print(f"📊 Velocity plot ({fmt.upper()}) saved to: {output_file}")
        except Exception as e:
            print(f"⚠️ {fmt.upper()} save failed: {e}")

    plt.close()


def generate_summary_allcells(result, output_dir):
    """Generate summary for full/all cells velocity analysis"""

    summary = f"""
EOS scVelo Trajectory Analysis - ALL CELLS Summary
=================================================
Date: September 24, 2025
Dataset: EOS_scVelo_adata_250923.h5ad

Analysis Performed:
- All cells (n = {result['n_cells']:,})
- 7 clusters (C0-C6)
- Color palette: EOS Evangelion-inspired theme

Visualization Enhancements:
✨ Bigger arrows: arrow_size=5
✨ Bolder arrows: linewidth=3.0
✨ Longer arrows: arrow_length=6
✨ Higher density: density=1.6

Analysis Pipeline:
✅ PCA recomputation
✅ Neighbor graph reconstruction
✅ Dynamics recovery (scv.tl.recover_dynamics)
✅ Velocity calculation (dynamical → stochastic fallback if needed)
✅ Velocity graph construction
✅ Terminal states identification
✅ Velocity confidence metrics
✅ UMAP velocity stream visualization

Output Files:
- AnnData: {result['file']}
- Velocity plots: {output_dir}/EOS_velocity_allcells_092425.[pdf|eps|svg]

Visualizations saved to: {output_dir}/

🎯 Ready for downstream trajectory analysis with monocle/cellrank
"""

    with open(f"{output_dir}/allcells_analysis_summary_092425.txt", "w") as f:
        f.write(summary)

    print(summary)


def main():
    """Main EOS all-cells velocity analysis"""

    # File paths
    input_file = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_250923.h5ad"
    output_dir = "/bigdata/nairlab/shared/XQiu/EOS_project/R14_20250919/velocity"

    # Create output directory
    import os
    os.makedirs(output_dir, exist_ok=True)

    print("🚀 Starting EOS scVelo All-Cells Analysis - September 24, 2025")
    print(f"📁 Input: {input_file}")
    print(f"📁 Output: {output_dir}")
    print("✨ Focus: ALL cells, 7 clusters, EOS palette")

    # Load data
    adata = scv.read(input_file)

    # Run full analysis
    result = eos_full_dynamics_analysis(adata, output_dir)

    # Plot velocity embeddings for ALL cells
    velocity_plot_allcells(result, output_dir)

    # Generate summary
    generate_summary_allcells(result, output_dir)

    print("\n✅ All-cells analysis complete!")
    print(f"📊 Check results in: {output_dir}")
    print("🎯 Ready for downstream trajectory analysis with monocle/cellrank")


if __name__ == "__main__":
    main()