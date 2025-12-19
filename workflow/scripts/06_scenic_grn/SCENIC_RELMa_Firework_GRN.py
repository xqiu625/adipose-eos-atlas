#!/usr/bin/env python3
"""
Firework-style GRN Network for RELMα Sex-Differential Regulons
===============================================================
Creates radial "firework" networks where:
- Center: TF/Regulon name
- Outer nodes: Target genes regulated by that TF
- Color: Expression change (Female vs Male) - Blue-White-Red
- WT: More spread/color difference (larger sex differences)
- KO: Less spread/color difference (diminished sex differences)

Style: Clean publication figure matching Frontiers Fig 3B
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
from pathlib import Path
import ast
import re
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Z-score colormap: Blue (Male-enriched) - White - Red (Female-enriched)
ZSCORE_CMAP = LinearSegmentedColormap.from_list('zscore', ['#4575b4', '#f7f7f7', '#d73027'])

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Paths
SCENIC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned")
CONDITION_AUC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/condition_analysis")
OUTPUT_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/relma_sex_effect")
H5AD_DIR = Path("/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/scenic_project/inputs")

# Key sex-differential TFs from previous analysis
KEY_TFS = ['Fli1', 'Tef', 'Rara', 'Rel', 'Irf4']


def load_regulon_targets(scenic_dir, cell_type='eos'):
    """
    Load target genes for each regulon from SCENIC regulons files.
    Returns dict: {tf_name: [list of target genes]}

    File format: tf,target_genes,regulon
    Example: Fli1,Gene1;Gene2;Gene3,Fli1(3)
    """
    logging.info("Loading regulon target genes...")

    tf_targets = {}

    clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']

    for cluster in clusters:
        # Try step2_regulons.csv first (the actual filename on cluster)
        regulon_file = scenic_dir / cell_type / cluster / 'step2_regulons.csv'
        if not regulon_file.exists():
            regulon_file = scenic_dir / cell_type / cluster / 'regulons.csv'
        if not regulon_file.exists():
            regulon_file = scenic_dir / cell_type / cluster / f'{cluster}_regulons.csv'
        if not regulon_file.exists():
            logging.warning(f"No regulons file found for {cluster}")
            continue

        try:
            df = pd.read_csv(regulon_file)
            logging.info(f"Loading {cluster}: {len(df)} regulons")

            for _, row in df.iterrows():
                tf = row['tf'] if 'tf' in row else row.iloc[0]
                if pd.isna(tf):
                    continue

                # Parse target genes (semicolon-separated)
                targets_str = row['target_genes'] if 'target_genes' in row else row.iloc[1]

                if pd.notna(targets_str) and str(targets_str) != 'nan':
                    # Split by semicolon
                    genes = [g.strip() for g in str(targets_str).split(';') if g.strip()]
                    # Remove any scientific notation artifacts like 'e-06'
                    genes = [g for g in genes if not g.startswith('e-') and not g.startswith('e+')]

                    if genes:
                        if tf not in tf_targets:
                            tf_targets[tf] = set()
                        tf_targets[tf].update(genes)

        except Exception as e:
            logging.warning(f"Error loading {regulon_file}: {e}")
            continue

    # Convert sets to lists
    tf_targets = {k: list(v) for k, v in tf_targets.items()}
    logging.info(f"Loaded targets for {len(tf_targets)} TFs")

    return tf_targets


def load_gene_expression(h5ad_dir, cell_type='eos'):
    """
    Load gene expression data from h5ad file.
    Returns AnnData object with cells x genes.
    """
    import scanpy as sc

    if cell_type == 'eos':
        h5ad_file = h5ad_dir / 'eos_all.h5ad'
    else:
        h5ad_file = h5ad_dir / 'mac_all.h5ad'

    logging.info(f"Loading expression from: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    logging.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def get_condition_from_obs(adata):
    """
    Create condition labels from adata.obs metadata.
    Format: {genotype}_{Sex}_HFD
    """
    conditions = []
    for idx in range(adata.n_obs):
        genotype = adata.obs['genotype'].iloc[idx] if 'genotype' in adata.obs else 'Unknown'
        sex = adata.obs['Sex'].iloc[idx] if 'Sex' in adata.obs else 'Unknown'
        conditions.append(f"{genotype}_{sex}_HFD")
    return conditions


def calculate_gene_sex_difference_adata(adata, genes, female_cond, male_cond):
    """
    Calculate sex difference (Female - Male) for each gene from AnnData.
    Returns dict: {gene: difference}

    female_cond/male_cond: dict with keys like 'genotype', 'Sex'
    Example: {'genotype': 'WT', 'Sex': 'Female'}
    """
    import scipy.sparse as sp

    # Build masks
    f_mask = np.ones(adata.n_obs, dtype=bool)
    for key, val in female_cond.items():
        if key in adata.obs.columns:
            f_mask &= (adata.obs[key] == val).values

    m_mask = np.ones(adata.n_obs, dtype=bool)
    for key, val in male_cond.items():
        if key in adata.obs.columns:
            m_mask &= (adata.obs[key] == val).values

    logging.info(f"Female cells: {f_mask.sum()}, Male cells: {m_mask.sum()}")

    gene_diffs = {}
    gene_names = list(adata.var_names)

    for gene in genes:
        if gene in gene_names:
            gene_idx = gene_names.index(gene)

            # Get expression values
            if sp.issparse(adata.X):
                f_expr = adata.X[f_mask, gene_idx].toarray().flatten()
                m_expr = adata.X[m_mask, gene_idx].toarray().flatten()
            else:
                f_expr = adata.X[f_mask, gene_idx]
                m_expr = adata.X[m_mask, gene_idx]

            f_mean = np.mean(f_expr) if len(f_expr) > 0 else 0
            m_mean = np.mean(m_expr) if len(m_expr) > 0 else 0
            gene_diffs[gene] = f_mean - m_mean

    return gene_diffs


def create_firework_network(tf_name, target_genes, gene_diffs, ax, title,
                            max_targets=50, vmin=-0.1, vmax=0.1):
    """
    Create a single firework-style network for one TF.

    - Center: TF name (large node)
    - Radiating: Target genes (smaller nodes)
    - Color: Sex difference (Blue=Male, Red=Female)
    """
    G = nx.Graph()

    # Add TF as center node
    G.add_node(tf_name, node_type='tf', diff=0)

    # Filter and sort target genes by absolute difference
    valid_targets = [(g, gene_diffs.get(g, 0)) for g in target_genes if g in gene_diffs]
    valid_targets.sort(key=lambda x: abs(x[1]), reverse=True)

    # Take top targets
    top_targets = valid_targets[:max_targets]

    # Add target genes as nodes
    for gene, diff in top_targets:
        G.add_node(gene, node_type='target', diff=diff)
        G.add_edge(tf_name, gene)

    if G.number_of_nodes() <= 1:
        ax.text(0.5, 0.5, f'No target data\nfor {tf_name}',
               ha='center', va='center', fontsize=8)
        ax.axis('off')
        return

    # Create radial layout - TF in center, targets around
    pos = {}
    pos[tf_name] = (0, 0)  # Center

    n_targets = len(top_targets)
    for i, (gene, _) in enumerate(top_targets):
        angle = 2 * np.pi * i / n_targets
        radius = 1.0
        pos[gene] = (radius * np.cos(angle), radius * np.sin(angle))

    # Normalize differences for coloring
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Node colors
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        diff = G.nodes[node]['diff']
        if G.nodes[node]['node_type'] == 'tf':
            node_colors.append('#333333')  # Dark gray for TF
            node_sizes.append(400)
        else:
            node_colors.append(ZSCORE_CMAP(norm(diff)))
            node_sizes.append(60)

    # Draw edges (thin gray lines radiating from center)
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.3, alpha=0.4,
                           edge_color='#cccccc')

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.3)

    # Draw labels - TF in center larger, targets smaller
    tf_labels = {tf_name: tf_name}
    nx.draw_networkx_labels(G, pos, labels=tf_labels, ax=ax,
                            font_size=10, font_weight='bold')

    # Target gene labels (only show top ones to avoid clutter)
    if n_targets <= 20:
        target_labels = {g: g for g, _ in top_targets}
        # Offset labels slightly outward
        label_pos = {}
        for gene, _ in top_targets:
            x, y = pos[gene]
            # Push label outward
            dist = np.sqrt(x**2 + y**2)
            if dist > 0:
                label_pos[gene] = (x * 1.15, y * 1.15)
            else:
                label_pos[gene] = (x, y)
        nx.draw_networkx_labels(G, label_pos, labels=target_labels, ax=ax,
                                font_size=5, font_weight='normal')

    ax.set_title(title, fontsize=8, pad=5)
    ax.axis('off')
    ax.set_aspect('equal')


def create_wt_ko_comparison(tf_name, target_genes, wt_gene_diffs, ko_gene_diffs,
                            output_path, max_targets=40):
    """
    Create side-by-side WT vs KO firework networks for one TF.
    """
    # Figure size: 18cm x 8cm
    fig_width = 18 / 2.54
    fig_height = 9 / 2.54

    fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    # Calculate vmin/vmax from both conditions for consistent coloring
    all_diffs = list(wt_gene_diffs.values()) + list(ko_gene_diffs.values())
    if all_diffs:
        max_abs = max(abs(min(all_diffs)), abs(max(all_diffs)))
        vmin, vmax = -max_abs, max_abs
    else:
        vmin, vmax = -0.1, 0.1

    # WT panel
    wt_diffs_subset = {g: wt_gene_diffs.get(g, 0) for g in target_genes}
    create_firework_network(tf_name, target_genes, wt_diffs_subset, axes[0],
                            f'WT: {tf_name}', max_targets=max_targets,
                            vmin=vmin, vmax=vmax)

    # KO panel
    ko_diffs_subset = {g: ko_gene_diffs.get(g, 0) for g in target_genes}
    create_firework_network(tf_name, target_genes, ko_diffs_subset, axes[1],
                            f'KO: {tf_name}', max_targets=max_targets,
                            vmin=vmin, vmax=vmax)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=ZSCORE_CMAP, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Female - Male', fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    plt.suptitle(f'{tf_name} Target Genes: Sex Difference (Female - Male)',
                fontsize=9, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_multi_tf_figure(tf_list, tf_targets, wt_gene_diffs, ko_gene_diffs,
                           output_path, max_targets=30):
    """
    Create multi-panel figure with all key TFs.
    Layout: 2 rows (WT, KO) x N columns (TFs)
    """
    n_tfs = len(tf_list)

    # Figure size
    fig_width = min(4 * n_tfs, 20) / 2.54 * 2.54  # Scale with number of TFs
    fig_height = 16 / 2.54

    fig, axes = plt.subplots(2, n_tfs, figsize=(fig_width, fig_height))

    if n_tfs == 1:
        axes = axes.reshape(2, 1)

    # Get global vmin/vmax
    all_diffs = list(wt_gene_diffs.values()) + list(ko_gene_diffs.values())
    if all_diffs:
        max_abs = np.percentile([abs(d) for d in all_diffs], 95)
        vmin, vmax = -max_abs, max_abs
    else:
        vmin, vmax = -0.1, 0.1

    for i, tf in enumerate(tf_list):
        targets = tf_targets.get(tf, [])

        if not targets:
            axes[0, i].text(0.5, 0.5, f'No targets\nfor {tf}',
                           ha='center', va='center', fontsize=8)
            axes[0, i].axis('off')
            axes[1, i].text(0.5, 0.5, f'No targets\nfor {tf}',
                           ha='center', va='center', fontsize=8)
            axes[1, i].axis('off')
            continue

        # WT (top row)
        wt_diffs = {g: wt_gene_diffs.get(g, 0) for g in targets}
        create_firework_network(tf, targets, wt_diffs, axes[0, i],
                                f'WT: {tf}', max_targets=max_targets,
                                vmin=vmin, vmax=vmax)

        # KO (bottom row)
        ko_diffs = {g: ko_gene_diffs.get(g, 0) for g in targets}
        create_firework_network(tf, targets, ko_diffs, axes[1, i],
                                f'KO: {tf}', max_targets=max_targets,
                                vmin=vmin, vmax=vmax)

    # Row labels
    fig.text(0.02, 0.75, 'WT', fontsize=12, fontweight='bold', rotation=90, va='center')
    fig.text(0.02, 0.25, 'KO', fontsize=12, fontweight='bold', rotation=90, va='center')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=ZSCORE_CMAP, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.93, 0.3, 0.015, 0.4])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Female - Male\n(AUC difference)', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    plt.suptitle('Sex-Differential Regulon Target Networks\n'
                 'RELMα maintains sex-specific GRN programs',
                fontsize=10, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0.04, 0, 0.92, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_summary_figure(tf_list, tf_targets, wt_gene_diffs, ko_gene_diffs, output_path):
    """
    Create summary figure showing variance/spread reduction in KO.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12/2.54, 6/2.54))

    # Calculate variance of sex differences for each TF's targets
    wt_variances = []
    ko_variances = []
    tf_labels = []

    for tf in tf_list:
        targets = tf_targets.get(tf, [])
        if not targets:
            continue

        wt_diffs = [wt_gene_diffs.get(g, 0) for g in targets if g in wt_gene_diffs]
        ko_diffs = [ko_gene_diffs.get(g, 0) for g in targets if g in ko_gene_diffs]

        if wt_diffs and ko_diffs:
            wt_variances.append(np.var(wt_diffs))
            ko_variances.append(np.var(ko_diffs))
            tf_labels.append(tf)

    if not tf_labels:
        logging.warning("No data for summary figure")
        plt.close()
        return

    # Bar plot of variances
    x = np.arange(len(tf_labels))
    width = 0.35

    axes[0].bar(x - width/2, wt_variances, width, label='WT', color='#1a9850', alpha=0.8)
    axes[0].bar(x + width/2, ko_variances, width, label='KO', color='#d73027', alpha=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(tf_labels, fontsize=7)
    axes[0].set_ylabel('Variance of sex diff', fontsize=7)
    axes[0].legend(fontsize=6)
    axes[0].set_title('Target gene sex-difference variance', fontsize=8)
    axes[0].tick_params(labelsize=6)

    # Reduction percentage
    reductions = [(wt - ko) / wt * 100 if wt > 0 else 0
                  for wt, ko in zip(wt_variances, ko_variances)]

    colors = ['#d73027' if r > 0 else '#4575b4' for r in reductions]
    axes[1].bar(x, reductions, color=colors, alpha=0.8)
    axes[1].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(tf_labels, fontsize=7)
    axes[1].set_ylabel('Reduction (%)', fontsize=7)
    axes[1].set_title('Sex-difference reduction after KO', fontsize=8)
    axes[1].tick_params(labelsize=6)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved: {output_path}")
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Firework GRN for RELMα Sex-Differential Regulons')
    parser.add_argument('--scenic-dir', type=str, default=str(SCENIC_DIR))
    parser.add_argument('--h5ad-dir', type=str, default=str(H5AD_DIR))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('--cell-type', type=str, default='eos')
    parser.add_argument('--tfs', type=str, nargs='+', default=KEY_TFS,
                        help='List of TFs to visualize')
    parser.add_argument('--max-targets', type=int, default=40,
                        help='Maximum target genes to show per TF')
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    h5ad_dir = Path(args.h5ad_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cell_type = args.cell_type
    tf_list = args.tfs

    logging.info("="*80)
    logging.info("FIREWORK GRN VISUALIZATION")
    logging.info("="*80)

    # Load regulon targets
    tf_targets = load_regulon_targets(scenic_dir, cell_type)

    # Check which TFs have targets
    available_tfs = [tf for tf in tf_list if tf in tf_targets and len(tf_targets[tf]) > 0]
    logging.info(f"TFs with targets: {available_tfs}")

    for tf in tf_list:
        if tf in tf_targets:
            logging.info(f"  {tf}: {len(tf_targets[tf])} targets")
        else:
            logging.info(f"  {tf}: No targets found")

    # Load expression data from h5ad
    logging.info("\nLoading expression data from h5ad...")
    adata = load_gene_expression(h5ad_dir, cell_type)

    # Collect all target genes
    all_genes = set()
    for tf in available_tfs:
        all_genes.update(tf_targets[tf])
    logging.info(f"Total target genes to analyze: {len(all_genes)}")

    # Calculate sex differences for WT
    logging.info("\nCalculating WT sex differences...")
    wt_gene_diffs = calculate_gene_sex_difference_adata(
        adata, list(all_genes),
        female_cond={'genotype': 'WT', 'Sex': 'Female'},
        male_cond={'genotype': 'WT', 'Sex': 'Male'}
    )
    logging.info(f"WT: {len(wt_gene_diffs)} genes with data")

    # Calculate sex differences for KO
    logging.info("Calculating KO sex differences...")
    ko_gene_diffs = calculate_gene_sex_difference_adata(
        adata, list(all_genes),
        female_cond={'genotype': 'KO', 'Sex': 'Female'},
        male_cond={'genotype': 'KO', 'Sex': 'Male'}
    )
    logging.info(f"KO: {len(ko_gene_diffs)} genes with data")

    # Create visualizations
    logging.info("\n=== Creating Visualizations ===")

    # 1. Individual TF comparisons
    for tf in available_tfs:
        create_wt_ko_comparison(
            tf, tf_targets[tf], wt_gene_diffs, ko_gene_diffs,
            output_dir / f'{cell_type}_firework_{tf}.pdf',
            max_targets=args.max_targets
        )

    # 2. Multi-TF figure (all TFs in one figure)
    if len(available_tfs) > 0:
        create_multi_tf_figure(
            available_tfs, tf_targets, wt_gene_diffs, ko_gene_diffs,
            output_dir / f'{cell_type}_firework_all_tfs.pdf',
            max_targets=30
        )

    # 3. Summary figure
    create_summary_figure(
        available_tfs, tf_targets, wt_gene_diffs, ko_gene_diffs,
        output_dir / f'{cell_type}_firework_summary.pdf'
    )

    logging.info("\n" + "="*80)
    logging.info("COMPLETE")
    logging.info("="*80)
    print(f"\nOutput: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir}" .')


if __name__ == "__main__":
    main()
