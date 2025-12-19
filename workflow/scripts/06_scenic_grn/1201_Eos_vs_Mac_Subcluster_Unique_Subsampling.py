#!/usr/bin/env python3
"""
Find Eos subcluster-UNIQUE regulons compared to ALL Macrophages using Random Subsampling

Method based on Nature Immunology (https://www.nature.com/articles/s41590-020-0736-z#Sec19):
Randomly subsample cells from each condition and compare.

This script finds regulons that are:
1. ENRICHED in Eos compared to Mac (Eos subcluster AUC >> Mac AUC)
2. UNIQUE to specific Eos subclusters (high in ONE Eos subcluster, low in OTHER Eos subclusters)

Output:
- eos_enriched_heatmap: Regulons enriched in Eos vs Mac (organized by which Eos subcluster is max)
- eos_subcluster_unique_heatmap: Regulons that are BOTH enriched vs Mac AND unique to one Eos subcluster
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import logging
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Custom colormap
CUSTOM_CMAP = LinearSegmentedColormap.from_list('custom_blue_red', ['#67a9cf', '#ef8a62'])

# Subsampling parameters
CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42

# Cutoff thresholds
CUTOFF_LEVELS = {
    'relax': {
        'min_eos_auc': 0.15,           # Min AUC in max Eos subcluster
        'min_fold_change_vs_mac': 1.2,  # Fold change vs Mac mean
        'min_uniqueness_score': 0.03,   # Difference from other Eos subclusters
    },
    'moderate': {
        'min_eos_auc': 0.20,
        'min_fold_change_vs_mac': 1.5,
        'min_uniqueness_score': 0.05,
    },
    'stringent': {
        'min_eos_auc': 0.25,
        'min_fold_change_vs_mac': 2.0,
        'min_uniqueness_score': 0.08,
    }
}


def load_auc_matrices(scenic_dir):
    """Load AUC matrices for Eos subclusters and Mac subtypes"""

    eos_clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
    eos_data = {}

    for cluster in eos_clusters:
        auc_file = scenic_dir / "eos" / cluster / "step3_auc_matrix.csv"
        if auc_file.exists():
            eos_data[cluster] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {cluster}: {eos_data[cluster].shape[0]} cells, {eos_data[cluster].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    mac_types = ['M1', 'M2_Functional', 'M2_Surface']
    mac_data = {}

    for mac_type in mac_types:
        auc_file = scenic_dir / "macrophages" / mac_type / "step3_auc_matrix.csv"
        if auc_file.exists():
            mac_data[mac_type] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {mac_type}: {mac_data[mac_type].shape[0]} cells, {mac_data[mac_type].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    eos_total = sum(df.shape[0] for df in eos_data.values())
    mac_total = sum(df.shape[0] for df in mac_data.values())
    logging.info(f"Total Eos cells: {eos_total}, Total Mac cells: {mac_total}")

    if eos_total == 0 or mac_total == 0:
        raise FileNotFoundError("Could not load AUC matrices.")

    loaded_eos = [c for c in eos_clusters if c in eos_data]
    loaded_mac = [m for m in mac_types if m in mac_data]

    return eos_data, mac_data, loaded_eos, loaded_mac


def load_gene_counts(scenic_dir, eos_clusters, mac_types):
    """Load gene counts from step2_regulons.csv files"""
    gene_counts = {}

    for cluster in eos_clusters:
        regulon_path = scenic_dir / "eos" / cluster / "step2_regulons.csv"
        if regulon_path.exists():
            try:
                regulons_df = pd.read_csv(regulon_path)
                if 'regulon' in regulons_df.columns:
                    for regulon_name in regulons_df['regulon'].unique():
                        if '(' in regulon_name and ')' in regulon_name:
                            tf_name = regulon_name.split('(')[0]
                            count = int(regulon_name.split('(')[1].split(')')[0])
                            if tf_name not in gene_counts or count > gene_counts[tf_name]:
                                gene_counts[tf_name] = count
            except:
                pass

    for mac_type in mac_types:
        regulon_path = scenic_dir / "macrophages" / mac_type / "step2_regulons.csv"
        if regulon_path.exists():
            try:
                regulons_df = pd.read_csv(regulon_path)
                if 'regulon' in regulons_df.columns:
                    for regulon_name in regulons_df['regulon'].unique():
                        if '(' in regulon_name and ')' in regulon_name:
                            tf_name = regulon_name.split('(')[0]
                            count = int(regulon_name.split('(')[1].split(')')[0])
                            if tf_name not in gene_counts or count > gene_counts[tf_name]:
                                gene_counts[tf_name] = count
            except:
                pass

    logging.info(f"Loaded gene counts for {len(gene_counts)} TFs")
    return gene_counts


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def get_all_regulons(eos_data, mac_data):
    """Get all unique regulons across Eos and Mac"""
    all_regulons = set()
    for df in eos_data.values():
        all_regulons.update(df.columns)
    for df in mac_data.values():
        all_regulons.update(df.columns)
    logging.info(f"Found {len(all_regulons)} total unique regulons")
    return list(all_regulons)


def calculate_mean_auc_subsampling(eos_data, mac_data, eos_clusters, mac_types,
                                    cells_per_cluster=200):
    """
    Calculate mean AUC using simple random subsampling for both Eos and Mac.
    Subsample cells_per_cluster cells from each condition once.
    """
    logging.info(f"\n=== Random Subsampling ===")
    logging.info(f"Cells per cluster: {cells_per_cluster}")

    np.random.seed(RANDOM_SEED)

    # Build TF mapping for Eos
    eos_tf_to_regulon = {}
    all_tfs = set()
    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        eos_tf_to_regulon[cluster] = {}
        for col in eos_data[cluster].columns:
            tf_name = extract_tf_name(col)
            eos_tf_to_regulon[cluster][tf_name] = col
            all_tfs.add(tf_name)

    # Build TF mapping for Mac
    mac_tf_to_regulon = {}
    for mac_type in mac_types:
        if mac_type not in mac_data:
            continue
        mac_tf_to_regulon[mac_type] = {}
        for col in mac_data[mac_type].columns:
            tf_name = extract_tf_name(col)
            mac_tf_to_regulon[mac_type][tf_name] = col
            all_tfs.add(tf_name)

    all_tfs = sorted(list(all_tfs))
    logging.info(f"Total unique TFs: {len(all_tfs)}")

    # Report cell counts
    cell_counts = {}
    logging.info("\nEos subclusters:")
    for cluster in eos_clusters:
        if cluster in eos_data:
            n_cells = eos_data[cluster].shape[0]
            cell_counts[cluster] = n_cells
            status = "OK" if n_cells >= cells_per_cluster else "will sample with replacement"
            logging.info(f"  {cluster}: {n_cells} cells ({status})")

    logging.info("\nMac subtypes:")
    for mac_type in mac_types:
        if mac_type in mac_data:
            n_cells = mac_data[mac_type].shape[0]
            cell_counts[mac_type] = n_cells
            status = "OK" if n_cells >= cells_per_cluster else "will sample with replacement"
            logging.info(f"  {mac_type}: {n_cells} cells ({status})")

    # Subsample each Eos cluster once
    eos_auc = {cluster: {} for cluster in eos_clusters}
    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        cluster_df = eos_data[cluster]
        n_cells = cluster_df.shape[0]

        if n_cells >= cells_per_cluster:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=False)
        else:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=True)

        sampled_df = cluster_df.iloc[sampled_idx]
        logging.info(f"  Subsampled {cluster}: {len(sampled_idx)} cells")

        for tf in all_tfs:
            if tf in eos_tf_to_regulon.get(cluster, {}):
                regulon_col = eos_tf_to_regulon[cluster][tf]
                eos_auc[cluster][tf] = sampled_df[regulon_col].mean()
            else:
                eos_auc[cluster][tf] = 0.0

    # Subsample each Mac subtype once
    mac_auc = {mac_type: {} for mac_type in mac_types}
    for mac_type in mac_types:
        if mac_type not in mac_data:
            continue
        mac_df = mac_data[mac_type]
        n_cells = mac_df.shape[0]

        if n_cells >= cells_per_cluster:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=False)
        else:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=True)

        sampled_df = mac_df.iloc[sampled_idx]
        logging.info(f"  Subsampled {mac_type}: {len(sampled_idx)} cells")

        for tf in all_tfs:
            if tf in mac_tf_to_regulon.get(mac_type, {}):
                regulon_col = mac_tf_to_regulon[mac_type][tf]
                mac_auc[mac_type][tf] = sampled_df[regulon_col].mean()
            else:
                mac_auc[mac_type][tf] = 0.0

    eos_mean_auc = pd.DataFrame(eos_auc).reindex(all_tfs)
    mac_mean_auc = pd.DataFrame(mac_auc).reindex(all_tfs)

    logging.info(f"Subsampling complete")

    return eos_mean_auc, mac_mean_auc, cell_counts, all_tfs


def find_eos_enriched_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types,
                                       cutoff_level='relax'):
    """
    Find regulons that are:
    1. ENRICHED in Eos vs Mac (max Eos subcluster >> Mac mean)
    2. UNIQUE to specific Eos subclusters (high in ONE subcluster, low in others)
    """

    criteria = CUTOFF_LEVELS[cutoff_level]
    min_eos_auc = criteria['min_eos_auc']
    min_fc_vs_mac = criteria['min_fold_change_vs_mac']
    min_uniqueness = criteria['min_uniqueness_score']

    all_tfs = list(eos_mean_auc.index)
    results = []

    for tf in all_tfs:
        # Get Eos subcluster AUCs
        eos_aucs = {cluster: eos_mean_auc.loc[tf, cluster] for cluster in eos_clusters
                    if cluster in eos_mean_auc.columns}

        # Get Mac AUCs
        mac_aucs = {mac_type: mac_mean_auc.loc[tf, mac_type] for mac_type in mac_types
                    if mac_type in mac_mean_auc.columns}

        if not eos_aucs or not mac_aucs:
            continue

        # Find max Eos subcluster
        max_eos_cluster = max(eos_aucs, key=eos_aucs.get)
        max_eos_auc = eos_aucs[max_eos_cluster]

        # Calculate Eos mean (across subclusters)
        eos_mean = np.mean(list(eos_aucs.values()))

        # Calculate Mac mean
        mac_mean = np.mean(list(mac_aucs.values()))

        # Calculate fold change vs Mac
        fold_change_vs_mac = max_eos_auc / mac_mean if mac_mean > 0 else 0

        # Calculate uniqueness score (max Eos - second max Eos)
        sorted_eos_aucs = sorted(eos_aucs.values(), reverse=True)
        if len(sorted_eos_aucs) >= 2:
            uniqueness_score = sorted_eos_aucs[0] - sorted_eos_aucs[1]
        else:
            uniqueness_score = sorted_eos_aucs[0]

        # Other Eos subclusters mean
        other_eos_aucs = [v for k, v in eos_aucs.items() if k != max_eos_cluster]
        max_other_eos = max(other_eos_aucs) if other_eos_aucs else 0
        mean_other_eos = np.mean(other_eos_aucs) if other_eos_aucs else 0

        results.append({
            'tf': tf,
            'max_eos_cluster': max_eos_cluster,
            'max_eos_auc': max_eos_auc,
            'eos_mean': eos_mean,
            'mac_mean': mac_mean,
            'fold_change_vs_mac': fold_change_vs_mac,
            'uniqueness_score': uniqueness_score,
            'max_other_eos_auc': max_other_eos,
            'mean_other_eos_auc': mean_other_eos,
            **{f'eos_{c}': eos_aucs.get(c, 0) for c in eos_clusters},
            **{f'mac_{m}': mac_aucs.get(m, 0) for m in mac_types}
        })

    results_df = pd.DataFrame(results)

    # Filter for Eos-enriched regulons
    eos_enriched = results_df[
        (results_df['max_eos_auc'] >= min_eos_auc) &
        (results_df['fold_change_vs_mac'] >= min_fc_vs_mac)
    ].copy()

    # Filter for Eos subcluster-UNIQUE regulons (enriched vs Mac AND unique within Eos)
    eos_unique = results_df[
        (results_df['max_eos_auc'] >= min_eos_auc) &
        (results_df['fold_change_vs_mac'] >= min_fc_vs_mac) &
        (results_df['uniqueness_score'] >= min_uniqueness)
    ].copy()

    # Sort by uniqueness score
    eos_enriched = eos_enriched.sort_values('fold_change_vs_mac', ascending=False)
    eos_unique = eos_unique.sort_values('uniqueness_score', ascending=False)

    logging.info(f"\n{cutoff_level.upper()} cutoff:")
    logging.info(f"  Eos-enriched vs Mac: {len(eos_enriched)} regulons")
    logging.info(f"  Eos subcluster-UNIQUE: {len(eos_unique)} regulons")

    # Summary by cluster
    if len(eos_unique) > 0:
        cluster_counts = eos_unique['max_eos_cluster'].value_counts()
        for cluster in eos_clusters:
            count = cluster_counts.get(cluster, 0)
            logging.info(f"    {cluster}: {count} unique regulons")

    return results_df, eos_enriched, eos_unique


def create_heatmap(regulon_df, eos_clusters, mac_types, output_path, gene_counts,
                   title="Eos-Enriched Regulons", group_by_cluster=True):
    """Create heatmap organized by Eos subcluster"""

    if len(regulon_df) == 0:
        logging.warning(f"No regulons to plot for {title}")
        return

    # Order by max cluster
    if group_by_cluster:
        ordered_dfs = []
        for cluster in eos_clusters:
            cluster_df = regulon_df[regulon_df['max_eos_cluster'] == cluster].copy()
            if len(cluster_df) > 0:
                cluster_df = cluster_df.sort_values('uniqueness_score', ascending=False)
                ordered_dfs.append(cluster_df)
        if ordered_dfs:
            regulon_df = pd.concat(ordered_dfs, ignore_index=True)

    # Build plot data
    all_cols = eos_clusters + mac_types
    data_cols = [f'eos_{c}' if c in eos_clusters else f'mac_{c}' for c in all_cols]

    plot_data = regulon_df[data_cols].values

    # Y-axis labels
    ylabels = []
    for _, row in regulon_df.iterrows():
        tf = row['tf']
        if tf in gene_counts:
            ylabels.append(f"{tf}({gene_counts[tf]})")
        else:
            ylabels.append(tf)

    # Create figure
    fig_height = max(8, len(regulon_df) * 0.35)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    im = ax.imshow(plot_data, cmap=CUSTOM_CMAP, aspect='auto', vmin=0, vmax=1.0)

    cbar = plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label('Mean AUC (subsampled)', fontsize=12, fontweight='bold')

    ax.set_xticks(np.arange(len(all_cols)))
    ax.set_xticklabels(all_cols, rotation=45, ha='right', fontsize=10)

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=8)

    # Add vertical line to separate Eos and Mac
    ax.axvline(x=len(eos_clusters) - 0.5, color='white', linewidth=2)

    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_ylabel('Regulon (gene count)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_summary_barplot(eos_unique_df, eos_clusters, output_path, cutoff_level):
    """Create bar plot showing unique regulons per Eos subcluster"""

    cluster_counts = eos_unique_df['max_eos_cluster'].value_counts()

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.Set2(np.linspace(0, 1, len(eos_clusters)))

    counts = [cluster_counts.get(c, 0) for c in eos_clusters]
    bars = ax.bar(eos_clusters, counts, color=colors, edgecolor='black')

    ax.set_xlabel('Eos Subcluster', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Unique Regulons', fontsize=12, fontweight='bold')
    ax.set_title(f'Eos Subcluster-Unique Regulons (vs Mac)\n({cutoff_level})',
                 fontsize=14, fontweight='bold')

    for bar, val in zip(bars, counts):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    str(val), ha='center', fontsize=11, fontweight='bold')

    max_val = max(counts) if counts else 1
    ax.set_ylim(0, max(max_val * 1.15, 5))
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_cell_counts_plot(cell_counts, cells_per_cluster, output_path):
    """Create cell counts summary plot"""

    fig, ax = plt.subplots(figsize=(12, 6))

    clusters = list(cell_counts.keys())
    counts = list(cell_counts.values())

    colors = ['#2ecc71' if c >= cells_per_cluster else '#e74c3c' for c in counts]

    bars = ax.bar(clusters, counts, color=colors, edgecolor='black')
    ax.axhline(y=cells_per_cluster, color='blue', linestyle='--', linewidth=2,
               label=f'Sampling size ({cells_per_cluster})')

    ax.set_xlabel('Cluster/Subtype', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Cell Counts vs Sampling Size', fontsize=14, fontweight='bold')

    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                str(val), ha='center', fontsize=9, fontweight='bold')

    ax.legend()
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Find Eos subcluster-unique regulons vs Mac using subsampling'
    )
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--output-dir', type=str, default='results_eos_vs_mac_unique',
                        help='Path to output directory')
    parser.add_argument('--cells-per-cluster', type=int, default=200,
                        help='Number of cells to subsample per cluster')
    parser.add_argument('--cutoff-levels', type=str, nargs='+',
                        default=['relax', 'moderate', 'stringent'],
                        help='Cutoff levels to use')
    args = parser.parse_args()

    global CELLS_PER_CLUSTER
    CELLS_PER_CLUSTER = args.cells_per_cluster

    scenic_dir = Path(args.scenic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("EOS SUBCLUSTER-UNIQUE REGULONS vs MAC")
    logging.info("Method: Random Subsampling (Nature Immunology)")
    logging.info("="*80)
    logging.info(f"Cells per cluster: {CELLS_PER_CLUSTER}")

    # Load data
    logging.info("\n=== Loading AUC Matrices ===")
    eos_data, mac_data, eos_clusters, mac_types = load_auc_matrices(scenic_dir)

    logging.info("\n=== Loading Gene Counts ===")
    gene_counts = load_gene_counts(scenic_dir, eos_clusters, mac_types)

    # Calculate mean AUC using subsampling
    logging.info("\n=== Calculating Mean AUC (Subsampling) ===")
    eos_mean_auc, mac_mean_auc, cell_counts, all_tfs = calculate_mean_auc_subsampling(
        eos_data, mac_data, eos_clusters, mac_types,
        cells_per_cluster=CELLS_PER_CLUSTER
    )

    # Save mean AUC matrices
    eos_mean_auc.to_csv(output_dir / 'eos_mean_auc_subsampling.csv')
    mac_mean_auc.to_csv(output_dir / 'mac_mean_auc_subsampling.csv')

    # Create cell counts plot
    create_cell_counts_plot(cell_counts, CELLS_PER_CLUSTER, output_dir / 'cell_counts_summary.pdf')

    # Process each cutoff level
    logging.info("\n=== Finding Eos-Enriched and Unique Regulons ===")

    for cutoff_level in args.cutoff_levels:
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {cutoff_level.upper()}")
        logging.info(f"{'='*60}")

        results_df, eos_enriched, eos_unique = find_eos_enriched_unique_regulons(
            eos_mean_auc, mac_mean_auc, eos_clusters, mac_types,
            cutoff_level=cutoff_level
        )

        # Save full results
        results_df.to_csv(output_dir / f'all_regulons_metrics_{cutoff_level}.csv', index=False)

        # Save Eos-enriched
        if len(eos_enriched) > 0:
            eos_enriched.to_csv(output_dir / f'eos_enriched_{cutoff_level}.csv', index=False)
            create_heatmap(
                eos_enriched, eos_clusters, mac_types,
                output_dir / f'eos_enriched_heatmap_{cutoff_level}.pdf',
                gene_counts,
                title=f'Eos-Enriched Regulons vs Mac (n={len(eos_enriched)}, {cutoff_level})'
            )

        # Save Eos subcluster-UNIQUE
        if len(eos_unique) > 0:
            eos_unique.to_csv(output_dir / f'eos_subcluster_unique_{cutoff_level}.csv', index=False)
            create_heatmap(
                eos_unique, eos_clusters, mac_types,
                output_dir / f'eos_subcluster_unique_heatmap_{cutoff_level}.pdf',
                gene_counts,
                title=f'Eos Subcluster-UNIQUE Regulons (n={len(eos_unique)}, {cutoff_level})'
            )
            create_summary_barplot(
                eos_unique, eos_clusters,
                output_dir / f'eos_subcluster_unique_barplot_{cutoff_level}.pdf',
                cutoff_level
            )

        # Print summary
        print(f"\n{'='*60}")
        print(f"SUMMARY ({cutoff_level.upper()})")
        print("="*60)
        print(f"Eos-enriched vs Mac: {len(eos_enriched)} regulons")
        print(f"Eos subcluster-UNIQUE: {len(eos_unique)} regulons")
        if len(eos_unique) > 0:
            print("\nUnique regulons per cluster:")
            cluster_counts = eos_unique['max_eos_cluster'].value_counts()
            for cluster in eos_clusters:
                print(f"  {cluster}: {cluster_counts.get(cluster, 0)}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)
    logging.info(f"Output directory: {output_dir}")

    # Print SCP command
    print("\n" + "="*80)
    print("SCP COMMAND TO DOWNLOAD RESULTS:")
    print("="*80)
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
