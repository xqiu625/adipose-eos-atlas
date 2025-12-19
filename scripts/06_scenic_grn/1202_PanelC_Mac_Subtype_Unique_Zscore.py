#!/usr/bin/env python3
"""
Panel C: Mac subtype-unique regulons
Compare each Mac subtype with ALL OTHER Mac subtypes to get subtype-unique regulons.

This is DIFFERENT from 1202_Mac_vs_Eos:
- Previous: Mac subtypes vs Eos (cross cell-type comparison)
- This: Mac subtypes vs OTHER Mac subtypes (within cell-type comparison)

Output: Regulons uniquely active in specific Mac subtypes (not active in other Mac subtypes)
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

# Blue to Red colormap
ZSCORE_CMAP = LinearSegmentedColormap.from_list('blue_white_red', ['#2166ac', '#f7f7f7', '#b2182b'])

# Subsampling parameters
CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42

# Cutoff thresholds
CUTOFF_LEVELS = {
    'relax': {
        'min_auc': 0.15,
        'min_uniqueness_score': 0.03,
    },
    'moderate': {
        'min_auc': 0.20,
        'min_uniqueness_score': 0.05,
    },
    'stringent': {
        'min_auc': 0.25,
        'min_uniqueness_score': 0.08,
    }
}


def load_auc_matrices(scenic_dir):
    """Load AUC matrices for Mac subtypes and Eos subclusters"""

    mac_types = ['M1', 'M2_Functional', 'M2_Surface']
    mac_data = {}

    for mac_type in mac_types:
        auc_file = scenic_dir / "macrophages" / mac_type / "step3_auc_matrix.csv"
        if auc_file.exists():
            mac_data[mac_type] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {mac_type}: {mac_data[mac_type].shape[0]} cells, {mac_data[mac_type].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    mac_total = sum(df.shape[0] for df in mac_data.values())
    logging.info(f"Total Mac cells: {mac_total}")

    if mac_total == 0:
        raise FileNotFoundError("Could not load Mac AUC matrices.")

    loaded_mac = [m for m in mac_types if m in mac_data]

    # Also load Eos subclusters
    eos_clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
    eos_data = {}

    for cluster in eos_clusters:
        auc_file = scenic_dir / "eos" / cluster / "step3_auc_matrix.csv"
        if auc_file.exists():
            eos_data[cluster] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {cluster}: {eos_data[cluster].shape[0]} cells, {eos_data[cluster].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    eos_total = sum(df.shape[0] for df in eos_data.values())
    logging.info(f"Total Eos cells: {eos_total}")

    loaded_eos = [c for c in eos_clusters if c in eos_data]

    return mac_data, loaded_mac, eos_data, loaded_eos


def load_gene_counts(scenic_dir, mac_types):
    """Load gene counts from step2_regulons.csv files"""
    gene_counts = {}

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


def calculate_mean_auc_subsampling(mac_data, mac_types, eos_data, eos_clusters, cells_per_cluster=200):
    """
    Calculate mean AUC using simple random subsampling for both Mac and Eos.
    """
    logging.info(f"\n=== Random Subsampling ===")
    logging.info(f"Cells per cluster: {cells_per_cluster}")

    np.random.seed(RANDOM_SEED)

    # Build TF mapping for Mac
    mac_tf_to_regulon = {}
    all_tfs = set()

    for mac_type in mac_types:
        if mac_type not in mac_data:
            continue
        mac_tf_to_regulon[mac_type] = {}
        for col in mac_data[mac_type].columns:
            tf_name = extract_tf_name(col)
            mac_tf_to_regulon[mac_type][tf_name] = col
            all_tfs.add(tf_name)

    # Build TF mapping for Eos
    eos_tf_to_regulon = {}
    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        eos_tf_to_regulon[cluster] = {}
        for col in eos_data[cluster].columns:
            tf_name = extract_tf_name(col)
            eos_tf_to_regulon[cluster][tf_name] = col
            all_tfs.add(tf_name)

    all_tfs = sorted(list(all_tfs))
    logging.info(f"Total unique TFs: {len(all_tfs)}")

    # Report cell counts
    cell_counts = {}
    for mac_type in mac_types:
        if mac_type in mac_data:
            n_cells = mac_data[mac_type].shape[0]
            cell_counts[mac_type] = n_cells
            status = "OK" if n_cells >= cells_per_cluster else "will sample with replacement"
            logging.info(f"  {mac_type}: {n_cells} cells ({status})")

    for cluster in eos_clusters:
        if cluster in eos_data:
            n_cells = eos_data[cluster].shape[0]
            cell_counts[cluster] = n_cells
            status = "OK" if n_cells >= cells_per_cluster else "will sample with replacement"
            logging.info(f"  {cluster}: {n_cells} cells ({status})")

    # Subsample each Mac subtype
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

    # Subsample each Eos cluster
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

    mac_mean_auc = pd.DataFrame(mac_auc).reindex(all_tfs)
    eos_mean_auc = pd.DataFrame(eos_auc).reindex(all_tfs)
    logging.info(f"Subsampling complete")

    return mac_mean_auc, eos_mean_auc, cell_counts, all_tfs


def find_mac_subtype_unique_regulons(mac_mean_auc, eos_mean_auc, mac_types, eos_clusters, cutoff_level='relax'):
    """
    Find regulons unique to specific Mac subtypes.
    Compare each subtype to ALL OTHER subtypes.
    Also include Eos subcluster data for visualization.
    """

    criteria = CUTOFF_LEVELS[cutoff_level]
    min_auc = criteria['min_auc']
    min_uniqueness = criteria['min_uniqueness_score']

    all_tfs = list(mac_mean_auc.index)
    results = []

    for tf in all_tfs:
        # Get AUC for each subtype
        mac_aucs = {mac_type: mac_mean_auc.loc[tf, mac_type] for mac_type in mac_types
                    if mac_type in mac_mean_auc.columns}

        # Get AUC for each Eos subcluster
        eos_aucs = {}
        for cluster in eos_clusters:
            if cluster in eos_mean_auc.columns:
                eos_aucs[cluster] = eos_mean_auc.loc[tf, cluster]

        if not mac_aucs:
            continue

        # Find max subtype
        max_subtype = max(mac_aucs, key=mac_aucs.get)
        max_auc = mac_aucs[max_subtype]

        # Calculate mean across all subtypes
        mac_mean = np.mean(list(mac_aucs.values()))

        # Calculate uniqueness score (max - second max)
        sorted_aucs = sorted(mac_aucs.values(), reverse=True)
        if len(sorted_aucs) >= 2:
            uniqueness_score = sorted_aucs[0] - sorted_aucs[1]
            second_max_auc = sorted_aucs[1]
        else:
            uniqueness_score = sorted_aucs[0]
            second_max_auc = 0

        # Other subtypes stats
        other_aucs = [v for k, v in mac_aucs.items() if k != max_subtype]
        max_other = max(other_aucs) if other_aucs else 0
        mean_other = np.mean(other_aucs) if other_aucs else 0

        results.append({
            'tf': tf,
            'max_subtype': max_subtype,
            'max_auc': max_auc,
            'second_max_auc': second_max_auc,
            'mac_mean': mac_mean,
            'uniqueness_score': uniqueness_score,
            'max_other_auc': max_other,
            'mean_other_auc': mean_other,
            **{f'mac_{m}': mac_aucs.get(m, 0) for m in mac_types},
            **{f'eos_{c}': eos_aucs.get(c, 0) for c in eos_clusters}
        })

    results_df = pd.DataFrame(results)

    # Filter for subtype-unique regulons
    mac_unique = results_df[
        (results_df['max_auc'] >= min_auc) &
        (results_df['uniqueness_score'] >= min_uniqueness)
    ].copy().sort_values('uniqueness_score', ascending=False)

    logging.info(f"\n{cutoff_level.upper()} cutoff:")
    logging.info(f"  Mac subtype-UNIQUE regulons: {len(mac_unique)}")

    # Summary by subtype
    if len(mac_unique) > 0:
        subtype_counts = mac_unique['max_subtype'].value_counts()
        for subtype in mac_types:
            count = subtype_counts.get(subtype, 0)
            logging.info(f"    {subtype}: {count} unique regulons")

    return results_df, mac_unique


def create_zscore_heatmap(regulon_df, mac_types, eos_clusters, output_path, gene_counts,
                           title="Mac Subtype-Unique Regulons (Z-score)",
                           group_by_subtype=True, zscore_cap=2.5):
    """Create heatmap with Z-score scaled values including Mac and Eos"""

    if len(regulon_df) == 0:
        logging.warning(f"No regulons to plot for {title}")
        return

    # Order by max subtype
    if group_by_subtype:
        ordered_dfs = []
        for subtype in mac_types:
            subtype_df = regulon_df[regulon_df['max_subtype'] == subtype].copy()
            if len(subtype_df) > 0:
                subtype_df = subtype_df.sort_values('uniqueness_score', ascending=False)
                ordered_dfs.append(subtype_df)
        if ordered_dfs:
            regulon_df = pd.concat(ordered_dfs, ignore_index=True)

    # Build plot data - include both Mac and Eos
    all_cols = mac_types + eos_clusters
    data_cols = [f'mac_{m}' for m in mac_types] + [f'eos_{c}' for c in eos_clusters]

    raw_data = regulon_df[data_cols].values

    # Calculate row-wise z-score
    row_means = raw_data.mean(axis=1, keepdims=True)
    row_stds = raw_data.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    zscore_data = (raw_data - row_means) / row_stds

    # Filter out regulons where ANY Mac subtype z-score falls between 0.3 and 0.7
    # (i.e., keep only regulons where ALL Mac z-scores are outside 0.3-0.7 range)
    n_mac = len(mac_types)
    keep_mask = []
    for i in range(zscore_data.shape[0]):
        mac_zscores = zscore_data[i, :n_mac]  # Get Mac columns only
        # Check if any Mac z-score falls in 0.3-0.7 range
        has_middle_value = any((0.3 <= z <= 0.7) for z in mac_zscores)
        keep_mask.append(not has_middle_value)

    keep_mask = np.array(keep_mask)
    zscore_data = zscore_data[keep_mask]
    regulon_df = regulon_df.iloc[keep_mask].reset_index(drop=True)

    logging.info(f"After z-score filtering (removing Mac 0.3-0.7): {len(regulon_df)} regulons remain")

    if len(regulon_df) == 0:
        logging.warning(f"No regulons remain after z-score filtering for {title}")
        return

    zscore_data_capped = np.clip(zscore_data, -zscore_cap, zscore_cap)

    # Y-axis labels
    ylabels = []
    for _, row in regulon_df.iterrows():
        tf = row['tf']
        if tf in gene_counts:
            ylabels.append(f"{tf}({gene_counts[tf]})")
        else:
            ylabels.append(tf)

    # Save z-score matrix as CSV
    zscore_df = pd.DataFrame(zscore_data, index=regulon_df['tf'].values, columns=all_cols)
    csv_path = str(output_path).replace('.pdf', '_matrix.csv')
    zscore_df.to_csv(csv_path)
    logging.info(f"Saved z-score matrix: {csv_path}")

    # Create figure
    fig_height = max(8, len(regulon_df) * 0.35)
    fig_width = max(12, len(all_cols) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    im = ax.imshow(zscore_data_capped, cmap=ZSCORE_CMAP, aspect='auto',
                   vmin=-zscore_cap, vmax=zscore_cap)

    cbar = plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label('Z-score (row-scaled)', fontsize=12, fontweight='bold')

    ax.set_xticks(np.arange(len(all_cols)))
    ax.set_xticklabels(all_cols, rotation=45, ha='right', fontsize=10)

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=8)

    # Add vertical line to separate Mac and Eos
    ax.axvline(x=len(mac_types) - 0.5, color='black', linewidth=2)

    ax.set_xlabel('Cell Type / Subcluster', fontsize=12, fontweight='bold')
    ax.set_ylabel('Regulon (gene count)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_summary_barplot(mac_unique_df, mac_types, output_path, cutoff_level):
    """Create bar plot showing unique regulons per Mac subtype"""

    subtype_counts = mac_unique_df['max_subtype'].value_counts()

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, Blue, Green for M1, M2_Func, M2_Surf

    counts = [subtype_counts.get(s, 0) for s in mac_types]
    bars = ax.bar(mac_types, counts, color=colors, edgecolor='black')

    ax.set_xlabel('Mac Subtype', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Unique Regulons', fontsize=12, fontweight='bold')
    ax.set_title(f'Panel C: Mac Subtype-Unique Regulons\n({cutoff_level})',
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

    fig, ax = plt.subplots(figsize=(8, 6))

    subtypes = list(cell_counts.keys())
    counts = list(cell_counts.values())

    colors = ['#2ecc71' if c >= cells_per_cluster else '#e74c3c' for c in counts]

    bars = ax.bar(subtypes, counts, color=colors, edgecolor='black')
    ax.axhline(y=cells_per_cluster, color='blue', linestyle='--', linewidth=2,
               label=f'Sampling size ({cells_per_cluster})')

    ax.set_xlabel('Mac Subtype', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Mac Cell Counts vs Sampling Size', fontsize=14, fontweight='bold')

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
        description='Panel C: Mac subtype-unique regulons (Z-score)'
    )
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--output-dir', type=str, default='results_panel_C_mac_unique',
                        help='Path to output directory')
    parser.add_argument('--cells-per-cluster', type=int, default=200,
                        help='Number of cells to subsample per cluster')
    parser.add_argument('--zscore-cap', type=float, default=2.5,
                        help='Cap z-scores at +/- this value')
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
    logging.info("PANEL C: MAC SUBTYPE-UNIQUE REGULONS (Z-SCORE)")
    logging.info("Compare each Mac subtype with OTHER Mac subtypes")
    logging.info("="*80)
    logging.info(f"Cells per cluster: {CELLS_PER_CLUSTER}")

    # Load data
    logging.info("\n=== Loading AUC Matrices ===")
    mac_data, mac_types, eos_data, eos_clusters = load_auc_matrices(scenic_dir)

    logging.info("\n=== Loading Gene Counts ===")
    gene_counts = load_gene_counts(scenic_dir, mac_types)

    # Calculate mean AUC using subsampling
    logging.info("\n=== Calculating Mean AUC (Subsampling) ===")
    mac_mean_auc, eos_mean_auc, cell_counts, all_tfs = calculate_mean_auc_subsampling(
        mac_data, mac_types, eos_data, eos_clusters, cells_per_cluster=CELLS_PER_CLUSTER
    )

    # Save mean AUC matrix
    mac_mean_auc.to_csv(output_dir / 'mac_mean_auc_subsampling.csv')
    eos_mean_auc.to_csv(output_dir / 'eos_mean_auc_subsampling.csv')

    # Create cell counts plot
    create_cell_counts_plot(cell_counts, CELLS_PER_CLUSTER, output_dir / 'cell_counts_summary.pdf')

    # Process each cutoff level
    logging.info("\n=== Finding Mac Subtype-Unique Regulons ===")

    for cutoff_level in args.cutoff_levels:
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {cutoff_level.upper()}")
        logging.info(f"{'='*60}")

        results_df, mac_unique = find_mac_subtype_unique_regulons(
            mac_mean_auc, eos_mean_auc, mac_types, eos_clusters, cutoff_level=cutoff_level
        )

        # Save results
        results_df.to_csv(output_dir / f'all_regulons_metrics_{cutoff_level}.csv', index=False)

        if len(mac_unique) > 0:
            mac_unique.to_csv(output_dir / f'mac_subtype_unique_{cutoff_level}.csv', index=False)

            create_zscore_heatmap(
                mac_unique, mac_types, eos_clusters,
                output_dir / f'panel_C_heatmap_zscore_{cutoff_level}.pdf',
                gene_counts,
                title=f'Panel C: Mac Subtype-Unique Regulons (n={len(mac_unique)}, {cutoff_level})\nZ-score scaled',
                zscore_cap=args.zscore_cap
            )

            create_summary_barplot(
                mac_unique, mac_types,
                output_dir / f'panel_C_barplot_{cutoff_level}.pdf',
                cutoff_level
            )

        # Print summary
        print(f"\n{'='*60}")
        print(f"PANEL C SUMMARY ({cutoff_level.upper()})")
        print("="*60)
        print(f"Mac subtype-UNIQUE regulons: {len(mac_unique)}")
        if len(mac_unique) > 0:
            print("\nUnique regulons per subtype:")
            subtype_counts = mac_unique['max_subtype'].value_counts()
            for subtype in mac_types:
                print(f"  {subtype}: {subtype_counts.get(subtype, 0)}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)

    print("\n" + "="*80)
    print("SCP COMMAND TO DOWNLOAD RESULTS:")
    print("="*80)
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
