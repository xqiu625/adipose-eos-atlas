#!/usr/bin/env python3
"""
Find Mac subtype-UNIQUE regulons compared to ALL Eosinophils using Random Subsampling
Plot heatmaps using Z-SCORE scaled mean AUC with blue-to-red colormap

Updated 2024-12-08:
- Heatmap size: 20cm (height) x 9cm (width)
- Font size: 6pt for all text
- Single smaller colorbar

This is the REVERSE of 1208_Eos_vs_Mac_Subcluster_Unique_Zscore.py:
- Previous: Eos subcluster-specific regulons vs Mac
- This: Mac subtype-specific regulons vs Eos

This script finds regulons that are:
1. ENRICHED in Mac compared to Eos (Mac subtype AUC >> Eos AUC)
2. UNIQUE to specific Mac subtypes (high in ONE Mac subtype, low in OTHER Mac subtypes)

Output:
- Heatmaps with Z-score scaled values (row-wise z-score)
- Blue (low) to Red (high) color scale
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

# Blue to Red colormap (blue=low, white=middle, red=high)
ZSCORE_CMAP = LinearSegmentedColormap.from_list('blue_white_red', ['#2166ac', '#f7f7f7', '#b2182b'])

# Subsampling parameters
CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42

# Cutoff thresholds
CUTOFF_LEVELS = {
    'relax': {
        'min_mac_auc': 0.15,
        'min_fold_change_vs_eos': 1.2,
        'min_uniqueness_score': 0.03,
    },
    'moderate': {
        'min_mac_auc': 0.20,
        'min_fold_change_vs_eos': 1.5,
        'min_uniqueness_score': 0.05,
    },
    'stringent': {
        'min_mac_auc': 0.25,
        'min_fold_change_vs_eos': 2.0,
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


def find_mac_enriched_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types,
                                       cutoff_level='relax'):
    """
    Find regulons that are:
    1. ENRICHED in Mac vs Eos (max Mac subtype >> Eos mean)
    2. UNIQUE to specific Mac subtypes (high in ONE subtype, low in others)
    """

    criteria = CUTOFF_LEVELS[cutoff_level]
    min_mac_auc = criteria['min_mac_auc']
    min_fc_vs_eos = criteria['min_fold_change_vs_eos']
    min_uniqueness = criteria['min_uniqueness_score']

    all_tfs = list(mac_mean_auc.index)
    results = []

    for tf in all_tfs:
        # Get Mac subtype AUCs
        mac_aucs = {mac_type: mac_mean_auc.loc[tf, mac_type] for mac_type in mac_types
                    if mac_type in mac_mean_auc.columns}

        # Get Eos AUCs
        eos_aucs = {cluster: eos_mean_auc.loc[tf, cluster] for cluster in eos_clusters
                    if cluster in eos_mean_auc.columns}

        if not mac_aucs or not eos_aucs:
            continue

        # Find max Mac subtype
        max_mac_subtype = max(mac_aucs, key=mac_aucs.get)
        max_mac_auc = mac_aucs[max_mac_subtype]

        # Calculate Mac mean (across subtypes)
        mac_mean = np.mean(list(mac_aucs.values()))

        # Calculate Eos mean (across all subclusters)
        eos_mean = np.mean(list(eos_aucs.values()))

        # Calculate fold change vs Eos
        fold_change_vs_eos = max_mac_auc / eos_mean if eos_mean > 0 else 0

        # Calculate uniqueness score (max Mac - second max Mac)
        sorted_mac_aucs = sorted(mac_aucs.values(), reverse=True)
        if len(sorted_mac_aucs) >= 2:
            uniqueness_score = sorted_mac_aucs[0] - sorted_mac_aucs[1]
        else:
            uniqueness_score = sorted_mac_aucs[0]

        # Other Mac subtypes
        other_mac_aucs = [v for k, v in mac_aucs.items() if k != max_mac_subtype]
        max_other_mac = max(other_mac_aucs) if other_mac_aucs else 0
        mean_other_mac = np.mean(other_mac_aucs) if other_mac_aucs else 0

        results.append({
            'tf': tf,
            'max_mac_subtype': max_mac_subtype,
            'max_mac_auc': max_mac_auc,
            'mac_mean': mac_mean,
            'eos_mean': eos_mean,
            'fold_change_vs_eos': fold_change_vs_eos,
            'uniqueness_score': uniqueness_score,
            'max_other_mac_auc': max_other_mac,
            'mean_other_mac_auc': mean_other_mac,
            **{f'mac_{m}': mac_aucs.get(m, 0) for m in mac_types},
            **{f'eos_{c}': eos_aucs.get(c, 0) for c in eos_clusters}
        })

    results_df = pd.DataFrame(results)

    # Filter for Mac-enriched regulons
    mac_enriched = results_df[
        (results_df['max_mac_auc'] >= min_mac_auc) &
        (results_df['fold_change_vs_eos'] >= min_fc_vs_eos)
    ].copy()

    # Filter for Mac subtype-UNIQUE regulons (enriched vs Eos AND unique within Mac)
    mac_unique = results_df[
        (results_df['max_mac_auc'] >= min_mac_auc) &
        (results_df['fold_change_vs_eos'] >= min_fc_vs_eos) &
        (results_df['uniqueness_score'] >= min_uniqueness)
    ].copy()

    # Sort by uniqueness score
    mac_enriched = mac_enriched.sort_values('fold_change_vs_eos', ascending=False)
    mac_unique = mac_unique.sort_values('uniqueness_score', ascending=False)

    # Keep only top 8 regulons per Mac subtype
    top_per_subtype = []
    for subtype in mac_types:
        subtype_df = mac_unique[mac_unique['max_mac_subtype'] == subtype].head(8)
        top_per_subtype.append(subtype_df)
    mac_unique = pd.concat(top_per_subtype, ignore_index=True)

    logging.info(f"\n{cutoff_level.upper()} cutoff:")
    logging.info(f"  Mac-enriched vs Eos: {len(mac_enriched)} regulons")
    logging.info(f"  Mac subtype-UNIQUE (top 8 per subtype): {len(mac_unique)} regulons")

    # Summary by subtype
    if len(mac_unique) > 0:
        subtype_counts = mac_unique['max_mac_subtype'].value_counts()
        for subtype in mac_types:
            count = subtype_counts.get(subtype, 0)
            logging.info(f"    {subtype}: {count} unique regulons")

    return results_df, mac_enriched, mac_unique


def create_zscore_heatmap(regulon_df, mac_types, eos_clusters, output_path, gene_counts,
                           title="Mac-Enriched Regulons (Z-score)", group_by_subtype=True,
                           zscore_cap=2.5):
    """Create heatmap with Z-score scaled values and blue-to-red colormap

    Updated sizing:
    - Heatmap: 20cm height x 9cm width
    - Font: 6pt for all text
    - Single smaller colorbar

    Note: Mac subtypes are shown FIRST (left), then Eos subclusters (right)
    """

    if len(regulon_df) == 0:
        logging.warning(f"No regulons to plot for {title}")
        return

    # Order by max Mac subtype
    if group_by_subtype:
        ordered_dfs = []
        for subtype in mac_types:
            subtype_df = regulon_df[regulon_df['max_mac_subtype'] == subtype].copy()
            if len(subtype_df) > 0:
                subtype_df = subtype_df.sort_values('uniqueness_score', ascending=False)
                ordered_dfs.append(subtype_df)
        if ordered_dfs:
            regulon_df = pd.concat(ordered_dfs, ignore_index=True)

    # Build plot data: Mac first, then Eos
    all_cols = mac_types + eos_clusters
    data_cols = [f'mac_{m}' if m in mac_types else f'eos_{m}' for m in all_cols]

    raw_data = regulon_df[data_cols].values

    # Calculate row-wise z-score
    row_means = raw_data.mean(axis=1, keepdims=True)
    row_stds = raw_data.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    zscore_data = (raw_data - row_means) / row_stds

    # Cap z-scores for better visualization
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

    # ========================================
    # NEW SIZING: 20cm height x 9cm width for heatmap
    # ========================================
    heatmap_width_cm = 9
    heatmap_height_cm = 20

    # Convert cm to inches (1 inch = 2.54 cm)
    heatmap_width_inches = heatmap_width_cm / 2.54
    heatmap_height_inches = heatmap_height_cm / 2.54

    # Add margins for labels and colorbar
    left_margin = 1.0   # inches for y-axis labels
    right_margin = 0.6  # inches for colorbar
    top_margin = 0.5    # inches for title
    bottom_margin = 0.8 # inches for x-axis labels

    fig_width = left_margin + heatmap_width_inches + right_margin
    fig_height = top_margin + heatmap_height_inches + bottom_margin

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set axis position to exact location for correct cell sizes
    ax_left = left_margin / fig_width
    ax_bottom = bottom_margin / fig_height
    ax_width = heatmap_width_inches / fig_width
    ax_height = heatmap_height_inches / fig_height
    ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

    # Plot with z-score and blue-to-red colormap
    im = ax.imshow(zscore_data_capped, cmap=ZSCORE_CMAP, aspect='auto',
                   vmin=-zscore_cap, vmax=zscore_cap)

    # Add single smaller colorbar
    cbar_width = 0.02  # narrower
    cbar_height = 0.3  # shorter (30% of heatmap height)
    cbar_left = ax_left + ax_width + 0.02
    cbar_bottom = ax_bottom + ax_height * 0.35  # centered vertically
    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height * ax_height / 0.3])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('Z-score', fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    # 6pt font for all text
    ax.set_xticks(np.arange(len(all_cols)))
    ax.set_xticklabels(all_cols, rotation=45, ha='right', fontsize=6)

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=6)

    # Add vertical line to separate Mac and Eos
    ax.axvline(x=len(mac_types) - 0.5, color='black', linewidth=1)

    ax.set_xlabel('Cell Type', fontsize=6)
    ax.set_ylabel('Regulon (gene count)', fontsize=6)
    ax.set_title(title, fontsize=6, pad=5)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Find Mac subtype-unique regulons vs Eos using subsampling (Z-score heatmaps)'
    )
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--output-dir', type=str, default='results_mac_vs_eos_zscore',
                        help='Path to output directory')
    parser.add_argument('--cells-per-cluster', type=int, default=200,
                        help='Number of cells to subsample per cluster')
    parser.add_argument('--zscore-cap', type=float, default=2.5,
                        help='Cap z-scores at +/- this value for visualization')
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
    logging.info("MAC SUBTYPE-UNIQUE REGULONS vs EOS (Z-SCORE)")
    logging.info("Method: Random Subsampling + Z-score scaling")
    logging.info("="*80)
    logging.info(f"Cells per cluster: {CELLS_PER_CLUSTER}")
    logging.info(f"Z-score cap: +/- {args.zscore_cap}")

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

    # Process each cutoff level
    logging.info("\n=== Finding Mac-Enriched and Unique Regulons ===")

    for cutoff_level in args.cutoff_levels:
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {cutoff_level.upper()}")
        logging.info(f"{'='*60}")

        results_df, mac_enriched, mac_unique = find_mac_enriched_unique_regulons(
            eos_mean_auc, mac_mean_auc, eos_clusters, mac_types,
            cutoff_level=cutoff_level
        )

        # Output heatmap PDF, CSV matrix, and regulon info for mac_unique
        if len(mac_unique) > 0:
            # Save regulon info CSV (needed for GRN network analysis)
            mac_unique.to_csv(output_dir / f'mac_subtype_unique_{cutoff_level}.csv', index=False)

            create_zscore_heatmap(
                mac_unique, mac_types, eos_clusters,
                output_dir / f'mac_subtype_unique_heatmap_{cutoff_level}.pdf',
                gene_counts,
                title=f'Mac Subtype-UNIQUE Regulons (n={len(mac_unique)}, {cutoff_level})\nZ-score scaled',
                zscore_cap=args.zscore_cap
            )

        # Print summary
        print(f"\n{'='*60}")
        print(f"SUMMARY ({cutoff_level.upper()})")
        print("="*60)
        print(f"Mac subtype-UNIQUE: {len(mac_unique)} regulons")
        if len(mac_unique) > 0:
            print("\nUnique regulons per subtype:")
            subtype_counts = mac_unique['max_mac_subtype'].value_counts()
            for subtype in mac_types:
                print(f"  {subtype}: {subtype_counts.get(subtype, 0)}")

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
