#!/usr/bin/env python3
"""
Combined heatmap: Eos subcluster-unique and Mac subtype-unique regulons
Side by side with single shared colorbar

Heatmap size: 20cm height x 9cm width each
Font: 6pt
Single shared colorbar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Blue to Red colormap
ZSCORE_CMAP = LinearSegmentedColormap.from_list('blue_white_red', ['#2166ac', '#f7f7f7', '#b2182b'])

CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42

CUTOFF_LEVELS = {
    'relax': {
        'min_eos_auc': 0.15,
        'min_mac_auc': 0.15,
        'min_fold_change': 1.2,
        'min_uniqueness_score': 0.03,
    },
}

EXCLUDE_REGULONS = {
    'Erf', 'Pura', 'Zbtb14', 'Rfx5', 'Smad3', 'Jund', 'Tef', 'Bclaf1',
    'Crem', 'Atf4', 'Etv5', 'Zbtb7b', 'Sp3', 'Elk1', 'Cebpb', 'Max',
    'Atf3', 'Ikzf1', 'Atf1', 'Cux1', 'Stat6', 'Elk3'
}


def load_auc_matrices(scenic_dir):
    eos_clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
    eos_data = {}
    for cluster in eos_clusters:
        auc_file = scenic_dir / "eos" / cluster / "step3_auc_matrix.csv"
        if auc_file.exists():
            eos_data[cluster] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {cluster}: {eos_data[cluster].shape[0]} cells")

    mac_types = ['M1', 'M2_Functional', 'M2_Surface']
    mac_data = {}
    for mac_type in mac_types:
        auc_file = scenic_dir / "macrophages" / mac_type / "step3_auc_matrix.csv"
        if auc_file.exists():
            mac_data[mac_type] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {mac_type}: {mac_data[mac_type].shape[0]} cells")

    loaded_eos = [c for c in eos_clusters if c in eos_data]
    loaded_mac = [m for m in mac_types if m in mac_data]
    return eos_data, mac_data, loaded_eos, loaded_mac


def load_gene_counts(scenic_dir, eos_clusters, mac_types):
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
    return gene_counts


def extract_tf_name(regulon_name):
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_mean_auc_subsampling(eos_data, mac_data, eos_clusters, mac_types, cells_per_cluster=200):
    np.random.seed(RANDOM_SEED)

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
        for tf in all_tfs:
            if tf in eos_tf_to_regulon.get(cluster, {}):
                regulon_col = eos_tf_to_regulon[cluster][tf]
                eos_auc[cluster][tf] = sampled_df[regulon_col].mean()
            else:
                eos_auc[cluster][tf] = 0.0

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
        for tf in all_tfs:
            if tf in mac_tf_to_regulon.get(mac_type, {}):
                regulon_col = mac_tf_to_regulon[mac_type][tf]
                mac_auc[mac_type][tf] = sampled_df[regulon_col].mean()
            else:
                mac_auc[mac_type][tf] = 0.0

    eos_mean_auc = pd.DataFrame(eos_auc).reindex(all_tfs)
    mac_mean_auc = pd.DataFrame(mac_auc).reindex(all_tfs)
    return eos_mean_auc, mac_mean_auc, all_tfs


def find_eos_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types):
    criteria = CUTOFF_LEVELS['relax']
    all_tfs = list(eos_mean_auc.index)
    results = []

    for tf in all_tfs:
        eos_aucs = {c: eos_mean_auc.loc[tf, c] for c in eos_clusters if c in eos_mean_auc.columns}
        mac_aucs = {m: mac_mean_auc.loc[tf, m] for m in mac_types if m in mac_mean_auc.columns}
        if not eos_aucs or not mac_aucs:
            continue

        max_eos_cluster = max(eos_aucs, key=eos_aucs.get)
        max_eos_auc = eos_aucs[max_eos_cluster]
        mac_mean = np.mean(list(mac_aucs.values()))
        fold_change = max_eos_auc / mac_mean if mac_mean > 0 else 0
        sorted_eos = sorted(eos_aucs.values(), reverse=True)
        uniqueness = sorted_eos[0] - sorted_eos[1] if len(sorted_eos) >= 2 else sorted_eos[0]

        results.append({
            'tf': tf, 'max_eos_cluster': max_eos_cluster, 'max_eos_auc': max_eos_auc,
            'fold_change': fold_change, 'uniqueness_score': uniqueness,
            **{f'eos_{c}': eos_aucs.get(c, 0) for c in eos_clusters},
            **{f'mac_{m}': mac_aucs.get(m, 0) for m in mac_types}
        })

    df = pd.DataFrame(results)
    eos_unique = df[
        (df['max_eos_auc'] >= criteria['min_eos_auc']) &
        (df['fold_change'] >= criteria['min_fold_change']) &
        (df['uniqueness_score'] >= criteria['min_uniqueness_score']) &
        (~df['tf'].isin(EXCLUDE_REGULONS))
    ].sort_values('uniqueness_score', ascending=False)
    return eos_unique


def find_mac_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types):
    criteria = CUTOFF_LEVELS['relax']
    all_tfs = list(mac_mean_auc.index)
    results = []

    for tf in all_tfs:
        mac_aucs = {m: mac_mean_auc.loc[tf, m] for m in mac_types if m in mac_mean_auc.columns}
        eos_aucs = {c: eos_mean_auc.loc[tf, c] for c in eos_clusters if c in eos_mean_auc.columns}
        if not mac_aucs or not eos_aucs:
            continue

        max_mac_subtype = max(mac_aucs, key=mac_aucs.get)
        max_mac_auc = mac_aucs[max_mac_subtype]
        eos_mean = np.mean(list(eos_aucs.values()))
        fold_change = max_mac_auc / eos_mean if eos_mean > 0 else 0
        sorted_mac = sorted(mac_aucs.values(), reverse=True)
        uniqueness = sorted_mac[0] - sorted_mac[1] if len(sorted_mac) >= 2 else sorted_mac[0]

        results.append({
            'tf': tf, 'max_mac_subtype': max_mac_subtype, 'max_mac_auc': max_mac_auc,
            'fold_change': fold_change, 'uniqueness_score': uniqueness,
            **{f'mac_{m}': mac_aucs.get(m, 0) for m in mac_types},
            **{f'eos_{c}': eos_aucs.get(c, 0) for c in eos_clusters}
        })

    df = pd.DataFrame(results)
    mac_unique = df[
        (df['max_mac_auc'] >= criteria['min_mac_auc']) &
        (df['fold_change'] >= criteria['min_fold_change']) &
        (df['uniqueness_score'] >= criteria['min_uniqueness_score'])
    ].sort_values('uniqueness_score', ascending=False)

    # Top 8 per subtype
    top_per_subtype = []
    for subtype in mac_types:
        subtype_df = mac_unique[mac_unique['max_mac_subtype'] == subtype].head(8)
        top_per_subtype.append(subtype_df)
    mac_unique = pd.concat(top_per_subtype, ignore_index=True)
    return mac_unique


def create_combined_heatmap(eos_df, mac_df, eos_clusters, mac_types, output_path, gene_counts, zscore_cap=2.5):
    """Create side-by-side heatmaps with single shared colorbar"""

    # Order Eos by cluster
    ordered_eos = []
    for cluster in eos_clusters:
        cluster_df = eos_df[eos_df['max_eos_cluster'] == cluster].copy()
        if len(cluster_df) > 0:
            cluster_df = cluster_df.sort_values('uniqueness_score', ascending=False)
            ordered_eos.append(cluster_df)
    if ordered_eos:
        eos_df = pd.concat(ordered_eos, ignore_index=True)

    # Order Mac by subtype
    ordered_mac = []
    for subtype in mac_types:
        subtype_df = mac_df[mac_df['max_mac_subtype'] == subtype].copy()
        if len(subtype_df) > 0:
            subtype_df = subtype_df.sort_values('uniqueness_score', ascending=False)
            ordered_mac.append(subtype_df)
    if ordered_mac:
        mac_df = pd.concat(ordered_mac, ignore_index=True)

    # Prepare Eos data (Eos clusters first, then Mac)
    eos_cols = eos_clusters + mac_types
    eos_data_cols = [f'eos_{c}' if c in eos_clusters else f'mac_{c}' for c in eos_cols]
    eos_raw = eos_df[eos_data_cols].values
    eos_means = eos_raw.mean(axis=1, keepdims=True)
    eos_stds = eos_raw.std(axis=1, keepdims=True)
    eos_stds[eos_stds == 0] = 1
    eos_zscore = np.clip((eos_raw - eos_means) / eos_stds, -zscore_cap, zscore_cap)

    # Prepare Mac data (Mac first, then Eos clusters)
    mac_cols = mac_types + eos_clusters
    mac_data_cols = [f'mac_{m}' if m in mac_types else f'eos_{m}' for m in mac_cols]
    mac_raw = mac_df[mac_data_cols].values
    mac_means = mac_raw.mean(axis=1, keepdims=True)
    mac_stds = mac_raw.std(axis=1, keepdims=True)
    mac_stds[mac_stds == 0] = 1
    mac_zscore = np.clip((mac_raw - mac_means) / mac_stds, -zscore_cap, zscore_cap)

    # Y-axis labels
    eos_ylabels = [f"{row['tf']}({gene_counts.get(row['tf'], '')})" if row['tf'] in gene_counts else row['tf']
                   for _, row in eos_df.iterrows()]
    mac_ylabels = [f"{row['tf']}({gene_counts.get(row['tf'], '')})" if row['tf'] in gene_counts else row['tf']
                   for _, row in mac_df.iterrows()]

    # Figure dimensions: 20cm height x 9cm width per heatmap
    heatmap_width_cm = 9
    heatmap_height_cm = 20
    heatmap_width = heatmap_width_cm / 2.54
    heatmap_height = heatmap_height_cm / 2.54

    # Total figure size - increased middle_gap to prevent label overlap
    left_margin = 1.0
    middle_gap = 1.5  # Increased from 0.3 to 1.5 inches for y-axis labels of right plot
    right_margin = 0.6  # for colorbar
    top_margin = 0.4
    bottom_margin = 0.8

    fig_width = left_margin + heatmap_width + middle_gap + heatmap_width + right_margin
    fig_height = top_margin + heatmap_height + bottom_margin

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Calculate positions
    ax1_left = left_margin / fig_width
    ax1_bottom = bottom_margin / fig_height
    ax1_width = heatmap_width / fig_width
    ax1_height = heatmap_height / fig_height

    ax2_left = (left_margin + heatmap_width + middle_gap) / fig_width
    ax2_bottom = bottom_margin / fig_height
    ax2_width = heatmap_width / fig_width
    ax2_height = heatmap_height / fig_height

    # Left heatmap: Eos unique
    ax1 = fig.add_axes([ax1_left, ax1_bottom, ax1_width, ax1_height])
    im1 = ax1.imshow(eos_zscore, cmap=ZSCORE_CMAP, aspect='auto', vmin=-zscore_cap, vmax=zscore_cap)
    ax1.set_xticks(np.arange(len(eos_cols)))
    ax1.set_xticklabels(eos_cols, rotation=45, ha='right', fontsize=6)
    ax1.set_yticks(np.arange(len(eos_ylabels)))
    ax1.set_yticklabels(eos_ylabels, fontsize=6)
    ax1.axvline(x=len(eos_clusters) - 0.5, color='black', linewidth=1)
    ax1.set_title(f'Eos Subcluster-Unique (n={len(eos_df)})', fontsize=6, pad=5)
    ax1.set_xlabel('Cell Type', fontsize=6)
    ax1.set_ylabel('Regulon', fontsize=6)

    # Right heatmap: Mac unique
    ax2 = fig.add_axes([ax2_left, ax2_bottom, ax2_width, ax2_height])
    im2 = ax2.imshow(mac_zscore, cmap=ZSCORE_CMAP, aspect='auto', vmin=-zscore_cap, vmax=zscore_cap)
    ax2.set_xticks(np.arange(len(mac_cols)))
    ax2.set_xticklabels(mac_cols, rotation=45, ha='right', fontsize=6)
    ax2.set_yticks(np.arange(len(mac_ylabels)))
    ax2.set_yticklabels(mac_ylabels, fontsize=6)
    ax2.axvline(x=len(mac_types) - 0.5, color='black', linewidth=1)
    ax2.set_title(f'Mac Subtype-Unique (n={len(mac_df)})', fontsize=6, pad=5)
    ax2.set_xlabel('Cell Type', fontsize=6)
    ax2.set_ylabel('Regulon', fontsize=6)

    # Single shared colorbar (positioned relative to figure, not axes)
    cbar_left = ax2_left + ax2_width + 0.01
    cbar_bottom = ax2_bottom + ax2_height * 0.35
    cbar_width = 0.015
    cbar_height = ax2_height * 0.3
    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
    cbar = plt.colorbar(im2, cax=cbar_ax)
    cbar.set_label('Z-score', fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    plt.savefig(output_path, dpi=300)
    logging.info(f"Saved: {output_path}")
    plt.close()

    # Also save as PNG
    png_path = str(output_path).replace('.pdf', '.png')
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = fig.add_axes([ax1_left, ax1_bottom, ax1_width, ax1_height])
    im1 = ax1.imshow(eos_zscore, cmap=ZSCORE_CMAP, aspect='auto', vmin=-zscore_cap, vmax=zscore_cap)
    ax1.set_xticks(np.arange(len(eos_cols)))
    ax1.set_xticklabels(eos_cols, rotation=45, ha='right', fontsize=6)
    ax1.set_yticks(np.arange(len(eos_ylabels)))
    ax1.set_yticklabels(eos_ylabels, fontsize=6)
    ax1.axvline(x=len(eos_clusters) - 0.5, color='black', linewidth=1)
    ax1.set_title(f'Eos Subcluster-Unique (n={len(eos_df)})', fontsize=6, pad=5)
    ax1.set_xlabel('Cell Type', fontsize=6)
    ax1.set_ylabel('Regulon', fontsize=6)
    ax2 = fig.add_axes([ax2_left, ax2_bottom, ax2_width, ax2_height])
    im2 = ax2.imshow(mac_zscore, cmap=ZSCORE_CMAP, aspect='auto', vmin=-zscore_cap, vmax=zscore_cap)
    ax2.set_xticks(np.arange(len(mac_cols)))
    ax2.set_xticklabels(mac_cols, rotation=45, ha='right', fontsize=6)
    ax2.set_yticks(np.arange(len(mac_ylabels)))
    ax2.set_yticklabels(mac_ylabels, fontsize=6)
    ax2.axvline(x=len(mac_types) - 0.5, color='black', linewidth=1)
    ax2.set_title(f'Mac Subtype-Unique (n={len(mac_df)})', fontsize=6, pad=5)
    ax2.set_xlabel('Cell Type', fontsize=6)
    ax2.set_ylabel('Regulon', fontsize=6)
    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
    cbar = plt.colorbar(im2, cax=cbar_ax)
    cbar.set_label('Z-score', fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    plt.savefig(png_path, dpi=300)
    logging.info(f"Saved: {png_path}")
    plt.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Combined Eos and Mac unique regulon heatmaps')
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned')
    parser.add_argument('--output-dir', type=str, default='results_combined_heatmap')
    parser.add_argument('--zscore-cap', type=float, default=2.5)
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading data...")
    eos_data, mac_data, eos_clusters, mac_types = load_auc_matrices(scenic_dir)
    gene_counts = load_gene_counts(scenic_dir, eos_clusters, mac_types)

    logging.info("Calculating mean AUC...")
    eos_mean_auc, mac_mean_auc, all_tfs = calculate_mean_auc_subsampling(
        eos_data, mac_data, eos_clusters, mac_types, CELLS_PER_CLUSTER
    )

    logging.info("Finding unique regulons...")
    eos_unique = find_eos_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types)
    mac_unique = find_mac_unique_regulons(eos_mean_auc, mac_mean_auc, eos_clusters, mac_types)

    logging.info(f"Eos unique: {len(eos_unique)}, Mac unique: {len(mac_unique)}")

    # Save CSVs
    eos_unique.to_csv(output_dir / 'eos_subcluster_unique_relax.csv', index=False)
    mac_unique.to_csv(output_dir / 'mac_subtype_unique_relax.csv', index=False)

    # Create combined figure
    create_combined_heatmap(
        eos_unique, mac_unique, eos_clusters, mac_types,
        output_dir / 'combined_eos_mac_unique_heatmap_relax.pdf',
        gene_counts, args.zscore_cap
    )

    print(f"\nOutput saved to: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
