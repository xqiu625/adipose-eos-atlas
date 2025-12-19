#!/usr/bin/env python3
"""
Panel A: Regulon activity between ALL Eos and ALL Macs
Plot Fold Change heatmap comparing pooled Eos vs pooled Mac

This combines all Eos subclusters into "Eos_all" and all Mac subtypes into "Mac_all"
to show regulons differentially active between the two cell types.
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

# Blue to Red colormap for fold change
FC_CMAP = LinearSegmentedColormap.from_list('blue_white_red', ['#2166ac', '#f7f7f7', '#b2182b'])

# Subsampling parameters
CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42

# Cutoff thresholds for differential regulons
CUTOFF_LEVELS = {
    'relax': {
        'min_auc': 0.15,
        'min_fold_change': 1.2,
        'min_diff': 0.05,
    },
    'moderate': {
        'min_auc': 0.20,
        'min_fold_change': 1.5,
        'min_diff': 0.08,
    },
    'stringent': {
        'min_auc': 0.25,
        'min_fold_change': 2.0,
        'min_diff': 0.10,
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


def calculate_pooled_mean_auc_subsampling(eos_data, mac_data, eos_clusters, mac_types,
                                           cells_per_cluster=200):
    """
    Calculate mean AUC for pooled Eos_all and Mac_all using subsampling.
    Also calculate per-subcluster means for reference.
    """
    logging.info(f"\n=== Random Subsampling ===")
    logging.info(f"Cells per cluster: {cells_per_cluster}")

    np.random.seed(RANDOM_SEED)

    # Build TF mapping
    eos_tf_to_regulon = {}
    mac_tf_to_regulon = {}
    all_tfs = set()

    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        eos_tf_to_regulon[cluster] = {}
        for col in eos_data[cluster].columns:
            tf_name = extract_tf_name(col)
            eos_tf_to_regulon[cluster][tf_name] = col
            all_tfs.add(tf_name)

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
    for cluster in eos_clusters:
        if cluster in eos_data:
            cell_counts[cluster] = eos_data[cluster].shape[0]
    for mac_type in mac_types:
        if mac_type in mac_data:
            cell_counts[mac_type] = mac_data[mac_type].shape[0]

    # Subsample and calculate per-cluster means
    eos_cluster_auc = {cluster: {} for cluster in eos_clusters}
    mac_subtype_auc = {mac_type: {} for mac_type in mac_types}

    # Collect all sampled cells for pooled calculation
    eos_all_sampled = []
    mac_all_sampled = []

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
        eos_all_sampled.append((cluster, sampled_df))
        logging.info(f"  Subsampled {cluster}: {len(sampled_idx)} cells")

        for tf in all_tfs:
            if tf in eos_tf_to_regulon.get(cluster, {}):
                regulon_col = eos_tf_to_regulon[cluster][tf]
                eos_cluster_auc[cluster][tf] = sampled_df[regulon_col].mean()
            else:
                eos_cluster_auc[cluster][tf] = 0.0

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
        mac_all_sampled.append((mac_type, sampled_df))
        logging.info(f"  Subsampled {mac_type}: {len(sampled_idx)} cells")

        for tf in all_tfs:
            if tf in mac_tf_to_regulon.get(mac_type, {}):
                regulon_col = mac_tf_to_regulon[mac_type][tf]
                mac_subtype_auc[mac_type][tf] = sampled_df[regulon_col].mean()
            else:
                mac_subtype_auc[mac_type][tf] = 0.0

    # Calculate pooled Eos_all and Mac_all means (average across subclusters/subtypes)
    eos_all_auc = {}
    mac_all_auc = {}

    for tf in all_tfs:
        eos_values = [eos_cluster_auc[c][tf] for c in eos_clusters if c in eos_cluster_auc]
        mac_values = [mac_subtype_auc[m][tf] for m in mac_types if m in mac_subtype_auc]
        eos_all_auc[tf] = np.mean(eos_values) if eos_values else 0.0
        mac_all_auc[tf] = np.mean(mac_values) if mac_values else 0.0

    # Create DataFrames
    eos_cluster_df = pd.DataFrame(eos_cluster_auc).reindex(all_tfs)
    mac_subtype_df = pd.DataFrame(mac_subtype_auc).reindex(all_tfs)

    pooled_df = pd.DataFrame({
        'Eos_all': eos_all_auc,
        'Mac_all': mac_all_auc
    }).reindex(all_tfs)

    logging.info(f"Subsampling complete")

    return eos_cluster_df, mac_subtype_df, pooled_df, cell_counts, all_tfs


def find_differential_regulons(pooled_df, eos_cluster_df, mac_subtype_df,
                                eos_clusters, mac_types, cutoff_level='relax'):
    """
    Find regulons differentially active between Eos and Mac.
    """
    criteria = CUTOFF_LEVELS[cutoff_level]
    min_auc = criteria['min_auc']
    min_fc = criteria['min_fold_change']
    min_diff = criteria['min_diff']

    all_tfs = list(pooled_df.index)
    results = []

    for tf in all_tfs:
        eos_auc = pooled_df.loc[tf, 'Eos_all']
        mac_auc = pooled_df.loc[tf, 'Mac_all']

        # Calculate fold change and difference
        fc_eos_vs_mac = eos_auc / mac_auc if mac_auc > 0 else 0
        fc_mac_vs_eos = mac_auc / eos_auc if eos_auc > 0 else 0
        diff = eos_auc - mac_auc

        # Determine which cell type has higher expression
        max_auc = max(eos_auc, mac_auc)
        enriched_in = 'Eos' if eos_auc > mac_auc else 'Mac'

        # Get per-cluster values
        eos_cluster_values = {c: eos_cluster_df.loc[tf, c] for c in eos_clusters if c in eos_cluster_df.columns}
        mac_subtype_values = {m: mac_subtype_df.loc[tf, m] for m in mac_types if m in mac_subtype_df.columns}

        results.append({
            'tf': tf,
            'eos_all_auc': eos_auc,
            'mac_all_auc': mac_auc,
            'diff_eos_minus_mac': diff,
            'fc_eos_vs_mac': fc_eos_vs_mac,
            'fc_mac_vs_eos': fc_mac_vs_eos,
            'max_auc': max_auc,
            'enriched_in': enriched_in,
            **{f'eos_{c}': eos_cluster_values.get(c, 0) for c in eos_clusters},
            **{f'mac_{m}': mac_subtype_values.get(m, 0) for m in mac_types}
        })

    results_df = pd.DataFrame(results)

    # Filter for Eos-enriched regulons
    eos_enriched = results_df[
        (results_df['eos_all_auc'] >= min_auc) &
        (results_df['fc_eos_vs_mac'] >= min_fc) &
        (results_df['diff_eos_minus_mac'] >= min_diff)
    ].copy().sort_values('fc_eos_vs_mac', ascending=False)

    # Filter for Mac-enriched regulons
    mac_enriched = results_df[
        (results_df['mac_all_auc'] >= min_auc) &
        (results_df['fc_mac_vs_eos'] >= min_fc) &
        (results_df['diff_eos_minus_mac'] <= -min_diff)
    ].copy().sort_values('fc_mac_vs_eos', ascending=False)

    logging.info(f"\n{cutoff_level.upper()} cutoff:")
    logging.info(f"  Eos-enriched regulons: {len(eos_enriched)}")
    logging.info(f"  Mac-enriched regulons: {len(mac_enriched)}")

    return results_df, eos_enriched, mac_enriched


def create_foldchange_heatmap_panel_a(regulon_df, eos_clusters, mac_types, output_path, gene_counts,
                                       title="Eos vs Mac Regulons (Fold Change)", fc_cap=5.0):
    """
    Create Panel A heatmap: Two columns showing log2(Eos/Mac) and log2(Mac/Eos)
    """
    if len(regulon_df) == 0:
        logging.warning(f"No regulons to plot for {title}")
        return

    # Sort by enrichment direction and fold change
    eos_enriched = regulon_df[regulon_df['enriched_in'] == 'Eos'].sort_values('fc_eos_vs_mac', ascending=False)
    mac_enriched = regulon_df[regulon_df['enriched_in'] == 'Mac'].sort_values('fc_mac_vs_eos', ascending=False)
    regulon_df = pd.concat([eos_enriched, mac_enriched], ignore_index=True)

    # Calculate log2 fold changes for both directions
    fc_eos_mac = []
    fc_mac_eos = []
    pseudo = 0.001

    for _, row in regulon_df.iterrows():
        eos_auc = row['eos_all_auc']
        mac_auc = row['mac_all_auc']
        # log2(Eos/Mac)
        log2_fc_eos = np.log2((eos_auc + pseudo) / (mac_auc + pseudo))
        # log2(Mac/Eos)
        log2_fc_mac = np.log2((mac_auc + pseudo) / (eos_auc + pseudo))
        fc_eos_mac.append(log2_fc_eos)
        fc_mac_eos.append(log2_fc_mac)

    # Create 2-column matrix: [Eos/Mac, Mac/Eos]
    fc_data = np.column_stack([fc_eos_mac, fc_mac_eos])
    fc_data_capped = np.clip(fc_data, -fc_cap, fc_cap)

    # Y-axis labels
    ylabels = []
    for _, row in regulon_df.iterrows():
        tf = row['tf']
        if tf in gene_counts:
            ylabels.append(f"{tf}({gene_counts[tf]})")
        else:
            ylabels.append(tf)

    # Save fold change matrix as CSV
    fc_df = pd.DataFrame({
        'TF': regulon_df['tf'].values,
        'log2_FC_Eos_vs_Mac': fc_eos_mac,
        'log2_FC_Mac_vs_Eos': fc_mac_eos,
        'Eos_all_AUC': regulon_df['eos_all_auc'].values,
        'Mac_all_AUC': regulon_df['mac_all_auc'].values,
        'enriched_in': regulon_df['enriched_in'].values
    })
    csv_path = str(output_path).replace('.pdf', '_matrix.csv')
    fc_df.to_csv(csv_path, index=False)
    logging.info(f"Saved fold change matrix: {csv_path}")

    # Create figure - width for 2 columns
    fig_height = max(8, len(regulon_df) * 0.35)
    fig_width = 6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    im = ax.imshow(fc_data_capped, cmap=FC_CMAP, aspect='auto',
                   vmin=-fc_cap, vmax=fc_cap)

    cbar = plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label('log2(Fold Change)', fontsize=12, fontweight='bold')

    # X-axis labels for two columns
    col_labels = ['log2(Eos/Mac)', 'log2(Mac/Eos)']
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=12)

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=8)

    # Add vertical line to separate the two columns
    ax.axvline(x=0.5, color='black', linewidth=2)

    # Add horizontal line between Eos-enriched and Mac-enriched
    n_eos_enriched = len(eos_enriched)
    if n_eos_enriched > 0 and n_eos_enriched < len(regulon_df):
        ax.axhline(y=n_eos_enriched - 0.5, color='black', linewidth=1, linestyle='--')

    ax.set_ylabel('Regulon (gene count)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_summary_barplot(eos_enriched, mac_enriched, output_path, cutoff_level):
    """Create bar plot showing number of enriched regulons"""

    fig, ax = plt.subplots(figsize=(8, 6))

    categories = ['Eos-enriched', 'Mac-enriched']
    counts = [len(eos_enriched), len(mac_enriched)]
    colors = ['#e41a1c', '#377eb8']

    bars = ax.bar(categories, counts, color=colors, edgecolor='black')

    ax.set_ylabel('Number of Regulons', fontsize=12, fontweight='bold')
    ax.set_title(f'Differential Regulons: Eos vs Mac\n({cutoff_level})',
                 fontsize=14, fontweight='bold')

    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=12, fontweight='bold')

    max_val = max(counts) if counts else 1
    ax.set_ylim(0, max_val * 1.15)

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
        description='Panel A: Regulon activity between ALL Eos and ALL Macs (Fold Change)'
    )
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--output-dir', type=str, default='results_panel_A_eos_vs_mac',
                        help='Path to output directory')
    parser.add_argument('--cells-per-cluster', type=int, default=200,
                        help='Number of cells to subsample per cluster')
    parser.add_argument('--fc-cap', type=float, default=5.0,
                        help='Cap log2 fold change at +/- this value')
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
    logging.info("PANEL A: EOS vs MAC REGULON ACTIVITY (FOLD CHANGE)")
    logging.info("="*80)
    logging.info(f"Cells per cluster: {CELLS_PER_CLUSTER}")

    # Load data
    logging.info("\n=== Loading AUC Matrices ===")
    eos_data, mac_data, eos_clusters, mac_types = load_auc_matrices(scenic_dir)

    logging.info("\n=== Loading Gene Counts ===")
    gene_counts = load_gene_counts(scenic_dir, eos_clusters, mac_types)

    # Calculate mean AUC using subsampling
    logging.info("\n=== Calculating Mean AUC (Subsampling) ===")
    eos_cluster_df, mac_subtype_df, pooled_df, cell_counts, all_tfs = calculate_pooled_mean_auc_subsampling(
        eos_data, mac_data, eos_clusters, mac_types,
        cells_per_cluster=CELLS_PER_CLUSTER
    )

    # Save mean AUC matrices
    eos_cluster_df.to_csv(output_dir / 'eos_cluster_mean_auc.csv')
    mac_subtype_df.to_csv(output_dir / 'mac_subtype_mean_auc.csv')
    pooled_df.to_csv(output_dir / 'pooled_eos_mac_mean_auc.csv')

    # Create cell counts plot
    create_cell_counts_plot(cell_counts, CELLS_PER_CLUSTER, output_dir / 'cell_counts_summary.pdf')

    # Process each cutoff level
    logging.info("\n=== Finding Differential Regulons ===")

    for cutoff_level in args.cutoff_levels:
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing {cutoff_level.upper()}")
        logging.info(f"{'='*60}")

        results_df, eos_enriched, mac_enriched = find_differential_regulons(
            pooled_df, eos_cluster_df, mac_subtype_df,
            eos_clusters, mac_types, cutoff_level=cutoff_level
        )

        # Save results
        results_df.to_csv(output_dir / f'all_regulons_metrics_{cutoff_level}.csv', index=False)

        # Combined differential regulons
        differential = pd.concat([eos_enriched, mac_enriched], ignore_index=True)

        if len(differential) > 0:
            differential.to_csv(output_dir / f'differential_regulons_{cutoff_level}.csv', index=False)

            # Heatmap with Eos_all vs Mac_all fold change
            create_foldchange_heatmap_panel_a(
                differential, eos_clusters, mac_types,
                output_dir / f'panel_A_heatmap_foldchange_{cutoff_level}.pdf',
                gene_counts,
                title=f'Panel A: Eos_all vs Mac_all (n={len(differential)}, {cutoff_level})\nlog2 Fold Change',
                fc_cap=args.fc_cap
            )

            # Summary barplot
            create_summary_barplot(
                eos_enriched, mac_enriched,
                output_dir / f'panel_A_barplot_{cutoff_level}.pdf',
                cutoff_level
            )

        # Save separate files
        if len(eos_enriched) > 0:
            eos_enriched.to_csv(output_dir / f'eos_enriched_{cutoff_level}.csv', index=False)
        if len(mac_enriched) > 0:
            mac_enriched.to_csv(output_dir / f'mac_enriched_{cutoff_level}.csv', index=False)

        # Print summary
        print(f"\n{'='*60}")
        print(f"PANEL A SUMMARY ({cutoff_level.upper()})")
        print("="*60)
        print(f"Eos-enriched regulons: {len(eos_enriched)}")
        print(f"Mac-enriched regulons: {len(mac_enriched)}")
        print(f"Total differential: {len(differential)}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)

    print("\n" + "="*80)
    print("SCP COMMAND TO DOWNLOAD RESULTS:")
    print("="*80)
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
