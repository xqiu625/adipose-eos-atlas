#!/usr/bin/env python3
"""
Identify Eos and Mac subcluster-UNIQUE regulons using AUC values

Supports two methodologies:
1. Bootstrap (original method) - uses all cells to calculate mean AUC
2. Random sampling (Nature Immunology method) - randomly samples N cells per cluster
   Reference: https://www.nature.com/articles/s41590-020-0736-z#Sec19

Finds regulons that are:
- HIGH in ONE cluster (AUC >= high_auc_thresh)
- LOW in ALL other clusters (AUC <= low_auc_thresh)

This creates a clear diagonal pattern in the heatmap where each regulon
is visually distinct to its cluster.

Panel B: Eos subcluster regulons - compare subcluster with all other eos subclusters
Panel C: Mac subcluster regulons - compare Mac subcluster with all other Mac subclusters
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import logging
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Create custom colormap: low=#67a9cf (blue), high=#ef8a62 (red)
CUSTOM_CMAP = LinearSegmentedColormap.from_list('custom_blue_red', ['#67a9cf', '#ef8a62'])

# Cutoff thresholds for identifying subcluster-UNIQUE regulons
CUTOFF_LEVELS = {
    'relax': {
        'min_target_auc': 0.20,
        'min_difference': 0.05,
        'max_per_cluster': None
    },
    'moderate': {
        'min_target_auc': 0.25,
        'min_difference': 0.08,
        'max_per_cluster': None
    },
    'stringent': {
        'min_target_auc': 0.30,
        'min_difference': 0.10,
        'max_per_cluster': None
    },
    'very_stringent': {
        'min_target_auc': 0.35,
        'min_difference': 0.12,
        'max_per_cluster': None
    },
    'ultra_stringent': {
        'min_target_auc': 0.40,
        'min_difference': 0.15,
        'max_per_cluster': None
    },
    'top3_per_cluster': {
        'min_target_auc': 0.15,
        'min_difference': 0.0,
        'max_per_cluster': 3
    },
    'top2_per_cluster': {
        'min_target_auc': 0.15,
        'min_difference': 0.0,
        'max_per_cluster': 2
    }
}

# Random sampling parameters (Nature Immunology method)
RANDOM_SAMPLING_CONFIG = {
    'cells_per_cluster': 300,
    'n_iterations': 100,
    'random_seed': 42
}


def load_auc_matrices(scenic_dir, cell_type='eos'):
    """Load AUC matrices for specified cell type subclusters

    Directory structure:
    scenic_results/scenic_cleaned/eos/Eos_C0/step3_auc_matrix.csv
    scenic_results/scenic_cleaned/macrophages/M1/step3_auc_matrix.csv
    """

    if cell_type == 'eos':
        clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
        subdir = 'eos'
    elif cell_type == 'mac':
        # Macrophage subtypes: M1, M2_Functional, M2_Surface
        clusters = ['M1', 'M2_Functional', 'M2_Surface']
        subdir = 'macrophages'
    else:
        raise ValueError(f"Unknown cell_type: {cell_type}. Use 'eos' or 'mac'.")

    data = {}

    for cluster in clusters:
        auc_file = scenic_dir / subdir / cluster / "step3_auc_matrix.csv"
        if auc_file.exists():
            data[cluster] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {cluster}: {data[cluster].shape}")
        else:
            logging.warning(f"File not found: {auc_file}")

    total_cells = sum(df.shape[0] for df in data.values())
    logging.info(f"Total {cell_type.upper()} cells loaded: {total_cells}")

    if total_cells == 0:
        raise FileNotFoundError(f"No {cell_type} data loaded. Check scenic_dir path.")

    loaded_clusters = [c for c in clusters if c in data]
    return data, loaded_clusters


def load_gene_counts(scenic_dir, clusters, cell_type='eos'):
    """Load gene counts from step2_regulons.csv files"""
    gene_counts = {}
    subdir = 'eos' if cell_type == 'eos' else 'macrophages'

    for cluster in clusters:
        regulon_path = scenic_dir / subdir / cluster / "step2_regulons.csv"
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
            except Exception as e:
                logging.warning(f"Failed to load gene counts from {regulon_path}: {e}")

    logging.info(f"Total unique TFs with gene counts: {len(gene_counts)}")
    return gene_counts


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format 'TF_name(gene_count)'"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_mean_auc_bootstrap(data, clusters):
    """
    Calculate mean AUC using all cells (original bootstrap method).
    """
    logging.info("Calculating mean AUC using BOOTSTRAP method (all cells)")

    tf_to_regulon = {}
    all_tfs = set()
    for cluster in clusters:
        if cluster not in data:
            continue
        tf_to_regulon[cluster] = {}
        for col in data[cluster].columns:
            tf_name = extract_tf_name(col)
            tf_to_regulon[cluster][tf_name] = col
            all_tfs.add(tf_name)

    all_tfs = sorted(list(all_tfs))

    mean_auc_results = {cluster: {} for cluster in clusters}
    for cluster in clusters:
        if cluster not in data:
            continue
        for tf in all_tfs:
            if tf in tf_to_regulon.get(cluster, {}):
                regulon_col = tf_to_regulon[cluster][tf]
                mean_auc_results[cluster][tf] = data[cluster][regulon_col].mean()
            else:
                mean_auc_results[cluster][tf] = 0.0

    mean_auc_df = pd.DataFrame(mean_auc_results)
    return mean_auc_df


def calculate_mean_auc_random_sampling(data, clusters, cells_per_cluster=300,
                                        n_iterations=100, random_seed=42):
    """
    Calculate mean AUC using random sampling (Nature Immunology method).

    Randomly samples N cells from each cluster, calculates mean AUC,
    and repeats for stability. Returns the average across iterations.

    Reference: https://www.nature.com/articles/s41590-020-0736-z#Sec19
    """
    logging.info(f"Calculating mean AUC using RANDOM SAMPLING method")
    logging.info(f"  Cells per cluster: {cells_per_cluster}")
    logging.info(f"  Iterations: {n_iterations}")

    np.random.seed(random_seed)

    tf_to_regulon = {}
    all_tfs = set()
    for cluster in clusters:
        if cluster not in data:
            continue
        tf_to_regulon[cluster] = {}
        for col in data[cluster].columns:
            tf_name = extract_tf_name(col)
            tf_to_regulon[cluster][tf_name] = col
            all_tfs.add(tf_name)

    all_tfs = sorted(list(all_tfs))

    iteration_results = []

    for iteration in range(n_iterations):
        iter_results = {cluster: {} for cluster in clusters}

        for cluster in clusters:
            if cluster not in data:
                continue

            cluster_df = data[cluster]
            n_cells = cluster_df.shape[0]

            if n_cells >= cells_per_cluster:
                sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=False)
            else:
                sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=True)
                if iteration == 0:
                    logging.warning(f"  {cluster} has only {n_cells} cells, sampling with replacement")

            sampled_df = cluster_df.iloc[sampled_idx]

            for tf in all_tfs:
                if tf in tf_to_regulon.get(cluster, {}):
                    regulon_col = tf_to_regulon[cluster][tf]
                    iter_results[cluster][tf] = sampled_df[regulon_col].mean()
                else:
                    iter_results[cluster][tf] = 0.0

        iteration_results.append(pd.DataFrame(iter_results))

    stacked = np.stack([df.values for df in iteration_results], axis=0)
    mean_values = np.mean(stacked, axis=0)
    std_values = np.std(stacked, axis=0)

    mean_auc_df = pd.DataFrame(mean_values, index=all_tfs, columns=clusters)
    std_auc_df = pd.DataFrame(std_values, index=all_tfs, columns=clusters)

    logging.info(f"  Completed {n_iterations} iterations")

    return mean_auc_df, std_auc_df


def compare_methods(bootstrap_df, sampling_df, output_dir, cell_type='eos'):
    """
    Compare bootstrap and random sampling methods.
    """
    logging.info("\n=== Comparing Bootstrap vs Random Sampling Methods ===")

    common_tfs = bootstrap_df.index.intersection(sampling_df.index)
    common_clusters = bootstrap_df.columns.intersection(sampling_df.columns)

    bootstrap_aligned = bootstrap_df.loc[common_tfs, common_clusters]
    sampling_aligned = sampling_df.loc[common_tfs, common_clusters]

    boot_flat = bootstrap_aligned.values.flatten()
    samp_flat = sampling_aligned.values.flatten()

    correlation, p_value = stats.pearsonr(boot_flat, samp_flat)
    spearman_corr, spearman_p = stats.spearmanr(boot_flat, samp_flat)

    logging.info(f"  Pearson correlation: {correlation:.4f} (p={p_value:.2e})")
    logging.info(f"  Spearman correlation: {spearman_corr:.4f} (p={spearman_p:.2e})")

    diff = bootstrap_aligned - sampling_aligned
    mean_abs_diff = np.abs(diff.values).mean()
    max_abs_diff = np.abs(diff.values).max()

    logging.info(f"  Mean absolute difference: {mean_abs_diff:.4f}")
    logging.info(f"  Max absolute difference: {max_abs_diff:.4f}")

    # Create comparison figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax1 = axes[0]
    ax1.scatter(boot_flat, samp_flat, alpha=0.3, s=10)
    ax1.plot([0, 1], [0, 1], 'r--', linewidth=2, label='y=x')
    ax1.set_xlabel('Bootstrap Mean AUC', fontsize=12)
    ax1.set_ylabel('Random Sampling Mean AUC', fontsize=12)
    ax1.set_title(f'{cell_type.upper()}: Method Comparison\nr={correlation:.3f}', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    ax2 = axes[1]
    ax2.hist(diff.values.flatten(), bins=50, edgecolor='black', alpha=0.7)
    ax2.axvline(0, color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('Difference (Bootstrap - Random Sampling)', fontsize=12)
    ax2.set_ylabel('Frequency', fontsize=12)
    ax2.set_title(f'Distribution of Differences\nMean={mean_abs_diff:.4f}', fontsize=14, fontweight='bold')

    ax3 = axes[2]
    cluster_diffs = diff.abs().mean(axis=0)
    colors = plt.cm.RdYlGn_r(cluster_diffs / cluster_diffs.max())
    ax3.bar(range(len(cluster_diffs)), cluster_diffs, color=colors, edgecolor='black')
    ax3.set_xticks(range(len(cluster_diffs)))
    ax3.set_xticklabels(cluster_diffs.index, rotation=45, ha='right')
    ax3.set_xlabel('Cluster', fontsize=12)
    ax3.set_ylabel('Mean Absolute Difference', fontsize=12)
    ax3.set_title('Difference by Cluster', fontsize=14, fontweight='bold')

    plt.tight_layout()
    output_path = output_dir / f'{cell_type}_method_comparison.pdf'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()

    stats_df = pd.DataFrame({
        'Metric': ['Pearson_r', 'Pearson_p', 'Spearman_r', 'Spearman_p',
                   'Mean_Abs_Diff', 'Max_Abs_Diff'],
        'Value': [correlation, p_value, spearman_corr, spearman_p,
                  mean_abs_diff, max_abs_diff]
    })
    stats_df.to_csv(output_dir / f'{cell_type}_method_comparison_stats.csv', index=False)

    return {
        'pearson_r': correlation,
        'spearman_r': spearman_corr,
        'mean_abs_diff': mean_abs_diff,
        'max_abs_diff': max_abs_diff
    }


def identify_subcluster_unique_regulons(mean_auc_df, clusters, cutoff_level='relax'):
    """
    Identify regulons that are UNIQUELY enriched in each subcluster.
    """
    criteria = CUTOFF_LEVELS[cutoff_level]
    min_target_auc = criteria['min_target_auc']
    min_difference = criteria['min_difference']
    max_per_cluster = criteria.get('max_per_cluster', None)

    unique_regulons = {}
    summary_rows = []
    all_tfs = list(mean_auc_df.index)

    for target_cluster in clusters:
        if target_cluster not in mean_auc_df.columns:
            continue

        other_clusters = [c for c in clusters if c != target_cluster and c in mean_auc_df.columns]
        cluster_results = []

        for tf in all_tfs:
            target_auc = mean_auc_df.loc[tf, target_cluster]
            other_aucs = [mean_auc_df.loc[tf, c] for c in other_clusters]
            max_other_auc = max(other_aucs) if other_aucs else 0
            mean_other_auc = np.mean(other_aucs) if other_aucs else 0
            uniqueness_score = target_auc - max_other_auc

            cluster_results.append({
                'tf': tf,
                'target_auc': target_auc,
                'max_other_auc': max_other_auc,
                'mean_other_auc': mean_other_auc,
                'uniqueness_score': uniqueness_score
            })

        cluster_results.sort(key=lambda x: x['uniqueness_score'], reverse=True)

        if max_per_cluster is not None:
            filtered = [r for r in cluster_results
                       if r['target_auc'] >= min_target_auc
                       and r['uniqueness_score'] > 0][:max_per_cluster]
        else:
            filtered = [r for r in cluster_results
                       if r['target_auc'] >= min_target_auc
                       and r['uniqueness_score'] >= min_difference]

        unique_regulons[target_cluster] = filtered

        summary_rows.append({
            'cluster': target_cluster,
            'n_unique_regulons': len(filtered),
            'regulons': ', '.join([r['tf'] for r in filtered[:10]]) + ('...' if len(filtered) > 10 else '')
        })

        logging.info(f"{target_cluster}: {len(filtered)} unique regulons")

    summary_df = pd.DataFrame(summary_rows)
    return unique_regulons, summary_df


def create_heatmap(unique_regulons, mean_auc_df, output_path, gene_counts, clusters,
                   title="Subcluster-Unique Regulons"):
    """Create heatmap showing mean AUC values for subcluster-unique regulons."""

    all_unique = []
    cluster_boundaries = []
    current_pos = 0

    for cluster in clusters:
        if cluster not in unique_regulons:
            continue

        cluster_regs = unique_regulons[cluster]
        if cluster_regs:
            start = current_pos
            for reg_info in cluster_regs:
                tf = reg_info['tf']
                if tf in mean_auc_df.index and tf not in all_unique:
                    all_unique.append(tf)
                    current_pos += 1
            if current_pos > start:
                cluster_boundaries.append((cluster, start, current_pos))

    if len(all_unique) == 0:
        logging.warning("No unique regulons found to plot")
        return

    logging.info(f"Total regulons in heatmap: {len(all_unique)}")

    plot_data = mean_auc_df.loc[all_unique]

    ylabels = []
    for tf in all_unique:
        if tf in gene_counts:
            ylabels.append(f"{tf}({gene_counts[tf]})")
        else:
            ylabels.append(tf)

    fig_height = max(6, len(all_unique) * 0.4)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    im = ax.imshow(plot_data.values, cmap=CUSTOM_CMAP, aspect='auto', vmin=0, vmax=1.0)

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('Mean AUC', fontsize=12, fontweight='bold')

    ax.set_xticks(np.arange(len(plot_data.columns)))
    ax.set_xticklabels(plot_data.columns, rotation=45, ha='right', fontsize=10)

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=9)

    ax.set_xlabel('Subcluster', fontsize=12, fontweight='bold')
    ax.set_ylabel('Regulon (gene count)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_barplot(summary_df, output_path, cell_type='eos'):
    """Create bar plot showing number of unique regulons per subcluster."""

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.Set2(np.linspace(0, 1, len(summary_df)))

    bars = ax.bar(summary_df['cluster'], summary_df['n_unique_regulons'],
                  color=colors, edgecolor='black')

    ax.set_xlabel(f'{cell_type.upper()} Subcluster', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Unique Regulons', fontsize=12, fontweight='bold')
    ax.set_title(f'Unique Regulons per {cell_type.upper()} Subcluster', fontsize=14, fontweight='bold')

    for bar, val in zip(bars, summary_df['n_unique_regulons']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha='center', fontsize=11, fontweight='bold')

    max_val = max(summary_df['n_unique_regulons'])
    ax.set_ylim(0, max(max_val * 1.15, 5))
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def run_analysis(scenic_dir, output_dir, cell_type='eos', method='both', cutoff_levels=None):
    """
    Run the complete analysis for a given cell type.
    """
    if cutoff_levels is None:
        cutoff_levels = ['relax', 'stringent', 'top3_per_cluster']

    logging.info("="*80)
    logging.info(f"ANALYSIS: {cell_type.upper()} SUBCLUSTER-UNIQUE REGULONS")
    logging.info(f"Method: {method}")
    logging.info("="*80)

    cell_output_dir = output_dir / cell_type
    cell_output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logging.info(f"\n=== Load {cell_type.upper()} AUC Matrices ===")
    data, clusters = load_auc_matrices(scenic_dir, cell_type=cell_type)

    logging.info(f"\n=== Load Gene Counts ===")
    gene_counts = load_gene_counts(scenic_dir, clusters, cell_type=cell_type)

    # Calculate mean AUC
    logging.info(f"\n=== Calculate Mean AUC ===")
    results = {}

    if method in ['bootstrap', 'both']:
        bootstrap_df = calculate_mean_auc_bootstrap(data, clusters)
        bootstrap_df.to_csv(cell_output_dir / f'{cell_type}_mean_auc_bootstrap.csv')
        results['bootstrap'] = bootstrap_df

    if method in ['random_sampling', 'both']:
        sampling_df, std_df = calculate_mean_auc_random_sampling(
            data, clusters,
            cells_per_cluster=RANDOM_SAMPLING_CONFIG['cells_per_cluster'],
            n_iterations=RANDOM_SAMPLING_CONFIG['n_iterations'],
            random_seed=RANDOM_SAMPLING_CONFIG['random_seed']
        )
        sampling_df.to_csv(cell_output_dir / f'{cell_type}_mean_auc_random_sampling.csv')
        std_df.to_csv(cell_output_dir / f'{cell_type}_std_auc_random_sampling.csv')
        results['random_sampling'] = sampling_df

    # Compare methods
    if method == 'both':
        logging.info(f"\n=== Compare Methods ===")
        comparison_stats = compare_methods(
            results['bootstrap'],
            results['random_sampling'],
            cell_output_dir,
            cell_type=cell_type
        )

        if comparison_stats['pearson_r'] > 0.95 and comparison_stats['mean_abs_diff'] < 0.02:
            logging.info("\nMethods are highly concordant (r > 0.95, diff < 0.02)")
            logging.info("Using BOOTSTRAP method for downstream analysis")
            primary_method = 'bootstrap'
        else:
            logging.info("\nMethods show some differences - running both")
            primary_method = 'both'
    else:
        primary_method = method

    # Find unique regulons
    logging.info(f"\n=== Find Unique Regulons ===")
    methods_to_run = [primary_method] if primary_method != 'both' else ['bootstrap', 'random_sampling']

    for method_name in methods_to_run:
        mean_auc_df = results[method_name]
        method_suffix = f"_{method_name}" if len(methods_to_run) > 1 else ""

        for cutoff_level in cutoff_levels:
            logging.info(f"\nProcessing {cutoff_level.upper()} ({method_name})")

            unique_regulons, summary_df = identify_subcluster_unique_regulons(
                mean_auc_df, clusters, cutoff_level=cutoff_level
            )

            # Save results
            summary_df.to_csv(
                cell_output_dir / f'{cell_type}_unique_summary_{cutoff_level}{method_suffix}.csv',
                index=False
            )

            for cluster, regs_info in unique_regulons.items():
                if regs_info:
                    cluster_df = pd.DataFrame(regs_info)
                    cluster_df.to_csv(
                        cell_output_dir / f'{cluster}_unique_regulons_{cutoff_level}{method_suffix}.csv',
                        index=False
                    )

            # Create visualizations
            if summary_df['n_unique_regulons'].sum() > 0:
                create_barplot(
                    summary_df,
                    cell_output_dir / f'{cell_type}_unique_barplot_{cutoff_level}{method_suffix}.pdf',
                    cell_type=cell_type
                )

                create_heatmap(
                    unique_regulons, mean_auc_df,
                    cell_output_dir / f'{cell_type}_unique_heatmap_{cutoff_level}{method_suffix}.pdf',
                    gene_counts, clusters,
                    title=f'{cell_type.upper()} Subcluster-Unique Regulons ({method_name})'
                )

            print(f"\n{'='*60}")
            print(f"{cell_type.upper()} UNIQUE REGULONS ({cutoff_level}, {method_name})")
            print("="*60)
            print(summary_df.to_string(index=False))

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Subcluster-unique regulon analysis')
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--output-dir', type=str, default='results',
                        help='Path to output directory')
    parser.add_argument('--cell-type', type=str, choices=['eos', 'mac', 'both'], default='both',
                        help='Cell type to analyze')
    parser.add_argument('--method', type=str, choices=['bootstrap', 'random_sampling', 'both'], default='both',
                        help='Method to use')
    parser.add_argument('--cells-per-cluster', type=int, default=300,
                        help='Cells to sample per cluster for random sampling')
    parser.add_argument('--n-iterations', type=int, default=100,
                        help='Number of iterations for random sampling')
    parser.add_argument('--cutoff-levels', type=str, nargs='+',
                        default=['relax', 'stringent', 'top3_per_cluster'],
                        help='Cutoff levels to use')
    args = parser.parse_args()

    RANDOM_SAMPLING_CONFIG['cells_per_cluster'] = args.cells_per_cluster
    RANDOM_SAMPLING_CONFIG['n_iterations'] = args.n_iterations

    scenic_dir = Path(args.scenic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("SUBCLUSTER-UNIQUE REGULON ANALYSIS")
    logging.info("Bootstrap vs Random Sampling (Nature Immunology) Comparison")
    logging.info("="*80)

    cell_types = ['eos', 'mac'] if args.cell_type == 'both' else [args.cell_type]

    for cell_type in cell_types:
        try:
            run_analysis(
                scenic_dir, output_dir,
                cell_type=cell_type,
                method=args.method,
                cutoff_levels=args.cutoff_levels
            )
        except FileNotFoundError as e:
            logging.error(f"Could not analyze {cell_type}: {e}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)


if __name__ == "__main__":
    main()
