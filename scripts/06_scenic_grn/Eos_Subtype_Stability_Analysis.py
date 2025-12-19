#!/usr/bin/env python3
"""
Eos Subtype Stability Analysis Across Conditions
=================================================
Quantify how similar/stable Eos subclusters are across different conditions:
- WT Female vs WT Male
- WT vs KO (within each sex)
- All pairwise comparisons

Metrics calculated:
1. Gene expression overlap (% shared expressed genes)
2. Correlation of mean expression profiles
3. Jaccard similarity of marker genes
4. Top expressed genes overlap

"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import cosine, correlation
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Paths
H5AD_FILE = Path("/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/scenic_project/inputs/eos_all.h5ad")
OUTPUT_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/subtype_stability")

# Eos subclusters (note: in h5ad they are stored as 'C0', 'C1', etc.)
EOS_CLUSTERS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

# Conditions
CONDITIONS = ['WT_Female', 'WT_Male', 'KO_Female', 'KO_Male']


def load_eos_data(h5ad_file):
    """Load Eos h5ad data"""
    logging.info(f"Loading data from: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    logging.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Print available metadata columns
    logging.info(f"Available obs columns: {list(adata.obs.columns)}")

    return adata


def get_condition_mask(adata, genotype, sex):
    """Get mask for cells matching genotype and sex"""
    genotype_mask = adata.obs['genotype'] == genotype
    sex_mask = adata.obs['Sex'] == sex
    return genotype_mask & sex_mask


def get_cluster_condition_cells(adata, cluster, genotype, sex):
    """Get cells for specific cluster and condition"""
    if 'eos_cluster' in adata.obs.columns:
        cluster_mask = adata.obs['eos_cluster'] == cluster
    elif 'cluster' in adata.obs.columns:
        cluster_mask = adata.obs['cluster'] == cluster
    else:
        # Try to infer from leiden or other clustering
        for col in ['leiden', 'louvain', 'cell_type']:
            if col in adata.obs.columns:
                cluster_mask = adata.obs[col] == cluster
                break
        else:
            logging.error("No cluster column found")
            return None

    condition_mask = get_condition_mask(adata, genotype, sex)
    combined_mask = cluster_mask & condition_mask

    return adata[combined_mask].copy()


def calculate_expressed_genes(adata, min_cells_pct=0.1, min_mean_expr=0.1):
    """
    Calculate set of expressed genes in a cell population.

    Args:
        adata: AnnData object
        min_cells_pct: Minimum percentage of cells expressing the gene
        min_mean_expr: Minimum mean expression level

    Returns:
        Set of expressed gene names
    """
    if adata.n_obs == 0:
        return set()

    # Convert to dense if sparse
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    # Calculate metrics
    pct_cells = (X > 0).sum(axis=0) / X.shape[0]
    mean_expr = X.mean(axis=0)

    # Handle different array types
    if hasattr(pct_cells, 'A1'):
        pct_cells = pct_cells.A1
    if hasattr(mean_expr, 'A1'):
        mean_expr = mean_expr.A1

    pct_cells = np.array(pct_cells).flatten()
    mean_expr = np.array(mean_expr).flatten()

    # Filter genes
    expressed_mask = (pct_cells >= min_cells_pct) & (mean_expr >= min_mean_expr)
    expressed_genes = set(adata.var_names[expressed_mask])

    return expressed_genes


def calculate_mean_expression(adata):
    """Calculate mean expression profile"""
    if adata.n_obs == 0:
        return pd.Series(dtype=float)

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    mean_expr = X.mean(axis=0)
    if hasattr(mean_expr, 'A1'):
        mean_expr = mean_expr.A1

    return pd.Series(np.array(mean_expr).flatten(), index=adata.var_names)


def get_top_genes(adata, n_top=100):
    """Get top N expressed genes by mean expression"""
    mean_expr = calculate_mean_expression(adata)
    return set(mean_expr.nlargest(n_top).index)


def get_marker_genes(adata, n_markers=50):
    """
    Get marker genes (highly variable + highly expressed).
    Uses a simple approach: top genes by expression * variance.
    """
    if adata.n_obs == 0:
        return set()

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    mean_expr = X.mean(axis=0)
    var_expr = X.var(axis=0)

    if hasattr(mean_expr, 'A1'):
        mean_expr = mean_expr.A1
    if hasattr(var_expr, 'A1'):
        var_expr = var_expr.A1

    mean_expr = np.array(mean_expr).flatten()
    var_expr = np.array(var_expr).flatten()

    # Score = mean * variance (to get highly expressed AND variable genes)
    score = mean_expr * var_expr
    top_idx = np.argsort(score)[-n_markers:]

    return set(adata.var_names[top_idx])


def jaccard_similarity(set1, set2):
    """Calculate Jaccard similarity between two sets"""
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    if len(set1) == 0 or len(set2) == 0:
        return 0.0

    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union


def overlap_percentage(set1, set2):
    """Calculate what percentage of set1 is also in set2"""
    if len(set1) == 0:
        return 0.0
    return len(set1 & set2) / len(set1) * 100


def analyze_cluster_stability(adata, cluster, output_dir):
    """
    Analyze stability of a single cluster across conditions.
    """
    logging.info(f"\n{'='*60}")
    logging.info(f"Analyzing {cluster}")
    logging.info(f"{'='*60}")

    results = {
        'cluster': cluster,
        'conditions': {},
        'pairwise': []
    }

    # Get data for each condition
    condition_data = {}
    for genotype in ['WT', 'KO']:
        for sex in ['Female', 'Male']:
            cond_name = f"{genotype}_{sex}"
            subset = get_cluster_condition_cells(adata, cluster, genotype, sex)

            if subset is None or subset.n_obs == 0:
                logging.warning(f"  No cells for {cond_name}")
                continue

            condition_data[cond_name] = {
                'adata': subset,
                'n_cells': subset.n_obs,
                'expressed_genes': calculate_expressed_genes(subset),
                'mean_expr': calculate_mean_expression(subset),
                'top_genes': get_top_genes(subset, n_top=100),
                'marker_genes': get_marker_genes(subset, n_markers=50)
            }

            results['conditions'][cond_name] = {
                'n_cells': subset.n_obs,
                'n_expressed_genes': len(condition_data[cond_name]['expressed_genes']),
                'n_top_genes': len(condition_data[cond_name]['top_genes']),
                'n_marker_genes': len(condition_data[cond_name]['marker_genes'])
            }

            logging.info(f"  {cond_name}: {subset.n_obs} cells, "
                        f"{len(condition_data[cond_name]['expressed_genes'])} expressed genes")

    # Pairwise comparisons
    cond_names = list(condition_data.keys())

    for i, cond1 in enumerate(cond_names):
        for cond2 in cond_names[i+1:]:
            data1 = condition_data[cond1]
            data2 = condition_data[cond2]

            # Gene overlap metrics
            expressed_jaccard = jaccard_similarity(data1['expressed_genes'], data2['expressed_genes'])
            expressed_overlap_1to2 = overlap_percentage(data1['expressed_genes'], data2['expressed_genes'])
            expressed_overlap_2to1 = overlap_percentage(data2['expressed_genes'], data1['expressed_genes'])

            top_jaccard = jaccard_similarity(data1['top_genes'], data2['top_genes'])
            top_overlap_1to2 = overlap_percentage(data1['top_genes'], data2['top_genes'])

            marker_jaccard = jaccard_similarity(data1['marker_genes'], data2['marker_genes'])
            marker_overlap_1to2 = overlap_percentage(data1['marker_genes'], data2['marker_genes'])

            # Expression correlation
            common_genes = list(set(data1['mean_expr'].index) & set(data2['mean_expr'].index))
            if len(common_genes) > 0:
                expr1 = data1['mean_expr'][common_genes].values
                expr2 = data2['mean_expr'][common_genes].values

                # Pearson correlation
                pearson_r, pearson_p = stats.pearsonr(expr1, expr2)

                # Spearman correlation
                spearman_r, spearman_p = stats.spearmanr(expr1, expr2)

                # Cosine similarity
                cosine_sim = 1 - cosine(expr1, expr2) if np.any(expr1) and np.any(expr2) else 0
            else:
                pearson_r = spearman_r = cosine_sim = 0
                pearson_p = spearman_p = 1

            comparison = {
                'condition1': cond1,
                'condition2': cond2,
                'n_cells_1': data1['n_cells'],
                'n_cells_2': data2['n_cells'],
                'expressed_genes_jaccard': expressed_jaccard,
                'expressed_genes_overlap_1to2': expressed_overlap_1to2,
                'expressed_genes_overlap_2to1': expressed_overlap_2to1,
                'top100_genes_jaccard': top_jaccard,
                'top100_genes_overlap': top_overlap_1to2,
                'marker_genes_jaccard': marker_jaccard,
                'marker_genes_overlap': marker_overlap_1to2,
                'pearson_correlation': pearson_r,
                'pearson_pvalue': pearson_p,
                'spearman_correlation': spearman_r,
                'spearman_pvalue': spearman_p,
                'cosine_similarity': cosine_sim
            }

            results['pairwise'].append(comparison)

            logging.info(f"\n  {cond1} vs {cond2}:")
            logging.info(f"    Expressed genes Jaccard: {expressed_jaccard:.3f}")
            logging.info(f"    Expressed genes overlap: {expressed_overlap_1to2:.1f}% / {expressed_overlap_2to1:.1f}%")
            logging.info(f"    Top 100 genes Jaccard: {top_jaccard:.3f}")
            logging.info(f"    Marker genes Jaccard: {marker_jaccard:.3f}")
            logging.info(f"    Pearson correlation: {pearson_r:.3f} (p={pearson_p:.2e})")
            logging.info(f"    Spearman correlation: {spearman_r:.3f}")
            logging.info(f"    Cosine similarity: {cosine_sim:.3f}")

    return results


def create_summary_figure(all_results, output_dir):
    """Create summary visualization"""

    # Collect pairwise data
    rows = []
    for cluster_results in all_results:
        cluster = cluster_results['cluster']
        for comp in cluster_results['pairwise']:
            rows.append({
                'cluster': cluster,
                'comparison': f"{comp['condition1']} vs {comp['condition2']}",
                'Expressed Genes Jaccard': comp['expressed_genes_jaccard'],
                'Top 100 Genes Jaccard': comp['top100_genes_jaccard'],
                'Marker Genes Jaccard': comp['marker_genes_jaccard'],
                'Pearson Correlation': comp['pearson_correlation'],
                'Spearman Correlation': comp['spearman_correlation'],
                'Cosine Similarity': comp['cosine_similarity']
            })

    df = pd.DataFrame(rows)

    if len(df) == 0:
        logging.warning("No data for summary figure")
        return

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    metrics = ['Expressed Genes Jaccard', 'Top 100 Genes Jaccard', 'Marker Genes Jaccard',
               'Pearson Correlation', 'Spearman Correlation', 'Cosine Similarity']

    for ax, metric in zip(axes.flat, metrics):
        pivot = df.pivot(index='cluster', columns='comparison', values=metric)
        sns.heatmap(pivot, ax=ax, annot=True, fmt='.2f', cmap='RdYlGn',
                    vmin=0, vmax=1, center=0.5, cbar_kws={'shrink': 0.5})
        ax.set_title(metric, fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        plt.setp(ax.get_yticklabels(), fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'subtype_stability_heatmaps.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'subtype_stability_heatmaps.png', dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_dir / 'subtype_stability_heatmaps.pdf'}")
    plt.close()

    # Create bar plot for key comparisons
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # WT Female vs WT Male
    wt_sex_df = df[df['comparison'] == 'WT_Female vs WT_Male'].copy()
    if len(wt_sex_df) > 0:
        ax = axes[0]
        x = np.arange(len(wt_sex_df))
        width = 0.15

        for i, metric in enumerate(['Expressed Genes Jaccard', 'Top 100 Genes Jaccard',
                                     'Pearson Correlation', 'Cosine Similarity']):
            ax.bar(x + i*width, wt_sex_df[metric], width, label=metric.replace(' ', '\n'))

        ax.set_xticks(x + 1.5*width)
        ax.set_xticklabels(wt_sex_df['cluster'], rotation=45, ha='right')
        ax.set_ylabel('Similarity Score')
        ax.set_title('WT Female vs WT Male\n(Sex Difference)', fontsize=12, fontweight='bold')
        ax.legend(loc='lower right', fontsize=7)
        ax.set_ylim(0, 1.1)
        ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5, label='80% threshold')

    # WT vs KO (averaged across sexes)
    wt_ko_df = df[df['comparison'].str.contains('WT.*KO|KO.*WT')].copy()
    if len(wt_ko_df) > 0:
        # Average by cluster
        wt_ko_avg = wt_ko_df.groupby('cluster').mean(numeric_only=True).reset_index()

        ax = axes[1]
        x = np.arange(len(wt_ko_avg))
        width = 0.15

        for i, metric in enumerate(['Expressed Genes Jaccard', 'Top 100 Genes Jaccard',
                                     'Pearson Correlation', 'Cosine Similarity']):
            ax.bar(x + i*width, wt_ko_avg[metric], width, label=metric.replace(' ', '\n'))

        ax.set_xticks(x + 1.5*width)
        ax.set_xticklabels(wt_ko_avg['cluster'], rotation=45, ha='right')
        ax.set_ylabel('Similarity Score')
        ax.set_title('WT vs KO\n(Genotype Difference)', fontsize=12, fontweight='bold')
        ax.legend(loc='lower right', fontsize=7)
        ax.set_ylim(0, 1.1)
        ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_dir / 'subtype_stability_barplots.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'subtype_stability_barplots.png', dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_dir / 'subtype_stability_barplots.pdf'}")
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Eos Subtype Stability Analysis')
    parser.add_argument('--h5ad', type=str, default=str(H5AD_FILE))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('--clusters', nargs='+', default=EOS_CLUSTERS)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("EOS SUBTYPE STABILITY ANALYSIS")
    logging.info("="*80)

    # Load data
    adata = load_eos_data(args.h5ad)

    # Analyze each cluster
    all_results = []
    for cluster in args.clusters:
        results = analyze_cluster_stability(adata, cluster, output_dir)
        all_results.append(results)

    # Save detailed results
    all_pairwise = []
    for cluster_results in all_results:
        for comp in cluster_results['pairwise']:
            comp['cluster'] = cluster_results['cluster']
            all_pairwise.append(comp)

    pairwise_df = pd.DataFrame(all_pairwise)
    pairwise_df.to_csv(output_dir / 'pairwise_comparisons.csv', index=False)
    logging.info(f"\nSaved: {output_dir / 'pairwise_comparisons.csv'}")

    # Create summary figures
    create_summary_figure(all_results, output_dir)

    # Print summary statistics
    logging.info("\n" + "="*80)
    logging.info("SUMMARY STATISTICS")
    logging.info("="*80)

    # Key comparisons
    key_comparisons = [
        ('WT_Female', 'WT_Male', 'Sex difference (WT)'),
        ('KO_Female', 'KO_Male', 'Sex difference (KO)'),
        ('WT_Female', 'KO_Female', 'Genotype difference (Female)'),
        ('WT_Male', 'KO_Male', 'Genotype difference (Male)')
    ]

    for cond1, cond2, label in key_comparisons:
        if len(pairwise_df) == 0 or 'condition1' not in pairwise_df.columns:
            continue
        subset = pairwise_df[((pairwise_df['condition1'] == cond1) & (pairwise_df['condition2'] == cond2)) |
                             ((pairwise_df['condition1'] == cond2) & (pairwise_df['condition2'] == cond1))]

        if len(subset) > 0:
            logging.info(f"\n{label}:")
            logging.info(f"  Mean Expressed Genes Jaccard: {subset['expressed_genes_jaccard'].mean():.3f} ± {subset['expressed_genes_jaccard'].std():.3f}")
            logging.info(f"  Mean Top 100 Genes Jaccard: {subset['top100_genes_jaccard'].mean():.3f} ± {subset['top100_genes_jaccard'].std():.3f}")
            logging.info(f"  Mean Pearson Correlation: {subset['pearson_correlation'].mean():.3f} ± {subset['pearson_correlation'].std():.3f}")
            logging.info(f"  Mean Cosine Similarity: {subset['cosine_similarity'].mean():.3f} ± {subset['cosine_similarity'].std():.3f}")

    # Overall stability
    logging.info(f"\n{'='*80}")
    logging.info("OVERALL CLUSTER STABILITY (across all conditions):")
    logging.info(f"{'='*80}")

    for cluster in args.clusters:
        if len(pairwise_df) == 0 or 'cluster' not in pairwise_df.columns:
            continue
        cluster_data = pairwise_df[pairwise_df['cluster'] == cluster]
        if len(cluster_data) > 0:
            logging.info(f"\n{cluster}:")
            logging.info(f"  Avg Expressed Genes Jaccard: {cluster_data['expressed_genes_jaccard'].mean():.3f}")
            logging.info(f"  Avg Top 100 Genes Jaccard: {cluster_data['top100_genes_jaccard'].mean():.3f}")
            logging.info(f"  Avg Pearson Correlation: {cluster_data['pearson_correlation'].mean():.3f}")
            logging.info(f"  Avg Cosine Similarity: {cluster_data['cosine_similarity'].mean():.3f}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)
    print(f"\nOutput directory: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir}" .')


if __name__ == "__main__":
    main()
