#!/usr/bin/env python3
"""
SCENIC Condition-Specific Analysis
===================================
Author: Xinru Qiu
Date: December 2024

Problem: SCENIC AUC matrices use Cell_0, Cell_1, etc. but we need to link
         back to original cell barcodes to compare across conditions.

Solution: Extract cell order from original h5ad input files (obs_names),
          then merge AUC data with condition metadata from Seurat.

Analyses:
1. Regulon activity comparison by condition (WT vs KO, Male vs Female, etc.)
2. Condition-specific GRN networks
3. Differential regulon activity testing
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import logging
import scanpy as sc

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# =====================================
# CONFIGURATION
# =====================================

SCENIC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic")
SEURAT_METADATA = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/cell_metadata_all.csv")
OUTPUT_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/condition_analysis")

# H5AD input files directories (contains original cell order)
# Different cell types have h5ad files in different locations
H5AD_INPUT_DIR_EOS = Path("/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/scenic_project/inputs")
H5AD_INPUT_DIR_MAC = Path("/bigdata/godziklab/shared/Xinru/302025/302025_git/scripts/R11_20250818/inputs/macrophage_subtypes")

# SCENIC output structure:
#   /scenic/eos_sub/Eos_C0/, Eos_C1/, ..., Eos_C6/  (7 clusters)
#   /scenic/macrophage_subtypes/M1_Macrophages/, M2_Macrophages_Functional/, M2_Macrophages_Surface_Only/
#   /scenic/eos_all/  (all eosinophils combined)
#   /scenic/non-eos/mac_all/  (all macrophages combined)

CELL_TYPE_CONFIG = {
    'eos': {
        # Eosinophil subclusters - SCENIC results in /scenic/eos_sub/
        # SCENIC folders: Eos_C0, Eos_C1, etc.
        # Seurat eos_cluster column: C0, C1, etc.
        'clusters': ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6'],
        'scenic_subdir': 'eos_sub',  # AUC files in /scenic/eos_sub/Eos_CX/
        'cluster_col': 'eos_cluster',
        'h5ad_dir': H5AD_INPUT_DIR_EOS,
        # Map SCENIC cluster name (Eos_CX) to Seurat cluster name (CX)
        'seurat_cluster_mapping': {
            'Eos_C0': 'C0', 'Eos_C1': 'C1', 'Eos_C2': 'C2', 'Eos_C3': 'C3',
            'Eos_C4': 'C4', 'Eos_C5': 'C5', 'Eos_C6': 'C6',
        },
        # Map cluster name to h5ad file (some may be missing)
        'h5ad_files': {
            'Eos_C0': 'Eos_C0.h5ad',
            'Eos_C1': 'Eos_C1.h5ad',
            'Eos_C2': 'Eos_C2.h5ad',
            'Eos_C3': 'Eos_C3.h5ad',
            'Eos_C4': 'Eos_C4.h5ad',
            'Eos_C5': 'Eos_C5.h5ad',
            'Eos_C6': 'Eos_C6.h5ad',
        }
    },
    'macrophages': {
        # Macrophage subtypes - SCENIC results in /scenic/macrophage_subtypes/
        'clusters': ['M1_Macrophages', 'M2_Macrophages_Functional', 'M2_Macrophages_Surface_Only'],
        'scenic_subdir': 'macrophage_subtypes',  # AUC files in /scenic/macrophage_subtypes/
        'cluster_col': 'macrophage_type',
        'h5ad_dir': H5AD_INPUT_DIR_MAC,
        # h5ad files are in R11 inputs folder with matching names
        'h5ad_files': {
            'M1_Macrophages': 'M1_Macrophages.h5ad',
            'M2_Macrophages_Functional': 'M2_Macrophages_Functional.h5ad',
            'M2_Macrophages_Surface_Only': 'M2_Macrophages_Surface_Only.h5ad',
        },
        'combine_for_condition_analysis': True
    }
}

# Condition order for plotting
CONDITION_ORDER = [
    'WT_Female_HFD',
    'WT_Male_HFD',
    'KO_Female_HFD',
    'KO_Male_HFD',
    'WT_Female_HFDInfected',
    'KO_Female_HFDInfected'
]

CONDITION_COLORS = {
    'WT_Female_HFD': '#1f77b4',
    'WT_Male_HFD': '#2ca02c',
    'KO_Female_HFD': '#d62728',
    'KO_Male_HFD': '#ff7f0e',
    'WT_Female_HFDInfected': '#9467bd',
    'KO_Female_HFDInfected': '#8c564b'
}


# =====================================
# FUNCTIONS
# =====================================

def load_cell_mapping_from_file(cell_file):
    """
    Load cell barcode order from a text file or h5ad file.

    SCENIC processes cells in the same order as the input,
    so line[0] = Cell_0, line[1] = Cell_1, etc.

    Returns: DataFrame with scenic_id -> cell_barcode mapping
    """
    cell_file = Path(cell_file)

    if not cell_file.exists():
        logging.warning(f"Cell file not found: {cell_file}")
        return None

    try:
        if cell_file.suffix == '.txt':
            # Read from text file (one barcode per line)
            with open(cell_file) as f:
                cell_barcodes = [line.strip() for line in f if line.strip()]
            logging.info(f"Loaded {len(cell_barcodes)} cell barcodes from {cell_file.name}")

        elif cell_file.suffix == '.h5ad':
            # Read from h5ad file
            adata = sc.read_h5ad(cell_file, backed='r')
            cell_barcodes = list(adata.obs_names)
            adata.file.close()
            logging.info(f"Extracted {len(cell_barcodes)} cell barcodes from {cell_file.name}")

        else:
            logging.warning(f"Unsupported file type: {cell_file.suffix}")
            return None

        mapping = pd.DataFrame({
            'scenic_id': [f'Cell_{i}' for i in range(len(cell_barcodes))],
            'cell_barcode': cell_barcodes
        })

        return mapping

    except Exception as e:
        logging.error(f"Failed to read cell file {cell_file}: {e}")
        return None


def get_cell_file_for_cluster(cell_type, cluster, config):
    """
    Get the cell list file path for a given cluster.
    Prefers .txt files over .h5ad files.
    """
    h5ad_dir = config.get('h5ad_dir')
    if h5ad_dir is None:
        h5ad_dir = H5AD_INPUT_DIR_EOS

    # First try _cells.txt file (preferred - faster to read)
    txt_file = h5ad_dir / f"{cluster}_cells.txt"
    if txt_file.exists():
        return txt_file

    # Fall back to h5ad file
    h5ad_mapping = config.get('h5ad_files', {})
    h5ad_name = h5ad_mapping.get(cluster, f"{cluster}.h5ad")
    h5ad_file = h5ad_dir / h5ad_name
    if h5ad_file.exists():
        return h5ad_file

    return None


def load_auc_with_barcodes(scenic_dir, cell_type, cluster, config):
    """
    Load AUC matrix. If cell IDs are Cell_0, Cell_1, etc., map to original barcodes.
    If cell IDs are already full barcodes, use them directly.
    """
    # Get the scenic subdirectory from config
    scenic_subdir = config.get('scenic_subdir', cell_type)

    # Find AUC file - try different path patterns
    auc_patterns = [
        scenic_dir / scenic_subdir / cluster / "step3_auc_matrix.csv",
        scenic_dir / cell_type / cluster / "step3_auc_matrix.csv",
        scenic_dir / "eos_sub" / cluster / "step3_auc_matrix.csv",
        scenic_dir / "macrophage_subtypes" / cluster / "step3_auc_matrix.csv",
        scenic_dir / "non-eos" / cluster / "step3_auc_matrix.csv",
    ]

    auc_file = None
    for pattern in auc_patterns:
        if pattern.exists():
            auc_file = pattern
            break

    if auc_file is None:
        logging.warning(f"AUC file not found for {cluster}")
        return None

    # Load AUC matrix
    auc_df = pd.read_csv(auc_file, index_col=0)
    logging.info(f"Loaded AUC matrix: {auc_df.shape} from {auc_file}")

    # Check if index already contains full barcodes or Cell_X IDs
    first_id = auc_df.index[0]
    if first_id.startswith('Cell_'):
        # Need to map Cell_0 -> barcode
        logging.info("Cell IDs are Cell_X format, mapping to barcodes...")
        cell_file = get_cell_file_for_cluster(cell_type, cluster, config)

        if cell_file:
            mapping = load_cell_mapping_from_file(cell_file)

            if mapping is not None:
                # Verify mapping size matches AUC matrix
                if len(mapping) != len(auc_df):
                    logging.warning(f"Mismatch! Cell file has {len(mapping)} cells, AUC has {len(auc_df)}")
                    logging.warning("Cells may have been filtered during SCENIC. Proceeding with available cells.")

                # Map Cell_0 -> barcode
                auc_df = auc_df.reset_index()
                auc_df = auc_df.rename(columns={'index': 'scenic_id'})
                auc_df = auc_df.merge(mapping, on='scenic_id', how='left')

                # Check for unmapped cells
                unmapped = auc_df['cell_barcode'].isna().sum()
                if unmapped > 0:
                    logging.warning(f"{unmapped} cells could not be mapped to barcodes")

                auc_df = auc_df.dropna(subset=['cell_barcode'])
                auc_df = auc_df.set_index('cell_barcode')
                auc_df = auc_df.drop(columns=['scenic_id'])
                logging.info(f"Mapped to barcodes: {auc_df.shape}")
        else:
            logging.warning(f"No cell file found for {cluster}, keeping Cell_X IDs")
    else:
        # Index already contains full barcodes - use directly
        logging.info(f"Cell IDs are already full barcodes (e.g., {first_id[:50]}...)")

    return auc_df


def load_metadata(metadata_file):
    """Load cell metadata with condition information"""
    if not metadata_file.exists():
        logging.error(f"Metadata file not found: {metadata_file}")
        return None

    metadata = pd.read_csv(metadata_file)

    # Set cell_barcode as index if it's a column
    if 'cell_barcode' in metadata.columns:
        metadata = metadata.set_index('cell_barcode')

    logging.info(f"Loaded metadata: {metadata.shape}")
    logging.info(f"Conditions: {metadata['condition'].unique()}")

    return metadata


def merge_auc_with_conditions(auc_df, metadata):
    """Merge AUC scores with condition metadata"""

    # Find common cells
    common_cells = auc_df.index.intersection(metadata.index)

    if len(common_cells) == 0:
        logging.warning("No common cells found between AUC and metadata!")
        logging.warning(f"AUC index examples: {list(auc_df.index[:5])}")
        logging.warning(f"Metadata index examples: {list(metadata.index[:5])}")
        return None

    logging.info(f"Common cells: {len(common_cells)} / {len(auc_df)} ({100*len(common_cells)/len(auc_df):.1f}%)")

    # Merge
    merged = auc_df.loc[common_cells].copy()
    merged['condition'] = metadata.loc[common_cells, 'condition']
    merged['genotype'] = metadata.loc[common_cells, 'genotype']
    merged['Sex'] = metadata.loc[common_cells, 'Sex']

    # Add infection status if available
    if 'infection_status' in metadata.columns:
        merged['infection_status'] = metadata.loc[common_cells, 'infection_status']

    return merged


def compare_regulon_activity_by_condition(merged_df, regulons=None, output_dir=None):
    """
    Compare regulon activity across conditions

    Returns: DataFrame with mean activity per condition
    """
    if regulons is None:
        # Use all regulon columns (exclude metadata columns)
        meta_cols = ['condition', 'genotype', 'Sex', 'cluster', 'infection_status',
                     'cell_type', 'macrophage_type', 'eos_cluster', 'cell_barcode']
        regulons = [c for c in merged_df.columns if c not in meta_cols]

    # Calculate mean activity per condition
    condition_means = merged_df.groupby('condition')[regulons].mean()

    # Reorder conditions
    ordered_conditions = [c for c in CONDITION_ORDER if c in condition_means.index]
    condition_means = condition_means.loc[ordered_conditions]

    if output_dir:
        condition_means.to_csv(output_dir / 'regulon_mean_activity_by_condition.csv')

    return condition_means


def statistical_tests_by_condition(merged_df, regulons=None, comparison='WT_vs_KO'):
    """
    Perform statistical tests comparing conditions

    Comparisons:
    - WT_vs_KO: Compare WT vs KO (all conditions pooled)
    - Male_vs_Female: Compare Male vs Female
    - Infected_vs_Uninfected: Compare infected vs uninfected
    """
    if regulons is None:
        meta_cols = ['condition', 'genotype', 'Sex', 'cluster', 'infection_status',
                     'cell_type', 'macrophage_type', 'eos_cluster', 'cell_barcode']
        regulons = [c for c in merged_df.columns if c not in meta_cols]

    results = []

    for regulon in regulons:
        if comparison == 'WT_vs_KO':
            group1 = merged_df[merged_df['genotype'] == 'WT'][regulon]
            group2 = merged_df[merged_df['genotype'] == 'KO'][regulon]
            group1_name, group2_name = 'WT', 'KO'
        elif comparison == 'Male_vs_Female':
            group1 = merged_df[merged_df['Sex'] == 'Male'][regulon]
            group2 = merged_df[merged_df['Sex'] == 'Female'][regulon]
            group1_name, group2_name = 'Male', 'Female'
        else:
            continue

        if len(group1) < 3 or len(group2) < 3:
            continue

        # Wilcoxon rank-sum test
        try:
            stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')

            results.append({
                'regulon': regulon,
                'comparison': comparison,
                f'{group1_name}_mean': group1.mean(),
                f'{group2_name}_mean': group2.mean(),
                'log2FC': np.log2(group1.mean() / group2.mean()) if group2.mean() > 0 else np.nan,
                'statistic': stat,
                'pvalue': pval,
                f'{group1_name}_n': len(group1),
                f'{group2_name}_n': len(group2)
            })
        except Exception as e:
            logging.warning(f"Test failed for {regulon}: {e}")

    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        # Adjust p-values (Benjamini-Hochberg)
        from statsmodels.stats.multitest import multipletests
        _, pvals_adj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df['padj'] = pvals_adj
        results_df['significant'] = results_df['padj'] < 0.05
        results_df = results_df.sort_values('pvalue')

    return results_df


def create_condition_heatmap(condition_means, output_path, title="Regulon Activity by Condition"):
    """Create heatmap of regulon activity by condition"""

    # Z-score normalize across conditions
    zscore_df = (condition_means - condition_means.mean()) / condition_means.std()

    # Select top variable regulons
    variance = zscore_df.var()
    top_regulons = variance.nlargest(30).index.tolist()

    plot_data = zscore_df[top_regulons].T

    fig, ax = plt.subplots(figsize=(10/2.54, 12/2.54))

    sns.heatmap(plot_data, cmap='RdBu_r', center=0, ax=ax,
                xticklabels=True, yticklabels=True,
                cbar_kws={'shrink': 0.5, 'label': 'Z-score'})

    ax.set_xlabel('Condition', fontsize=8)
    ax.set_ylabel('Regulon', fontsize=8)
    ax.set_title(title, fontsize=10, fontweight='bold')

    plt.xticks(fontsize=6, rotation=45, ha='right')
    plt.yticks(fontsize=6)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"Saved heatmap: {output_path}")


def create_volcano_plot(results_df, output_path, title="Differential Regulon Activity"):
    """Create volcano plot of differential regulon activity"""

    if len(results_df) == 0:
        return

    fig, ax = plt.subplots(figsize=(8/2.54, 8/2.54))

    # -log10(pvalue) vs log2FC
    results_df['neglog10p'] = -np.log10(results_df['pvalue'])

    # Color by significance
    colors = ['red' if sig else 'gray' for sig in results_df['significant']]

    ax.scatter(results_df['log2FC'], results_df['neglog10p'],
               c=colors, alpha=0.6, s=20)

    # Add labels for significant regulons
    sig_regulons = results_df[results_df['significant']]
    for _, row in sig_regulons.iterrows():
        ax.annotate(row['regulon'],
                    (row['log2FC'], row['neglog10p']),
                    fontsize=5, alpha=0.8)

    ax.axhline(-np.log10(0.05), ls='--', color='gray', alpha=0.5)
    ax.axvline(0, ls='--', color='gray', alpha=0.5)

    ax.set_xlabel('log2(Fold Change)', fontsize=8)
    ax.set_ylabel('-log10(p-value)', fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"Saved volcano: {output_path}")


# =====================================
# MAIN ANALYSIS
# =====================================

def main():
    import argparse

    parser = argparse.ArgumentParser(description='SCENIC Condition-Specific Analysis')
    parser.add_argument('--cell-type', type=str, default='eos',
                        choices=['eos', 'macrophages'],
                        help='Cell type to analyze')
    parser.add_argument('--scenic-dir', type=str, default=str(SCENIC_DIR))
    parser.add_argument('--metadata', type=str, default=str(SEURAT_METADATA))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    metadata_file = Path(args.metadata)
    output_dir = Path(args.output_dir) / args.cell_type
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("=" * 70)
    logging.info(f"SCENIC CONDITION ANALYSIS: {args.cell_type.upper()}")
    logging.info("=" * 70)

    # Load metadata
    logging.info("\n=== Loading metadata ===")
    metadata = load_metadata(metadata_file)
    if metadata is None:
        logging.error("Failed to load metadata!")
        return 1

    # Get cell type config
    config = CELL_TYPE_CONFIG[args.cell_type]

    # Process each cluster
    all_merged = []

    for cluster in config['clusters']:
        logging.info(f"\n=== Processing {cluster} ===")

        # Load AUC matrix
        auc_df = load_auc_with_barcodes(
            scenic_dir,
            args.cell_type,
            cluster,
            config
        )

        if auc_df is None:
            continue

        # Merge with conditions
        merged = merge_auc_with_conditions(auc_df, metadata)
        if merged is not None:
            merged['cluster'] = cluster
            all_merged.append(merged)

    if not all_merged:
        logging.error("No data merged successfully!")
        return 1

    # Combine all clusters
    # Note: Different subtypes may have different regulons detected
    # pd.concat will create NaN for regulons not found in a subtype
    combined_df = pd.concat(all_merged)
    logging.info(f"\nCombined data: {combined_df.shape}")

    # Report on regulon overlap between subtypes
    meta_cols = ['condition', 'genotype', 'Sex', 'cluster', 'infection_status',
                 'cell_type', 'macrophage_type', 'eos_cluster', 'cell_barcode']
    all_regulons = set()
    subtype_regulons = {}
    for merged_df in all_merged:
        cluster_name = merged_df['cluster'].iloc[0]
        regulons = set([c for c in merged_df.columns if c not in meta_cols])
        subtype_regulons[cluster_name] = regulons
        all_regulons.update(regulons)

    logging.info(f"\nTotal unique regulons across all subtypes: {len(all_regulons)}")
    for cluster_name, regulons in subtype_regulons.items():
        logging.info(f"  {cluster_name}: {len(regulons)} regulons")

    # Find shared regulons (present in all subtypes)
    if len(subtype_regulons) > 1:
        shared_regulons = set.intersection(*subtype_regulons.values())
        logging.info(f"\nShared regulons (in all subtypes): {len(shared_regulons)}")

    # Save combined data
    combined_df.to_csv(output_dir / 'auc_with_conditions.csv')

    # =========================================
    # ANALYSIS 1: Mean activity by condition
    # =========================================
    logging.info("\n=== Analysis 1: Mean activity by condition ===")

    condition_means = compare_regulon_activity_by_condition(
        combined_df, output_dir=output_dir
    )

    create_condition_heatmap(
        condition_means,
        output_dir / 'regulon_activity_heatmap.pdf',
        title=f'{args.cell_type.title()} Regulon Activity by Condition'
    )

    # =========================================
    # ANALYSIS 1b: Mean activity by subtype AND condition
    # =========================================
    if len(subtype_regulons) > 1:
        logging.info("\n=== Analysis 1b: Mean activity by subtype and condition ===")

        # Group by both cluster and condition
        subtype_condition_means = combined_df.groupby(['cluster', 'condition'])[list(all_regulons)].mean()
        subtype_condition_means.to_csv(output_dir / 'regulon_mean_activity_by_subtype_condition.csv')
        logging.info(f"Saved subtype x condition means: {subtype_condition_means.shape}")

    # =========================================
    # ANALYSIS 2: Statistical tests
    # =========================================
    logging.info("\n=== Analysis 2: Statistical tests ===")

    # WT vs KO
    wt_ko_results = statistical_tests_by_condition(combined_df, comparison='WT_vs_KO')
    wt_ko_results.to_csv(output_dir / 'differential_regulons_WT_vs_KO.csv', index=False)

    n_sig = wt_ko_results['significant'].sum() if len(wt_ko_results) > 0 else 0
    logging.info(f"WT vs KO: {n_sig} significant regulons")

    if len(wt_ko_results) > 0:
        create_volcano_plot(
            wt_ko_results,
            output_dir / 'volcano_WT_vs_KO.pdf',
            title=f'{args.cell_type.title()}: WT vs KO'
        )

    # Male vs Female
    mf_results = statistical_tests_by_condition(combined_df, comparison='Male_vs_Female')
    mf_results.to_csv(output_dir / 'differential_regulons_Male_vs_Female.csv', index=False)

    n_sig = mf_results['significant'].sum() if len(mf_results) > 0 else 0
    logging.info(f"Male vs Female: {n_sig} significant regulons")

    if len(mf_results) > 0:
        create_volcano_plot(
            mf_results,
            output_dir / 'volcano_Male_vs_Female.pdf',
            title=f'{args.cell_type.title()}: Male vs Female'
        )

    # =========================================
    # SUMMARY
    # =========================================
    logging.info("\n" + "=" * 70)
    logging.info("ANALYSIS COMPLETE")
    logging.info("=" * 70)

    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    print("  - auc_with_conditions.csv")
    print("  - regulon_mean_activity_by_condition.csv")
    print("  - regulon_activity_heatmap.pdf/png")
    print("  - differential_regulons_WT_vs_KO.csv")
    print("  - differential_regulons_Male_vs_Female.csv")
    print("  - volcano_WT_vs_KO.pdf/png")
    print("  - volcano_Male_vs_Female.pdf/png")

    return 0


if __name__ == "__main__":
    exit(main())
