#!/usr/bin/env python3
"""
SCENIC Condition-Specific Comparison with Z-score Subsampling Method
=====================================================================
Author: Xinru Qiu
Date: December 2024

Comparisons:
1. WT_F vs KO_F (not infected females)
2. WT_M vs KO_M (males)
3. WT_F vs WT_M (sex effect in WT)
4. KO_F vs KO_M (sex effect in KO)
5. Infected_F vs Not_Infected_F (infection effect in females)

Method based on:
- 1202_Eos_vs_Mac_Subcluster_Unique_Zscore.py
- 1202_Mac_vs_Eos_Subtype_Unique_Zscore.py

Uses random subsampling to normalize cell counts and Z-score scaling for visualization.
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
# Note: actual subsample size is determined by minimum group size
RANDOM_SEED = 42

# Cutoff thresholds (adapted from Eos vs Mac scripts)
CUTOFF_LEVELS = {
    'relax': {
        'min_auc': 0.10,
        'min_fold_change': 1.2,
        'min_diff': 0.02,
    },
    'moderate': {
        'min_auc': 0.15,
        'min_fold_change': 1.3,
        'min_diff': 0.03,
    },
    'stringent': {
        'min_auc': 0.20,
        'min_fold_change': 1.5,
        'min_diff': 0.05,
    }
}

# =====================================
# CONFIGURATION
# =====================================

CONDITION_AUC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/condition_analysis")
OUTPUT_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/condition_comparison_zscore")

# Condition comparisons to perform
COMPARISONS = {
    'WT_F_vs_KO_F': {
        'group1_name': 'WT_Female',
        'group2_name': 'KO_Female',
        'group1_conditions': ['WT_Female_HFD'],
        'group2_conditions': ['KO_Female_HFD'],
        'description': 'WT vs KO in Females (not infected)'
    },
    'WT_M_vs_KO_M': {
        'group1_name': 'WT_Male',
        'group2_name': 'KO_Male',
        'group1_conditions': ['WT_Male_HFD'],
        'group2_conditions': ['KO_Male_HFD'],
        'description': 'WT vs KO in Males'
    },
    'WT_F_vs_WT_M': {
        'group1_name': 'WT_Female',
        'group2_name': 'WT_Male',
        'group1_conditions': ['WT_Female_HFD'],
        'group2_conditions': ['WT_Male_HFD'],
        'description': 'Female vs Male in WT (not infected)'
    },
    'KO_F_vs_KO_M': {
        'group1_name': 'KO_Female',
        'group2_name': 'KO_Male',
        'group1_conditions': ['KO_Female_HFD'],
        'group2_conditions': ['KO_Male_HFD'],
        'description': 'Female vs Male in KO (not infected)'
    },
    'Infected_vs_NotInfected_F': {
        'group1_name': 'Infected_Female',
        'group2_name': 'NotInfected_Female',
        'group1_conditions': ['WT_Female_HFDInfected', 'KO_Female_HFDInfected'],
        'group2_conditions': ['WT_Female_HFD', 'KO_Female_HFD'],
        'description': 'Infected vs Not Infected in Females'
    }
}


# =====================================
# FUNCTIONS
# =====================================

def load_auc_with_conditions(cell_type):
    """Load AUC matrix with condition metadata"""
    auc_file = CONDITION_AUC_DIR / cell_type / 'auc_with_conditions.csv'

    if not auc_file.exists():
        logging.error(f"AUC file not found: {auc_file}")
        return None

    df = pd.read_csv(auc_file, index_col=0)
    logging.info(f"Loaded {cell_type}: {df.shape}")
    logging.info(f"Conditions: {df['condition'].unique()}")

    return df


def get_regulon_columns(df):
    """Get regulon column names (exclude metadata columns)"""
    meta_cols = ['condition', 'genotype', 'Sex', 'cluster', 'infection_status',
                 'cell_type', 'macrophage_type', 'eos_cluster', 'cell_barcode']
    return [c for c in df.columns if c not in meta_cols]


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_mean_auc_subsampling(df, comparison_config, max_cells_per_group=None):
    """
    Calculate mean AUC using random subsampling to normalize cell counts.

    The subsample size is determined by the minimum cell count between the two groups,
    ensuring fair comparison when groups are imbalanced.

    Args:
        df: DataFrame with AUC values and condition metadata
        comparison_config: Configuration dict with group names and conditions
        max_cells_per_group: Optional maximum cells to sample (if min group size is larger)

    Returns mean AUC for each regulon in each group.
    """
    logging.info(f"\n=== Random Subsampling ===")

    np.random.seed(RANDOM_SEED)

    group1_name = comparison_config['group1_name']
    group2_name = comparison_config['group2_name']
    group1_conditions = comparison_config['group1_conditions']
    group2_conditions = comparison_config['group2_conditions']

    # Filter cells by condition
    group1_mask = df['condition'].isin(group1_conditions)
    group2_mask = df['condition'].isin(group2_conditions)

    group1_df = df[group1_mask]
    group2_df = df[group2_mask]

    n_group1 = len(group1_df)
    n_group2 = len(group2_df)

    logging.info(f"{group1_name}: {n_group1} cells")
    logging.info(f"{group2_name}: {n_group2} cells")

    if n_group1 == 0 or n_group2 == 0:
        logging.warning("Empty group, cannot proceed")
        return None, None, {}

    # Determine subsample size based on minimum group count
    subsample_size = min(n_group1, n_group2)

    # Apply optional maximum limit
    if max_cells_per_group is not None and subsample_size > max_cells_per_group:
        subsample_size = max_cells_per_group

    logging.info(f"Subsample size (min of groups): {subsample_size} cells")

    # Get regulon columns
    regulons = get_regulon_columns(df)

    # Subsample group 1
    sampled_idx_1 = np.random.choice(n_group1, subsample_size, replace=False)
    sampled_g1 = group1_df.iloc[sampled_idx_1]
    logging.info(f"Subsampled {group1_name}: {len(sampled_idx_1)} cells")

    # Subsample group 2
    sampled_idx_2 = np.random.choice(n_group2, subsample_size, replace=False)
    sampled_g2 = group2_df.iloc[sampled_idx_2]
    logging.info(f"Subsampled {group2_name}: {len(sampled_idx_2)} cells")

    # Calculate mean AUC for each regulon
    group1_auc = {}
    group2_auc = {}

    for regulon in regulons:
        group1_auc[regulon] = sampled_g1[regulon].mean()
        group2_auc[regulon] = sampled_g2[regulon].mean()

    cell_counts = {
        group1_name: n_group1,
        group2_name: n_group2,
        'subsample_size': subsample_size
    }

    return group1_auc, group2_auc, cell_counts


def find_enriched_regulons(group1_auc, group2_auc, group1_name, group2_name, cutoff_level='relax'):
    """
    Find regulons enriched in each group.

    Returns:
        results_df: All regulons with metrics
        group1_enriched: Regulons enriched in group1
        group2_enriched: Regulons enriched in group2
    """
    criteria = CUTOFF_LEVELS[cutoff_level]
    min_auc = criteria['min_auc']
    min_fc = criteria['min_fold_change']
    min_diff = criteria['min_diff']

    results = []

    for regulon in group1_auc.keys():
        g1_auc = group1_auc[regulon]
        g2_auc = group2_auc[regulon]

        # Calculate metrics
        diff = g1_auc - g2_auc

        if g2_auc > 0:
            fold_change_g1 = g1_auc / g2_auc
        else:
            fold_change_g1 = g1_auc / 0.001 if g1_auc > 0 else 0

        if g1_auc > 0:
            fold_change_g2 = g2_auc / g1_auc
        else:
            fold_change_g2 = g2_auc / 0.001 if g2_auc > 0 else 0

        # Log2 fold change (group1 / group2)
        if g2_auc > 0 and g1_auc > 0:
            log2fc = np.log2(g1_auc / g2_auc)
        else:
            log2fc = np.nan

        results.append({
            'tf': extract_tf_name(regulon),
            'regulon': regulon,
            f'{group1_name}_mean': g1_auc,
            f'{group2_name}_mean': g2_auc,
            'diff': diff,
            'fold_change_g1': fold_change_g1,
            'fold_change_g2': fold_change_g2,
            'log2FC': log2fc,
            'abs_diff': abs(diff)
        })

    results_df = pd.DataFrame(results)

    # Filter for group1-enriched regulons
    group1_enriched = results_df[
        (results_df[f'{group1_name}_mean'] >= min_auc) &
        (results_df['fold_change_g1'] >= min_fc) &
        (results_df['diff'] >= min_diff)
    ].copy()
    group1_enriched['enrichment'] = f'{group1_name}_enriched'
    group1_enriched = group1_enriched.sort_values('diff', ascending=False)

    # Filter for group2-enriched regulons
    group2_enriched = results_df[
        (results_df[f'{group2_name}_mean'] >= min_auc) &
        (results_df['fold_change_g2'] >= min_fc) &
        (results_df['diff'] <= -min_diff)
    ].copy()
    group2_enriched['enrichment'] = f'{group2_name}_enriched'
    group2_enriched = group2_enriched.sort_values('diff', ascending=True)

    logging.info(f"\n{cutoff_level.upper()} cutoff:")
    logging.info(f"  {group1_name}-enriched: {len(group1_enriched)} regulons")
    logging.info(f"  {group2_name}-enriched: {len(group2_enriched)} regulons")

    return results_df, group1_enriched, group2_enriched


def create_zscore_heatmap(group1_enriched, group2_enriched, group1_name, group2_name,
                          output_path, title, zscore_cap=2.5, max_regulons=50):
    """Create heatmap with Z-score scaled values and blue-to-red colormap"""

    # Combine enriched regulons
    combined = pd.concat([group1_enriched, group2_enriched])

    if len(combined) == 0:
        logging.warning(f"No regulons to plot for {title}")
        return

    # Limit to top regulons if too many
    if len(combined) > max_regulons:
        # Take top from each group
        n_each = max_regulons // 2
        g1_top = group1_enriched.head(n_each)
        g2_top = group2_enriched.head(n_each)
        combined = pd.concat([g1_top, g2_top])
        logging.info(f"Limiting to top {max_regulons} regulons ({n_each} per group)")

    # Build plot data (raw AUC values)
    cols = [f'{group1_name}_mean', f'{group2_name}_mean']
    raw_data = combined[cols].values

    # Calculate row-wise z-score
    row_means = raw_data.mean(axis=1, keepdims=True)
    row_stds = raw_data.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1  # Avoid division by zero
    zscore_data = (raw_data - row_means) / row_stds

    # Cap z-scores for better visualization
    zscore_data_capped = np.clip(zscore_data, -zscore_cap, zscore_cap)

    # Y-axis labels
    ylabels = combined['tf'].values

    # Save z-score matrix as CSV
    zscore_df = pd.DataFrame(zscore_data, index=ylabels, columns=[group1_name, group2_name])
    csv_path = str(output_path).replace('.pdf', '_matrix.csv')
    zscore_df.to_csv(csv_path)
    logging.info(f"Saved z-score matrix: {csv_path}")

    # Create figure with exact 1.1cm x 1.1cm squares
    cell_size_cm = 1.1
    cell_size_inches = cell_size_cm / 2.54
    n_rows = len(combined)
    n_cols = 2  # Two groups

    # Calculate exact heatmap dimensions
    heatmap_width = n_cols * cell_size_inches
    heatmap_height = n_rows * cell_size_inches

    # Add margins for labels, title, colorbar
    left_margin = 2.0  # inches for y-axis labels
    right_margin = 1.2  # inches for colorbar
    top_margin = 1.0   # inches for title
    bottom_margin = 1.0  # inches for x-axis labels

    fig_width = left_margin + heatmap_width + right_margin
    fig_height = top_margin + heatmap_height + bottom_margin

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set axis position to exact location for correct cell sizes
    ax_left = left_margin / fig_width
    ax_bottom = bottom_margin / fig_height
    ax_width = heatmap_width / fig_width
    ax_height = heatmap_height / fig_height
    ax.set_position([ax_left, ax_bottom, ax_width, ax_height])

    # Plot with z-score and blue-to-red colormap
    im = ax.imshow(zscore_data_capped, cmap=ZSCORE_CMAP, aspect='equal',
                   vmin=-zscore_cap, vmax=zscore_cap)

    # Add colorbar
    cbar_ax = fig.add_axes([ax_left + ax_width + 0.02, ax_bottom, 0.03, ax_height])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('Z-score (row-scaled)', fontsize=10, fontweight='bold')

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([group1_name, group2_name], rotation=45, ha='right', fontsize=11, fontweight='bold')

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=8)

    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_ylabel('Regulon', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_volcano_plot(results_df, group1_name, group2_name, output_path,
                        comparison_name, cutoff_level='moderate'):
    """Create volcano plot showing fold change vs differential expression"""

    criteria = CUTOFF_LEVELS[cutoff_level]
    min_auc = criteria['min_auc']
    min_fc = criteria['min_fold_change']
    min_diff = criteria['min_diff']

    fig, ax = plt.subplots(figsize=(10, 8))

    # Determine significance
    g1_sig = (
        (results_df[f'{group1_name}_mean'] >= min_auc) &
        (results_df['fold_change_g1'] >= min_fc) &
        (results_df['diff'] >= min_diff)
    )
    g2_sig = (
        (results_df[f'{group2_name}_mean'] >= min_auc) &
        (results_df['fold_change_g2'] >= min_fc) &
        (results_df['diff'] <= -min_diff)
    )

    # Colors
    colors = []
    for i in range(len(results_df)):
        if g1_sig.iloc[i]:
            colors.append('#b2182b')  # red for group1
        elif g2_sig.iloc[i]:
            colors.append('#2166ac')  # blue for group2
        else:
            colors.append('gray')

    ax.scatter(results_df['log2FC'], results_df['abs_diff'],
               c=colors, alpha=0.6, s=50, edgecolors='none')

    # Add labels for significant regulons
    for i, (idx, row) in enumerate(results_df.iterrows()):
        if g1_sig.iloc[i] or g2_sig.iloc[i]:
            ax.annotate(row['tf'],
                       (row['log2FC'], row['abs_diff']),
                       fontsize=7, alpha=0.8,
                       xytext=(5, 5), textcoords='offset points')

    ax.axhline(min_diff, ls='--', color='gray', alpha=0.5)
    ax.axvline(0, ls='--', color='gray', alpha=0.5)
    ax.axvline(np.log2(min_fc), ls='--', color='lightgray', alpha=0.5)
    ax.axvline(-np.log2(min_fc), ls='--', color='lightgray', alpha=0.5)

    ax.set_xlabel(f'log2FC ({group1_name} / {group2_name})', fontsize=12)
    ax.set_ylabel('|Difference in mean AUC|', fontsize=12)
    ax.set_title(f'{comparison_name}\n{group1_name} (red) vs {group2_name} (blue) [{cutoff_level}]',
                 fontsize=12, fontweight='bold')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#b2182b', label=f'{group1_name}-enriched ({g1_sig.sum()})'),
        Patch(facecolor='#2166ac', label=f'{group2_name}-enriched ({g2_sig.sum()})'),
        Patch(facecolor='gray', label=f'Not significant ({(~(g1_sig | g2_sig)).sum()})')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_comparison_barplot(all_results, output_path, cell_type):
    """Create bar plot showing number of enriched regulons per comparison and cutoff level"""

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)

    comparisons = list(all_results.keys())
    x = np.arange(len(comparisons))
    width = 0.35

    for ax_idx, cutoff_level in enumerate(['relax', 'moderate', 'stringent']):
        ax = axes[ax_idx]

        g1_counts = []
        g2_counts = []
        g1_names = []
        g2_names = []

        for comp_name in comparisons:
            results_df, config = all_results[comp_name]
            group1_name = config['group1_name']
            group2_name = config['group2_name']

            _, g1_enriched, g2_enriched = find_enriched_regulons(
                {r: results_df.loc[results_df['regulon'] == r, f'{group1_name}_mean'].values[0]
                 for r in results_df['regulon']},
                {r: results_df.loc[results_df['regulon'] == r, f'{group2_name}_mean'].values[0]
                 for r in results_df['regulon']},
                group1_name, group2_name, cutoff_level=cutoff_level
            )

            g1_counts.append(len(g1_enriched))
            g2_counts.append(len(g2_enriched))
            g1_names.append(group1_name)
            g2_names.append(group2_name)

        bars1 = ax.bar(x - width/2, g1_counts, width, label='Group 1', color='#b2182b', alpha=0.7)
        bars2 = ax.bar(x + width/2, g2_counts, width, label='Group 2', color='#2166ac', alpha=0.7)

        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}', ha='center', va='bottom', fontsize=9)
        for bar in bars2:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}', ha='center', va='bottom', fontsize=9)

        ax.set_xlabel('Comparison', fontsize=11)
        ax.set_ylabel('Number of enriched regulons', fontsize=11)
        ax.set_title(f'{cutoff_level.upper()}', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace('_vs_', '\nvs\n') for c in comparisons],
                          fontsize=8, rotation=0)
        ax.legend()

    plt.suptitle(f'{cell_type.upper()} - Enriched Regulons per Comparison (Z-score Method)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_summary_table(all_results, output_dir, cell_type):
    """Create summary table showing counts for each comparison and cutoff level"""

    summary_rows = []

    for comparison_name, (results_df, config) in all_results.items():
        group1_name = config['group1_name']
        group2_name = config['group2_name']

        # Reconstruct AUC dicts from results_df
        group1_auc = {r: results_df.loc[results_df['regulon'] == r, f'{group1_name}_mean'].values[0]
                      for r in results_df['regulon']}
        group2_auc = {r: results_df.loc[results_df['regulon'] == r, f'{group2_name}_mean'].values[0]
                      for r in results_df['regulon']}

        for level_name in ['relax', 'moderate', 'stringent']:
            _, g1_enriched, g2_enriched = find_enriched_regulons(
                group1_auc, group2_auc, group1_name, group2_name, cutoff_level=level_name
            )

            summary_rows.append({
                'cell_type': cell_type,
                'comparison': comparison_name,
                'cutoff_level': level_name.upper(),
                f'{group1_name}_enriched': len(g1_enriched),
                f'{group2_name}_enriched': len(g2_enriched),
                'total': len(g1_enriched) + len(g2_enriched)
            })

    summary_df = pd.DataFrame(summary_rows)
    output_path = output_dir / f'{cell_type}_comparison_summary.csv'
    summary_df.to_csv(output_path, index=False)
    logging.info(f"Saved summary: {output_path}")

    print("\n" + "="*80)
    print(f"{cell_type.upper()} COMPARISON SUMMARY (Z-SCORE METHOD)")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80 + "\n")

    return summary_df


def process_cell_type(cell_type, output_base_dir, max_cells_per_group=None):
    """Process all comparisons for a cell type

    Args:
        cell_type: 'eos' or 'macrophages'
        output_base_dir: Base output directory
        max_cells_per_group: Optional maximum cells to sample per group
    """

    logging.info(f"\n{'='*60}")
    logging.info(f"Processing {cell_type.upper()}")
    logging.info(f"{'='*60}")

    # Load AUC data with conditions
    df = load_auc_with_conditions(cell_type)
    if df is None:
        return None

    output_dir = output_base_dir / cell_type
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    # Process each comparison
    for comparison_name, config in COMPARISONS.items():
        logging.info(f"\n=== {comparison_name}: {config['description']} ===")

        # Calculate mean AUC using subsampling (subsample size = min of two groups)
        group1_auc, group2_auc, cell_counts = calculate_mean_auc_subsampling(
            df, config, max_cells_per_group=max_cells_per_group
        )

        if group1_auc is None:
            continue

        group1_name = config['group1_name']
        group2_name = config['group2_name']

        # Find enriched regulons for each cutoff level
        for cutoff_level in ['relax', 'moderate', 'stringent']:
            results_df, g1_enriched, g2_enriched = find_enriched_regulons(
                group1_auc, group2_auc, group1_name, group2_name,
                cutoff_level=cutoff_level
            )

            # Save enriched regulons CSVs
            if len(g1_enriched) > 0:
                g1_enriched.to_csv(
                    output_dir / f'{comparison_name}_{group1_name}_enriched_{cutoff_level}.csv',
                    index=False
                )
            if len(g2_enriched) > 0:
                g2_enriched.to_csv(
                    output_dir / f'{comparison_name}_{group2_name}_enriched_{cutoff_level}.csv',
                    index=False
                )

            # Create heatmap
            if len(g1_enriched) > 0 or len(g2_enriched) > 0:
                create_zscore_heatmap(
                    g1_enriched, g2_enriched, group1_name, group2_name,
                    output_dir / f'{comparison_name}_heatmap_{cutoff_level}.pdf',
                    f'{config["description"]}\n({cutoff_level.upper()}: {len(g1_enriched)} + {len(g2_enriched)} regulons)'
                )

        # Save full results (using moderate cutoff for all_results)
        results_df_full, _, _ = find_enriched_regulons(
            group1_auc, group2_auc, group1_name, group2_name, cutoff_level='relax'
        )
        results_df_full.to_csv(output_dir / f'{comparison_name}_all_metrics.csv', index=False)

        all_results[comparison_name] = (results_df_full, config)

        # Create volcano plot
        create_volcano_plot(
            results_df_full, group1_name, group2_name,
            output_dir / f'{comparison_name}_volcano.pdf',
            comparison_name,
            cutoff_level='moderate'
        )

    # Create summary table
    if all_results:
        create_summary_table(all_results, output_dir, cell_type)
        create_comparison_barplot(all_results, output_dir / f'{cell_type}_barplot_summary.pdf', cell_type)

    return all_results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='SCENIC Condition Comparison with Z-score Subsampling Method'
    )
    parser.add_argument('--cell-type', type=str, default='all',
                        choices=['eos', 'macrophages', 'all'],
                        help='Cell type to analyze')
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('--max-cells-per-group', type=int, default=None,
                        help='Optional maximum cells to subsample per group (default: use min of two groups)')
    parser.add_argument('--zscore-cap', type=float, default=2.5,
                        help='Cap z-scores at +/- this value for visualization')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("SCENIC CONDITION COMPARISON WITH Z-SCORE SUBSAMPLING")
    logging.info("="*80)
    logging.info(f"Comparisons: {list(COMPARISONS.keys())}")
    if args.max_cells_per_group:
        logging.info(f"Max cells per group: {args.max_cells_per_group}")
    else:
        logging.info("Subsample size: min of two groups (adaptive)")
    logging.info(f"Z-score cap: +/- {args.zscore_cap}")

    cell_types = ['eos', 'macrophages'] if args.cell_type == 'all' else [args.cell_type]

    all_summaries = []

    for cell_type in cell_types:
        results = process_cell_type(cell_type, output_dir, max_cells_per_group=args.max_cells_per_group)
        if results:
            # Create summary for this cell type
            summary_file = output_dir / cell_type / f'{cell_type}_comparison_summary.csv'
            if summary_file.exists():
                all_summaries.append(pd.read_csv(summary_file))

    # Combine all summaries
    if all_summaries:
        combined_summary = pd.concat(all_summaries, ignore_index=True)
        combined_summary.to_csv(output_dir / 'all_comparison_summary.csv', index=False)

        print("\n" + "="*80)
        print("COMBINED SUMMARY (ALL CELL TYPES)")
        print("="*80)
        print(combined_summary.to_string(index=False))
        print("="*80)

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
