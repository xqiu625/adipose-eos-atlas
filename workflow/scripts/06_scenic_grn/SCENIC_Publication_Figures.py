#!/usr/bin/env python3
"""
Publication-Ready SCENIC Figures for RELMα Sex-Differential Analysis
=====================================================================
Figure A: Sex-Differential Regulon Heatmap
Figure B: Target Gene Sex-Difference Variance

Color scheme:
- Female: #bf70ae (pink/purple)
- Male: #8dc63f (green)
- WT: unfilled circle
- KO: filled circle
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from pathlib import Path
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Color scheme
FEMALE_COLOR = '#bf70ae'  # Pink/purple for female
MALE_COLOR = '#8dc63f'    # Green for male
WT_COLOR = '#1a9850'      # Green for WT (for bars if needed)
KO_COLOR = '#d73027'      # Red for KO (for bars if needed)

# Z-score colormap for heatmap: Blue (low, -2) - Gray (mid, 0) - Red (high, +2)
# Matching other SCENIC heatmaps
ZSCORE_CMAP = LinearSegmentedColormap.from_list('scenic_zscore', ['#2a6daf', '#f2f3f5', '#b62437'])

# Paths
SCENIC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned")
CONDITION_AUC_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/condition_analysis")
OUTPUT_DIR = Path("/bigdata/godziklab/shared/Xinru/302025/scenic_results/relma_sex_effect")
H5AD_DIR = Path("/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/scenic_project/inputs")

# Key sex-differential TFs
KEY_TFS = ['Fli1', 'Tef', 'Rara', 'Rel', 'Irf4']

# Cell size in cm for heatmap
CELL_HEIGHT_CM = 0.81
CELL_WIDTH_CM = 0.87


def load_auc_with_conditions(cell_type='eos'):
    """Load AUC matrix with condition metadata"""
    auc_file = CONDITION_AUC_DIR / cell_type / 'auc_with_conditions.csv'

    if not auc_file.exists():
        logging.error(f"AUC file not found: {auc_file}")
        return None

    df = pd.read_csv(auc_file, index_col=0)
    logging.info(f"Loaded {cell_type}: {df.shape}")
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


def calculate_condition_means(df, conditions, regulons):
    """Calculate mean AUC for each regulon in each condition"""
    means = {}
    for cond in conditions:
        cond_df = df[df['condition'] == cond]
        means[cond] = {reg: cond_df[reg].mean() for reg in regulons}
    return means


def find_sex_differential_regulons(df, min_fc=1.2, min_diff=0.02, min_auc=0.10, top_n=20):
    """
    Find regulons that show sex-differential activity in WT.

    Args:
        df: DataFrame with AUC values (already filtered by cluster if needed)
        min_fc: Minimum fold change
        min_diff: Minimum absolute difference
        min_auc: Minimum AUC threshold
        top_n: Maximum number of TFs to return

    Returns:
        List of TF names sorted by absolute difference
    """
    regulons = get_regulon_columns(df)
    conditions = ['WT_Female_HFD', 'WT_Male_HFD']

    # Calculate means for WT conditions
    means = calculate_condition_means(df, conditions, regulons)

    results = []
    for reg in regulons:
        wt_f = means['WT_Female_HFD'].get(reg, 0)
        wt_m = means['WT_Male_HFD'].get(reg, 0)

        diff = wt_f - wt_m
        abs_diff = abs(diff)

        # Calculate fold change
        if wt_m > 0 and wt_f > 0:
            fc = max(wt_f / wt_m, wt_m / wt_f)
        else:
            fc = 0

        max_auc = max(wt_f, wt_m)

        # Filter by criteria
        if max_auc >= min_auc and fc >= min_fc and abs_diff >= min_diff:
            tf_name = extract_tf_name(reg)
            results.append({
                'tf': tf_name,
                'regulon': reg,
                'wt_f': wt_f,
                'wt_m': wt_m,
                'diff': diff,
                'abs_diff': abs_diff,
                'fc': fc,
                'enrichment': 'Female' if diff > 0 else 'Male'
            })

    # Sort by absolute difference
    results.sort(key=lambda x: x['abs_diff'], reverse=True)

    # Get unique TFs (take first if multiple regulons for same TF)
    seen_tfs = set()
    unique_results = []
    for r in results:
        if r['tf'] not in seen_tfs:
            seen_tfs.add(r['tf'])
            unique_results.append(r)

    # Take top N
    top_results = unique_results[:top_n]

    logging.info(f"Found {len(top_results)} sex-differential TFs:")
    for r in top_results:
        logging.info(f"  {r['tf']}: {r['enrichment']}-enriched (diff={r['diff']:.4f}, FC={r['fc']:.2f})")

    return [r['tf'] for r in top_results], top_results


def create_figure_a_heatmap(df, output_path, tfs=KEY_TFS, cluster_filter=None, title_suffix="", auto_discover_tfs=False, top_n_tfs=20):
    """
    Create Figure A: Sex-Differential Regulon Heatmap

    Features:
    - 4 conditions: WT_F, WT_M, KO_F, KO_M
    - Side annotation bar showing Female/Male enrichment in WT
    - Square cells: 0.81cm (height) x 0.87cm (width)
    - Color scale with legend
    - New colors: Female=#bf70ae, Male=#8dc63f

    Args:
        df: DataFrame with AUC values
        output_path: Path to save the figure
        tfs: List of TFs to include (ignored if auto_discover_tfs=True)
        cluster_filter: If provided, filter to specific cluster (e.g., 'Eos_C3')
        title_suffix: Additional text to add to title
        auto_discover_tfs: If True, automatically find sex-differential TFs
        top_n_tfs: Maximum number of TFs to include when auto-discovering
    """
    cluster_info = f" ({cluster_filter})" if cluster_filter else " (All Clusters)"
    logging.info(f"Creating Figure A: Sex-Differential Regulon Heatmap{cluster_info}")

    # Filter by cluster if specified
    if cluster_filter is not None:
        if 'eos_cluster' in df.columns:
            df = df[df['eos_cluster'] == cluster_filter].copy()
            logging.info(f"Filtered to {cluster_filter}: {len(df)} cells")
        elif 'cluster' in df.columns:
            df = df[df['cluster'] == cluster_filter].copy()
            logging.info(f"Filtered to {cluster_filter}: {len(df)} cells")
        else:
            logging.warning("No cluster column found, using all cells")

    if len(df) == 0:
        logging.error("No cells after filtering")
        return

    # Auto-discover sex-differential TFs if requested
    if auto_discover_tfs:
        logging.info(f"Auto-discovering sex-differential TFs (top {top_n_tfs})...")
        tfs, tf_results = find_sex_differential_regulons(df, top_n=top_n_tfs)
        if not tfs:
            logging.warning("No sex-differential TFs found with current thresholds")
            return

        # Save the discovered TFs to CSV
        csv_path = str(output_path).replace('.pdf', '_discovered_tfs.csv')
        import pandas as pd
        pd.DataFrame(tf_results).to_csv(csv_path, index=False)
        logging.info(f"Saved discovered TFs to: {csv_path}")

    # Get regulon columns
    all_regulons = get_regulon_columns(df)

    # Find matching regulons for our TFs
    tf_regulons = {}
    for tf in tfs:
        for reg in all_regulons:
            if extract_tf_name(reg) == tf:
                tf_regulons[tf] = reg
                break

    if not tf_regulons:
        logging.error("No matching regulons found for key TFs")
        return

    logging.info(f"Found regulons: {tf_regulons}")

    # Define conditions (non-infected only)
    conditions = ['WT_Female_HFD', 'WT_Male_HFD', 'KO_Female_HFD', 'KO_Male_HFD']
    condition_labels = ['WT_F', 'WT_M', 'KO_F', 'KO_M']

    # Calculate mean AUC for each condition
    means = calculate_condition_means(df, conditions, list(tf_regulons.values()))

    # Build data matrix (TFs x conditions)
    n_tfs = len(tfs)
    n_conds = len(conditions)

    data_matrix = np.zeros((n_tfs, n_conds))
    for i, tf in enumerate(tfs):
        reg = tf_regulons.get(tf)
        if reg:
            for j, cond in enumerate(conditions):
                data_matrix[i, j] = means[cond].get(reg, 0)

    # Calculate row-wise z-scores
    row_means = data_matrix.mean(axis=1, keepdims=True)
    row_stds = data_matrix.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    zscore_matrix = (data_matrix - row_means) / row_stds

    # Cap z-scores
    zscore_cap = 2.5
    zscore_matrix = np.clip(zscore_matrix, -zscore_cap, zscore_cap)

    # Determine sex enrichment in WT (for annotation bar)
    # Positive z-score in WT_F vs WT_M means Female-enriched
    sex_enrichment = []
    for i, tf in enumerate(tfs):
        wt_f_zscore = zscore_matrix[i, 0]  # WT_F
        wt_m_zscore = zscore_matrix[i, 1]  # WT_M
        if wt_f_zscore > wt_m_zscore:
            sex_enrichment.append('Female')
        else:
            sex_enrichment.append('Male')

    # Figure dimensions based on cell size
    cell_width_in = CELL_WIDTH_CM / 2.54
    cell_height_in = CELL_HEIGHT_CM / 2.54

    heatmap_width = n_conds * cell_width_in
    heatmap_height = n_tfs * cell_height_in

    # Add margins
    left_margin = 0.8  # for TF labels
    right_margin = 2.5  # for colorbar and legend
    top_margin = 1.0   # for title
    bottom_margin = 0.8  # for condition labels
    annot_width = 0.15  # annotation bar width

    fig_width = left_margin + annot_width + 0.05 + heatmap_width + right_margin
    fig_height = top_margin + heatmap_height + bottom_margin

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create axes for annotation bar
    ax_annot_left = left_margin / fig_width
    ax_annot_bottom = bottom_margin / fig_height
    ax_annot_width = annot_width / fig_width
    ax_annot_height = heatmap_height / fig_height

    ax_annot = fig.add_axes([ax_annot_left, ax_annot_bottom, ax_annot_width, ax_annot_height])

    # Create axes for heatmap
    ax_heat_left = (left_margin + annot_width + 0.05) / fig_width
    ax_heat_width = heatmap_width / fig_width

    ax_heat = fig.add_axes([ax_heat_left, ax_annot_bottom, ax_heat_width, ax_annot_height])

    # Draw annotation bar
    annot_colors = [FEMALE_COLOR if e == 'Female' else MALE_COLOR for e in sex_enrichment]
    for i, color in enumerate(annot_colors):
        ax_annot.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5))

    ax_annot.set_xlim(0, 1)
    ax_annot.set_ylim(0, n_tfs)
    ax_annot.set_xticks([])
    ax_annot.set_yticks(np.arange(n_tfs) + 0.5)
    ax_annot.set_yticklabels(tfs, fontsize=9)
    ax_annot.invert_yaxis()
    ax_annot.spines['top'].set_visible(False)
    ax_annot.spines['right'].set_visible(False)
    ax_annot.spines['bottom'].set_visible(False)
    ax_annot.spines['left'].set_visible(False)

    # Draw heatmap
    im = ax_heat.imshow(zscore_matrix, cmap=ZSCORE_CMAP, aspect='equal',
                        vmin=-zscore_cap, vmax=zscore_cap)

    # Set ticks
    ax_heat.set_xticks(np.arange(n_conds))
    ax_heat.set_xticklabels(condition_labels, fontsize=9, fontweight='bold')
    ax_heat.set_yticks([])

    # Add group labels (WT, KO) below
    ax_heat.text(0.5, n_tfs + 0.7, 'WT', ha='center', va='top', fontsize=10, fontweight='bold',
                 transform=ax_heat.get_xaxis_transform())
    ax_heat.text(2.5, n_tfs + 0.7, 'KO', ha='center', va='top', fontsize=10, fontweight='bold',
                 transform=ax_heat.get_xaxis_transform())

    # Add vertical line between WT and KO
    ax_heat.axvline(x=1.5, color='black', linewidth=1.5)

    # Add grid lines
    for i in range(n_tfs + 1):
        ax_heat.axhline(y=i - 0.5, color='white', linewidth=0.5)
    for j in range(n_conds + 1):
        ax_heat.axvline(x=j - 0.5, color='white', linewidth=0.5)

    # Add colorbar
    cbar_left = ax_heat_left + ax_heat_width + 0.02
    cbar_ax = fig.add_axes([cbar_left, ax_annot_bottom + 0.1, 0.02, ax_annot_height - 0.2])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('Z-score', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Add legend for annotation colors
    legend_elements = [
        mpatches.Patch(facecolor=FEMALE_COLOR, edgecolor='black', label='Female-enriched in WT'),
        mpatches.Patch(facecolor=MALE_COLOR, edgecolor='black', label='Male-enriched in WT')
    ]

    legend_ax = fig.add_axes([cbar_left + 0.08, ax_annot_bottom, 0.3, ax_annot_height])
    legend_ax.axis('off')
    legend_ax.legend(handles=legend_elements, loc='center left', frameon=True,
                     fontsize=7, title_fontsize=8)

    # Title
    # Title with optional cluster info
    title_text = f'EOS: Sex-Differential Regulons{title_suffix}\nWT vs KO Comparison'
    if cluster_filter:
        title_text = f'EOS {cluster_filter}: Sex-Differential Regulons{title_suffix}\nWT vs KO Comparison'
    fig.suptitle(title_text,
                 fontsize=11, fontweight='bold', y=0.95)

    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved Figure A: {output_path}")
    plt.close()


def load_regulon_targets(scenic_dir, cell_type='eos'):
    """Load target genes for each regulon"""
    logging.info("Loading regulon target genes...")

    tf_targets = {}
    clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']

    for cluster in clusters:
        regulon_file = scenic_dir / cell_type / cluster / 'step2_regulons.csv'
        if not regulon_file.exists():
            regulon_file = scenic_dir / cell_type / cluster / 'regulons.csv'
        if not regulon_file.exists():
            continue

        try:
            df = pd.read_csv(regulon_file)
            for _, row in df.iterrows():
                tf = row['tf'] if 'tf' in row else row.iloc[0]
                if pd.isna(tf):
                    continue

                targets_str = row['target_genes'] if 'target_genes' in row else row.iloc[1]
                if pd.notna(targets_str) and str(targets_str) != 'nan':
                    genes = [g.strip() for g in str(targets_str).split(';') if g.strip()]
                    genes = [g for g in genes if not g.startswith('e-') and not g.startswith('e+')]
                    if genes:
                        if tf not in tf_targets:
                            tf_targets[tf] = set()
                        tf_targets[tf].update(genes)
        except Exception as e:
            logging.warning(f"Error loading {regulon_file}: {e}")
            continue

    tf_targets = {k: list(v) for k, v in tf_targets.items()}
    logging.info(f"Loaded targets for {len(tf_targets)} TFs")
    return tf_targets


def load_gene_expression(h5ad_dir, cell_type='eos'):
    """Load gene expression data from h5ad file"""
    import scanpy as sc

    if cell_type == 'eos':
        h5ad_file = h5ad_dir / 'eos_all.h5ad'
    else:
        h5ad_file = h5ad_dir / 'mac_all.h5ad'

    logging.info(f"Loading expression from: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)
    logging.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def calculate_gene_sex_differences(adata, genes, genotype):
    """
    Calculate sex difference (Female - Male) for each gene for a given genotype.
    Returns dict: {gene: {'diff': value, 'female_vals': [...], 'male_vals': [...], 'sem': value}}
    """
    import scipy.sparse as sp

    f_mask = (adata.obs['genotype'] == genotype) & (adata.obs['Sex'] == 'Female')
    m_mask = (adata.obs['genotype'] == genotype) & (adata.obs['Sex'] == 'Male')

    f_mask = f_mask.values
    m_mask = m_mask.values

    gene_data = {}
    gene_names = list(adata.var_names)

    for gene in genes:
        if gene in gene_names:
            gene_idx = gene_names.index(gene)

            if sp.issparse(adata.X):
                f_expr = adata.X[f_mask, gene_idx].toarray().flatten()
                m_expr = adata.X[m_mask, gene_idx].toarray().flatten()
            else:
                f_expr = adata.X[f_mask, gene_idx]
                m_expr = adata.X[m_mask, gene_idx]

            f_mean = np.mean(f_expr) if len(f_expr) > 0 else 0
            m_mean = np.mean(m_expr) if len(m_expr) > 0 else 0

            gene_data[gene] = {
                'diff': f_mean - m_mean,
                'female_mean': f_mean,
                'male_mean': m_mean,
                'female_vals': f_expr,
                'male_vals': m_expr
            }

    return gene_data


def create_figure_b_variance(tf_targets, wt_gene_data, ko_gene_data, output_path, tfs=KEY_TFS):
    """
    Create Figure B: Target Gene Sex-Difference Variance

    Features:
    - Bar plot of variance with symbols for genotype
    - WT = unfilled circle, KO = filled circle
    - Female/Male colors: #bf70ae / #8dc63f
    - Individual data points
    - SEM error bars
    """
    logging.info("Creating Figure B: Target Gene Sex-Difference Variance")

    # Calculate variance of sex differences for each TF's targets
    wt_variances = []
    ko_variances = []
    wt_diffs_all = []  # All individual differences for each TF
    ko_diffs_all = []
    tf_labels = []

    for tf in tfs:
        targets = tf_targets.get(tf, [])
        if not targets:
            continue

        wt_diffs = [wt_gene_data.get(g, {}).get('diff', 0) for g in targets if g in wt_gene_data]
        ko_diffs = [ko_gene_data.get(g, {}).get('diff', 0) for g in targets if g in ko_gene_data]

        if wt_diffs and ko_diffs:
            wt_variances.append(np.var(wt_diffs))
            ko_variances.append(np.var(ko_diffs))
            wt_diffs_all.append(wt_diffs)
            ko_diffs_all.append(ko_diffs)
            tf_labels.append(tf)

    if not tf_labels:
        logging.warning("No data for Figure B")
        return

    n_tfs = len(tf_labels)

    # Calculate SEM for variances (using bootstrap or standard error of variance)
    wt_sem = []
    ko_sem = []
    for i in range(n_tfs):
        # SEM of variance can be approximated as sqrt(2*var^2 / (n-1))
        n_wt = len(wt_diffs_all[i])
        n_ko = len(ko_diffs_all[i])
        wt_sem.append(np.sqrt(2 * wt_variances[i]**2 / max(n_wt - 1, 1)))
        ko_sem.append(np.sqrt(2 * ko_variances[i]**2 / max(n_ko - 1, 1)))

    # Figure size
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    x = np.arange(n_tfs)
    width = 0.35

    # Panel 1: Variance bar plot with symbols
    ax1 = axes[0]

    # Plot bars
    bars_wt = ax1.bar(x - width/2, wt_variances, width, color='white', edgecolor='black', linewidth=1.5, label='WT')
    bars_ko = ax1.bar(x + width/2, ko_variances, width, color='gray', edgecolor='black', linewidth=1.5, label='KO')

    # Add error bars
    ax1.errorbar(x - width/2, wt_variances, yerr=wt_sem, fmt='none', ecolor='black', capsize=3, capthick=1)
    ax1.errorbar(x + width/2, ko_variances, yerr=ko_sem, fmt='none', ecolor='black', capsize=3, capthick=1)

    # Add individual data points (scatter of variance estimates from subsets)
    # For individual points, we'll bootstrap variance estimates
    np.random.seed(42)
    for i in range(n_tfs):
        # Bootstrap variance for WT
        wt_vals = wt_diffs_all[i]
        ko_vals = ko_diffs_all[i]

        n_bootstrap = 20
        wt_boot_vars = []
        ko_boot_vars = []

        for _ in range(n_bootstrap):
            if len(wt_vals) > 1:
                sample_wt = np.random.choice(wt_vals, size=len(wt_vals), replace=True)
                wt_boot_vars.append(np.var(sample_wt))
            if len(ko_vals) > 1:
                sample_ko = np.random.choice(ko_vals, size=len(ko_vals), replace=True)
                ko_boot_vars.append(np.var(sample_ko))

        # Plot individual points
        if wt_boot_vars:
            ax1.scatter([x[i] - width/2] * len(wt_boot_vars), wt_boot_vars,
                       marker='o', facecolors='none', edgecolors='black', s=15, alpha=0.5, zorder=3)
        if ko_boot_vars:
            ax1.scatter([x[i] + width/2] * len(ko_boot_vars), ko_boot_vars,
                       marker='o', facecolors='black', edgecolors='black', s=15, alpha=0.5, zorder=3)

    ax1.set_xticks(x)
    ax1.set_xticklabels(tf_labels, fontsize=9)
    ax1.set_ylabel('Variance of sex diff', fontsize=10)
    ax1.set_title('Target gene sex-difference variance', fontsize=11, fontweight='bold')

    # Custom legend with symbols
    legend_elements = [
        Line2D([0], [0], marker='o', color='white', markerfacecolor='white',
               markeredgecolor='black', markersize=8, label='WT', linestyle='None'),
        Line2D([0], [0], marker='o', color='white', markerfacecolor='black',
               markeredgecolor='black', markersize=8, label='KO', linestyle='None')
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=8)

    # Panel 2: Reduction percentage with sex colors
    ax2 = axes[1]

    reductions = [(wt - ko) / wt * 100 if wt > 0 else 0
                  for wt, ko in zip(wt_variances, ko_variances)]

    # Color bars based on direction: reduction = female color (good), increase = male color
    colors = [FEMALE_COLOR if r > 0 else MALE_COLOR for r in reductions]

    bars = ax2.bar(x, reductions, color=colors, edgecolor='black', linewidth=1)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

    # Add individual data points for reduction (using bootstrap)
    for i in range(n_tfs):
        wt_vals = wt_diffs_all[i]
        ko_vals = ko_diffs_all[i]

        n_bootstrap = 20
        reduction_points = []

        for _ in range(n_bootstrap):
            if len(wt_vals) > 1 and len(ko_vals) > 1:
                sample_wt = np.random.choice(wt_vals, size=len(wt_vals), replace=True)
                sample_ko = np.random.choice(ko_vals, size=len(ko_vals), replace=True)
                var_wt = np.var(sample_wt)
                var_ko = np.var(sample_ko)
                if var_wt > 0:
                    reduction_points.append((var_wt - var_ko) / var_wt * 100)

        if reduction_points:
            # Use symbol based on whether positive or negative
            for rp in reduction_points:
                marker_face = 'none' if rp > 0 else 'black'
                ax2.scatter(x[i], rp, marker='o', facecolors=marker_face,
                           edgecolors='black', s=15, alpha=0.5, zorder=3)

    # Calculate and add SEM for reductions
    reduction_sems = []
    for i in range(n_tfs):
        wt_vals = wt_diffs_all[i]
        ko_vals = ko_diffs_all[i]

        reduction_boot = []
        for _ in range(100):
            if len(wt_vals) > 1 and len(ko_vals) > 1:
                sample_wt = np.random.choice(wt_vals, size=len(wt_vals), replace=True)
                sample_ko = np.random.choice(ko_vals, size=len(ko_vals), replace=True)
                var_wt = np.var(sample_wt)
                var_ko = np.var(sample_ko)
                if var_wt > 0:
                    reduction_boot.append((var_wt - var_ko) / var_wt * 100)

        reduction_sems.append(np.std(reduction_boot) if reduction_boot else 0)

    ax2.errorbar(x, reductions, yerr=reduction_sems, fmt='none', ecolor='black', capsize=3, capthick=1)

    ax2.set_xticks(x)
    ax2.set_xticklabels(tf_labels, fontsize=9)
    ax2.set_ylabel('Reduction (%)', fontsize=10)
    ax2.set_title('Sex-difference reduction after KO', fontsize=11, fontweight='bold')

    # Add legend for colors
    legend_elements_2 = [
        mpatches.Patch(facecolor=FEMALE_COLOR, edgecolor='black', label='Reduction (Female)'),
        mpatches.Patch(facecolor=MALE_COLOR, edgecolor='black', label='Increase (Male)')
    ]
    ax2.legend(handles=legend_elements_2, loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved Figure B: {output_path}")
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Publication-Ready SCENIC Figures')
    parser.add_argument('--scenic-dir', type=str, default=str(SCENIC_DIR))
    parser.add_argument('--auc-dir', type=str, default=str(CONDITION_AUC_DIR))
    parser.add_argument('--h5ad-dir', type=str, default=str(H5AD_DIR))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('--cell-type', type=str, default='eos')
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    h5ad_dir = Path(args.h5ad_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("PUBLICATION-READY SCENIC FIGURES")
    logging.info("="*80)

    # ===================
    # Figure A: Heatmap (All Clusters)
    # ===================
    df = load_auc_with_conditions(args.cell_type)
    if df is not None:
        # All clusters
        create_figure_a_heatmap(
            df.copy(),
            output_dir / 'Figure_A_Sex_Differential_Regulon_Heatmap.pdf',
            tfs=KEY_TFS
        )

        # Eos C3 only (terminally differentiated) - with same 5 TFs for comparison
        create_figure_a_heatmap(
            df.copy(),
            output_dir / 'Figure_A_Sex_Differential_Regulon_Heatmap_C3_5TFs.pdf',
            tfs=KEY_TFS,
            cluster_filter='Eos_C3',
            title_suffix=' (Terminally Differentiated)'
        )

        # Eos C3 only - auto-discover ALL sex-differential TFs
        create_figure_a_heatmap(
            df.copy(),
            output_dir / 'Figure_A_Sex_Differential_Regulon_Heatmap_C3_AllTFs.pdf',
            cluster_filter='Eos_C3',
            title_suffix=' (Terminally Differentiated)',
            auto_discover_tfs=True,
            top_n_tfs=20
        )

    # ===================
    # Figure B: Variance
    # ===================
    # Load regulon targets
    tf_targets = load_regulon_targets(scenic_dir, args.cell_type)

    # Load expression data
    adata = load_gene_expression(h5ad_dir, args.cell_type)

    # Collect all target genes
    all_genes = set()
    for tf in KEY_TFS:
        if tf in tf_targets:
            all_genes.update(tf_targets[tf])

    # Calculate sex differences for WT and KO
    logging.info("Calculating WT sex differences...")
    wt_gene_data = calculate_gene_sex_differences(adata, list(all_genes), 'WT')

    logging.info("Calculating KO sex differences...")
    ko_gene_data = calculate_gene_sex_differences(adata, list(all_genes), 'KO')

    # Create Figure B
    create_figure_b_variance(
        tf_targets,
        {g: d for g, d in wt_gene_data.items()},
        {g: d for g, d in ko_gene_data.items()},
        output_dir / 'Figure_B_Target_Gene_Sex_Difference_Variance.pdf',
        tfs=KEY_TFS
    )

    logging.info("\n" + "="*80)
    logging.info("COMPLETE")
    logging.info("="*80)
    print(f"\nOutput: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir}" .')


if __name__ == "__main__":
    main()
