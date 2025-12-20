#!/usr/bin/env python3
"""
GRN Network Visualization for Eos-Enriched Regulons
Updated 2024-12-09:
- Manuscript colors for Eos subclusters
- Figure size: 8cm height x 9cm width
- Font: 6pt minimum
- Three visualization modes: AUC colored, Subcluster colored, Hybrid (color by cluster + size by AUC)
- Correlation thresholds: 0.3 and 0.6
- Aesthetic style matching Frontiers paper Fig 3B

Output:
- Network colored by AUC (white-red gradient, 0-1 range)
- Network colored by subcluster
- Hybrid network: colored by subcluster with node size proportional to AUC
- Both at correlation >= 0.3 and >= 0.6
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from pathlib import Path
from scipy import stats
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Manuscript Eos subcluster colors
EOS_COLORS = {
    'Eos_C0': '#6ca6cd',
    'Eos_C1': '#e84a60',
    'Eos_C2': '#f37026',
    'Eos_C3': '#ccabcb',
    'Eos_C4': '#b13b49',
    'Eos_C5': '#5ba4a4',
    'Eos_C6': '#2f3d7f',
}

# Mac subtype colors (for reference)
MAC_COLORS = {
    'M1': '#cc426c',
    'M2_Functional': '#8fc43e',
    'M2_Surface': '#58429a',
}

# AUC colormap (white to red)
AUC_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'auc', ['#f7f7f7', '#fee0d2', '#fc9272', '#de2d26']
)

# AUC scale for Eos (0 to 0.9)
AUC_VMIN = 0
AUC_VMAX = 0.9

CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42


def load_auc_matrices(scenic_dir):
    """Load AUC matrices for Eos subclusters"""
    eos_clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
    eos_data = {}

    for cluster in eos_clusters:
        auc_file = scenic_dir / "eos" / cluster / "step3_auc_matrix.csv"
        if auc_file.exists():
            eos_data[cluster] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {cluster}: {eos_data[cluster].shape[0]} cells")

    loaded_eos = [c for c in eos_clusters if c in eos_data]
    return eos_data, loaded_eos


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_regulon_coactivity(eos_data, eos_clusters, target_regulons, cells_per_cluster=200):
    """Calculate co-activity correlation between regulons across all Eos cells."""
    logging.info("Calculating regulon co-activity correlations...")
    np.random.seed(RANDOM_SEED)

    tf_to_regulon = {}
    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        for col in eos_data[cluster].columns:
            tf_name = extract_tf_name(col)
            if tf_name not in tf_to_regulon:
                tf_to_regulon[tf_name] = {}
            tf_to_regulon[tf_name][cluster] = col

    all_cells_auc = []
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

        cell_data = {}
        for tf in target_regulons:
            if tf in tf_to_regulon and cluster in tf_to_regulon[tf]:
                regulon_col = tf_to_regulon[tf][cluster]
                cell_data[tf] = sampled_df[regulon_col].values
            else:
                cell_data[tf] = np.zeros(cells_per_cluster)

        cluster_auc_df = pd.DataFrame(cell_data)
        cluster_auc_df['cluster'] = cluster
        all_cells_auc.append(cluster_auc_df)

    combined_auc = pd.concat(all_cells_auc, ignore_index=True)
    auc_matrix = combined_auc[target_regulons]
    corr_matrix = auc_matrix.corr(method='pearson')

    logging.info(f"Calculated {len(target_regulons)}x{len(target_regulons)} correlation matrix")
    return corr_matrix, combined_auc


def build_grn_network(corr_matrix, regulon_info, correlation_threshold=0.3):
    """Build NetworkX graph from correlation matrix."""
    logging.info(f"Building network with correlation threshold >= {correlation_threshold}")

    G = nx.Graph()

    for tf in corr_matrix.columns:
        if tf in regulon_info['tf'].values:
            info = regulon_info[regulon_info['tf'] == tf].iloc[0]

            # Get z-score columns (eos_Eos_C0, etc.)
            zscore_cols = [c for c in regulon_info.columns if c.startswith('eos_')]
            if zscore_cols:
                zscores = [info[c] for c in zscore_cols]
                mean_zscore = np.mean(zscores)
                max_zscore = np.max(zscores)
            else:
                mean_zscore = 0
                max_zscore = 0

            G.add_node(tf,
                       cluster=info['max_eos_cluster'],
                       auc=info['max_eos_auc'],
                       uniqueness=info['uniqueness_score'],
                       zscore=max_zscore,
                       mean_zscore=mean_zscore,
                       color=EOS_COLORS.get(info['max_eos_cluster'], '#808080'))
        else:
            G.add_node(tf, cluster='Unknown', auc=0, uniqueness=0, zscore=0, mean_zscore=0, color='#808080')

    regulons = list(corr_matrix.columns)
    n_edges = 0

    for i in range(len(regulons)):
        for j in range(i + 1, len(regulons)):
            corr = corr_matrix.iloc[i, j]
            if abs(corr) >= correlation_threshold:
                G.add_edge(regulons[i], regulons[j], weight=corr, abs_weight=abs(corr))
                n_edges += 1

    logging.info(f"Created network with {G.number_of_nodes()} nodes and {n_edges} edges")
    return G


def create_auc_network(G, output_path, title="Eos Regulon Network (AUC)",
                       show_labels=True):
    """
    Create network visualization with nodes colored by AUC value.
    AUC range: 0-1, using white-to-red gradient.

    Figure size: 8cm height x 9cm width
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Figure size: 9cm width x 8cm height
    fig_width = 9 / 2.54  # inches
    fig_height = 8 / 2.54  # inches

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate layout - use kamada_kawai for cleaner appearance
    # Use absolute weights to avoid negative weight issues
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Get AUC values and normalize for coloring (0-1 range)
    auc_values = np.array([G.nodes[n].get('auc', 0) for n in G.nodes()])

    # Normalize to AUC scale for colormap
    norm = mcolors.Normalize(vmin=AUC_VMIN, vmax=AUC_VMAX)
    node_colors = [AUC_CMAP(norm(a)) for a in auc_values]

    # Node sizes - smaller and more uniform for cleaner look
    node_sizes = [80 for _ in G.nodes()]  # Fixed size for cleaner appearance

    # Edge styling - thin gray lines
    edge_widths = [0.3 for _ in G.edges()]  # Thin lines
    edge_colors = ['#cccccc' for _ in G.edges()]  # Light gray

    # Draw edges first (behind nodes)
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.5,
                           edge_color=edge_colors)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.5)

    # Draw labels
    if show_labels:
        all_labels = {n: n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, ax=ax, labels=all_labels,
                                font_size=6, font_weight='normal')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=AUC_CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, aspect=15, pad=0.02)
    cbar.set_label('AUC', fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    ax.set_title(title, fontsize=6, pad=5)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_hybrid_network(G, output_path, title="Eos Regulon Network (Hybrid)",
                          show_labels=True, min_size=40, max_size=200):
    """
    Create network visualization with:
    - Node COLOR by subcluster
    - Node SIZE proportional to AUC value

    Figure size: 8cm height x 9cm width
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Figure size: 9cm width x 8cm height
    fig_width = 9 / 2.54  # inches
    fig_height = 8 / 2.54  # inches

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate layout
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Node colors by subcluster
    node_colors = [G.nodes[n].get('color', '#808080') for n in G.nodes()]

    # Node sizes proportional to AUC (scale from min_size to max_size)
    auc_values = np.array([G.nodes[n].get('auc', 0) for n in G.nodes()])
    # Scale AUC to size range (normalize by AUC_VMAX)
    auc_normalized = np.clip(auc_values / AUC_VMAX, 0, 1)
    node_sizes = min_size + (auc_normalized * (max_size - min_size))

    # Edge styling
    edge_widths = [0.3 for _ in G.edges()]
    edge_colors = ['#cccccc' for _ in G.edges()]

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.5,
                           edge_color=edge_colors)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.5)

    # Draw labels
    if show_labels:
        all_labels = {n: n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, ax=ax, labels=all_labels,
                                font_size=6, font_weight='normal')

    # Create legend for Eos subclusters
    legend_elements = []
    for cluster in ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']:
        color = EOS_COLORS.get(cluster, '#808080')
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            legend_elements.append(plt.scatter([], [], c=color, s=30,
                                               label=f'{cluster} ({count})'))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper left', fontsize=5,
                  title='Subcluster', title_fontsize=5, framealpha=0.8,
                  markerscale=0.8)

    # Add size legend annotation
    ax.annotate('Node size ∝ AUC', xy=(0.98, 0.02), xycoords='axes fraction',
                fontsize=5, ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.set_title(title, fontsize=6, pad=5)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_unified_network(G, output_path, title="Eos Regulon Network",
                           show_labels=True, min_size=20, max_size=300):
    """
    Create UNIFIED network visualization with:
    - Node COLOR by subcluster (using manuscript colors)
    - Node SIZE proportional to AUC value (bigger = AUC closer to 1, smaller = AUC closer to 0)
    - Proper size legend showing AUC scale

    This is the main figure combining both subcluster and AUC information in one plot.

    Figure size: 9cm width x 10cm height (extra space for legends)
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Figure size: 9cm width x 10cm height
    fig_width = 9 / 2.54
    fig_height = 10 / 2.54

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate layout
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Node colors by subcluster
    node_colors = [G.nodes[n].get('color', '#808080') for n in G.nodes()]

    # Node sizes proportional to AUC (normalize by AUC_VMAX)
    auc_values = np.array([G.nodes[n].get('auc', 0) for n in G.nodes()])
    auc_normalized = np.clip(auc_values / AUC_VMAX, 0, 1)
    node_sizes = min_size + (auc_normalized * (max_size - min_size))

    # Edge styling
    edge_widths = [0.3 for _ in G.edges()]
    edge_colors = ['#cccccc' for _ in G.edges()]

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.5,
                           edge_color=edge_colors)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.5)

    # Draw labels
    if show_labels:
        all_labels = {n: n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, ax=ax, labels=all_labels,
                                font_size=5, font_weight='normal')

    # === COLOR LEGEND (Subclusters) ===
    color_legend_elements = []
    for cluster in ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']:
        color = EOS_COLORS.get(cluster, '#808080')
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            color_legend_elements.append(plt.scatter([], [], c=color, s=40,
                                                      label=f'{cluster}'))

    if color_legend_elements:
        legend1 = ax.legend(handles=color_legend_elements, loc='upper left', fontsize=5,
                            title='Subcluster', title_fontsize=5, framealpha=0.9,
                            markerscale=0.8, labelspacing=0.3)
        ax.add_artist(legend1)

    # === SIZE LEGEND (AUC) ===
    # Create size legend showing AUC = 0, 0.45, 0.9 (based on Eos AUC scale)
    size_legend_elements = []
    for auc_val, label in [(0.0, 'AUC=0'), (0.45, 'AUC=0.45'), (0.9, 'AUC=0.9')]:
        size = min_size + ((auc_val / AUC_VMAX) * (max_size - min_size))
        size_legend_elements.append(plt.scatter([], [], c='gray', s=size,
                                                 label=label, alpha=0.7))

    legend2 = ax.legend(handles=size_legend_elements, loc='lower right', fontsize=5,
                        title='Node Size', title_fontsize=5, framealpha=0.9,
                        markerscale=1.0, labelspacing=1.0, handletextpad=1.5)

    ax.set_title(title, fontsize=7, fontweight='bold', pad=8)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved UNIFIED network: {output_path}")
    plt.close()


def create_subcluster_network(G, output_path, title="Eos Regulon Network (Subcluster)",
                               show_labels=True):
    """
    Create network visualization with nodes colored by Eos subcluster.

    Figure size: 8cm height x 9cm width
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Figure size: 9cm width x 8cm height
    fig_width = 9 / 2.54
    fig_height = 8 / 2.54

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Calculate layout - use absolute weights to avoid negative weight issues
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Node colors by subcluster
    node_colors = [G.nodes[n].get('color', '#808080') for n in G.nodes()]

    # Node sizes - fixed for cleaner look
    node_sizes = [80 for _ in G.nodes()]

    # Edge styling
    edge_widths = [0.3 for _ in G.edges()]
    edge_colors = ['#cccccc' for _ in G.edges()]

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.5,
                           edge_color=edge_colors)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.5)

    # Draw labels
    if show_labels:
        all_labels = {n: n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, ax=ax, labels=all_labels,
                                font_size=6, font_weight='normal')

    # Create legend for Eos subclusters
    legend_elements = []
    for cluster in ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']:
        color = EOS_COLORS.get(cluster, '#808080')
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            legend_elements.append(plt.scatter([], [], c=color, s=30,
                                               label=f'{cluster} ({count})'))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper left', fontsize=5,
                  title='Subcluster', title_fontsize=5, framealpha=0.8,
                  markerscale=0.8)

    ax.set_title(title, fontsize=6, pad=5)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_combined_figure(G, output_path, title_prefix="Eos Regulon Network", min_size=40, max_size=200):
    """
    Create combined figure with:
    - Top panel: AUC colored network
    - Bottom panel: Hybrid network (color by subcluster, size by AUC)

    Total figure size: 9cm width x 16cm height (8cm each)
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Figure size: 9cm width x 16cm height
    fig_width = 9 / 2.54
    fig_height = 16 / 2.54

    fig, axes = plt.subplots(2, 1, figsize=(fig_width, fig_height))

    # Calculate layout once for consistency - use absolute weights to avoid negative weight issues
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # ========== TOP PANEL: AUC colored ==========
    ax1 = axes[0]

    auc_values = np.array([G.nodes[n].get('auc', 0) for n in G.nodes()])
    norm = mcolors.Normalize(vmin=AUC_VMIN, vmax=AUC_VMAX)
    node_colors_auc = [AUC_CMAP(norm(a)) for a in auc_values]

    node_sizes_fixed = [60 for _ in G.nodes()]
    edge_widths = [0.3 for _ in G.edges()]
    edge_colors = ['#cccccc' for _ in G.edges()]

    nx.draw_networkx_edges(G, pos, ax=ax1, width=edge_widths, alpha=0.5, edge_color=edge_colors)
    nx.draw_networkx_nodes(G, pos, ax=ax1, node_color=node_colors_auc,
                           node_size=node_sizes_fixed, alpha=0.9, edgecolors='white', linewidths=0.3)

    all_labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, ax=ax1, labels=all_labels, font_size=5)

    sm = plt.cm.ScalarMappable(cmap=AUC_CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax1, shrink=0.5, aspect=10, pad=0.02)
    cbar.set_label('AUC', fontsize=5)
    cbar.ax.tick_params(labelsize=5)

    ax1.set_title(f'{title_prefix} (AUC)', fontsize=6, pad=3)
    ax1.axis('off')

    # ========== BOTTOM PANEL: Hybrid (color by cluster, size by AUC) ==========
    ax2 = axes[1]

    node_colors_cluster = [G.nodes[n].get('color', '#808080') for n in G.nodes()]
    # Node sizes proportional to AUC (normalize by AUC_VMAX)
    auc_normalized = np.clip(auc_values / AUC_VMAX, 0, 1)
    node_sizes_auc = min_size + (auc_normalized * (max_size - min_size))

    nx.draw_networkx_edges(G, pos, ax=ax2, width=edge_widths, alpha=0.5, edge_color=edge_colors)
    nx.draw_networkx_nodes(G, pos, ax=ax2, node_color=node_colors_cluster,
                           node_size=node_sizes_auc, alpha=0.9, edgecolors='white', linewidths=0.3)
    nx.draw_networkx_labels(G, pos, ax=ax2, labels=all_labels, font_size=5)

    # Legend for subclusters
    legend_elements = []
    for cluster in ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']:
        color = EOS_COLORS.get(cluster, '#808080')
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            legend_elements.append(plt.scatter([], [], c=color, s=20, label=f'{cluster} ({count})'))

    if legend_elements:
        ax2.legend(handles=legend_elements, loc='upper left', fontsize=4,
                   title='Subcluster', title_fontsize=4, framealpha=0.8, markerscale=0.6)

    # Add size legend annotation
    ax2.annotate('Node size ∝ AUC', xy=(0.98, 0.02), xycoords='axes fraction',
                fontsize=4, ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax2.set_title(f'{title_prefix} (Subcluster + AUC size)', fontsize=6, pad=3)
    ax2.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def save_network_data(G, corr_matrix, output_dir, prefix="eos_grn"):
    """Save network data for downstream analysis."""
    edges = []
    for u, v, data in G.edges(data=True):
        edges.append({
            'source': u,
            'target': v,
            'correlation': data.get('weight', 0),
            'abs_correlation': data.get('abs_weight', 0)
        })
    edge_df = pd.DataFrame(edges)
    edge_df.to_csv(output_dir / f'{prefix}_edges.csv', index=False)

    nodes = []
    for node, data in G.nodes(data=True):
        nodes.append({
            'regulon': node,
            'cluster': data.get('cluster', 'Unknown'),
            'auc': data.get('auc', 0),
            'uniqueness': data.get('uniqueness', 0),
            'zscore': data.get('zscore', 0),
            'degree': G.degree(node)
        })
    node_df = pd.DataFrame(nodes)
    node_df.to_csv(output_dir / f'{prefix}_nodes.csv', index=False)

    corr_matrix.to_csv(output_dir / f'{prefix}_correlation_matrix.csv')
    logging.info(f"Saved network data: {len(edges)} edges, {len(nodes)} nodes")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='GRN Network Visualization for Eos-Enriched Regulons')
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned')
    parser.add_argument('--regulon-file', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/302025_git/scripts/R19_20251126/results_eos_vs_mac_zscore/eos_subcluster_unique_relax.csv',
                        help='Path to Eos-enriched regulon CSV')
    parser.add_argument('--output-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/eos_grn_network')
    parser.add_argument('--cells-per-cluster', type=int, default=200)
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("EOS-ENRICHED REGULON GRN NETWORK (Updated)")
    logging.info("="*80)

    # Load regulon info
    logging.info(f"\nLoading regulon info from: {args.regulon_file}")
    regulon_info = pd.read_csv(args.regulon_file)
    target_regulons = regulon_info['tf'].tolist()
    logging.info(f"Found {len(target_regulons)} Eos-enriched regulons")

    # Extract cutoff level
    input_filename = Path(args.regulon_file).stem
    cutoff_level = "unknown"
    for level in ['relax', 'moderate', 'stringent']:
        if level in input_filename:
            cutoff_level = level
            break

    # Load AUC matrices
    logging.info("\n=== Loading AUC Matrices ===")
    eos_data, eos_clusters = load_auc_matrices(scenic_dir)

    # Calculate co-activity correlations
    logging.info("\n=== Calculating Co-activity Correlations ===")
    corr_matrix, combined_auc = calculate_regulon_coactivity(
        eos_data, eos_clusters, target_regulons, cells_per_cluster=args.cells_per_cluster
    )

    # Build networks at different thresholds
    for corr_thresh in [0.3, 0.6]:
        logging.info(f"\n=== Building Network (corr >= {corr_thresh}) ===")
        G = build_grn_network(corr_matrix, regulon_info, correlation_threshold=corr_thresh)

        thresh_str = str(corr_thresh).replace('.', '')

        # Create individual visualizations
        create_auc_network(
            G, output_dir / f'eos_grn_auc_corr{thresh_str}_{cutoff_level}.pdf',
            title=f'Eos GRN (n={G.number_of_nodes()}, corr≥{corr_thresh})'
        )

        create_subcluster_network(
            G, output_dir / f'eos_grn_subcluster_corr{thresh_str}_{cutoff_level}.pdf',
            title=f'Eos GRN (n={G.number_of_nodes()}, corr≥{corr_thresh})'
        )

        create_hybrid_network(
            G, output_dir / f'eos_grn_hybrid_corr{thresh_str}_{cutoff_level}.pdf',
            title=f'Eos GRN (n={G.number_of_nodes()}, corr≥{corr_thresh})'
        )

        # Create UNIFIED figure (single figure with color=subcluster, size=AUC)
        create_unified_network(
            G, output_dir / f'eos_grn_unified_corr{thresh_str}_{cutoff_level}.pdf',
            title=f'Eos GRN: Subcluster + AUC (n={G.number_of_nodes()}, corr≥{corr_thresh})'
        )

        # Create combined figure
        create_combined_figure(
            G, output_dir / f'eos_grn_combined_corr{thresh_str}_{cutoff_level}.pdf',
            title_prefix=f'Eos GRN (n={G.number_of_nodes()}, corr≥{corr_thresh})'
        )

        # Save network data
        save_network_data(G, corr_matrix, output_dir, prefix=f"eos_grn_corr{thresh_str}_{cutoff_level}")

        # Print summary
        print(f"\n--- Correlation >= {corr_thresh} ---")
        print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
        if G.number_of_nodes() > 0:
            print(f"Density: {nx.density(G):.4f}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)
    print(f"\nOutput directory: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
