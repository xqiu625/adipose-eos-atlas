#!/usr/bin/env python3
"""
Quick GRN visualization from pre-computed CSV files.
Uses updated AUC scales:
- Eos: 0-0.9
- Mac: 0-0.7
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from pathlib import Path
import argparse

# Color schemes
EOS_COLORS = {
    'Eos_C0': '#6ca6cd', 'Eos_C1': '#e84a60', 'Eos_C2': '#f37026',
    'Eos_C3': '#ccabcb', 'Eos_C4': '#b13b49', 'Eos_C5': '#5ba4a4', 'Eos_C6': '#2f3d7f',
}
MAC_COLORS = {
    'M1': '#cc426c', 'M2_Functional': '#8fc43e', 'M2_Surface': '#58429a',
}

AUC_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'auc', ['#f7f7f7', '#fee0d2', '#fc9272', '#de2d26']
)


def create_unified_network(nodes_df, edges_df, output_path, cell_type='eos',
                           title="GRN Network", min_size=20, max_size=300):
    """Create unified network with color by cluster/subtype and size by AUC."""

    # Set AUC scale based on cell type
    if cell_type == 'eos':
        auc_vmax = 0.9
        colors = EOS_COLORS
        cluster_col = 'cluster'
        cluster_order = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']
    else:
        auc_vmax = 0.7
        colors = MAC_COLORS
        cluster_col = 'subtype'
        cluster_order = ['M1', 'M2_Functional', 'M2_Surface']

    # Build network
    G = nx.Graph()
    for _, row in nodes_df.iterrows():
        cluster = row[cluster_col]
        G.add_node(row['regulon'],
                   cluster=cluster,
                   auc=row['auc'],
                   color=colors.get(cluster, '#808080'))

    for _, row in edges_df.iterrows():
        if row['source'] in G.nodes() and row['target'] in G.nodes():
            G.add_edge(row['source'], row['target'],
                       weight=row['correlation'],
                       abs_weight=abs(row['correlation']))

    if G.number_of_nodes() == 0:
        print("No nodes in network!")
        return

    # Figure setup
    fig_width = 9 / 2.54
    fig_height = 10 / 2.54
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Layout
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G, weight='abs_weight')
    else:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=42)

    # Node colors and sizes
    node_colors = [G.nodes[n].get('color', '#808080') for n in G.nodes()]
    auc_values = np.array([G.nodes[n].get('auc', 0) for n in G.nodes()])
    auc_normalized = np.clip(auc_values / auc_vmax, 0, 1)
    node_sizes = min_size + (auc_normalized * (max_size - min_size))

    # Draw
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.3, alpha=0.5, edge_color='#cccccc')
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           edgecolors='white', linewidths=0.5)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=5)

    # Color legend
    color_legend = []
    for cluster in cluster_order:
        color = colors.get(cluster, '#808080')
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            label = cluster.replace('_', ' ') if cell_type == 'mac' else cluster
            color_legend.append(plt.scatter([], [], c=color, s=40, label=label))

    if color_legend:
        legend1 = ax.legend(handles=color_legend, loc='upper left', fontsize=5,
                            title='Cluster' if cell_type == 'eos' else 'Subtype',
                            title_fontsize=5, framealpha=0.9)
        ax.add_artist(legend1)

    # Size legend
    size_legend = []
    mid_val = auc_vmax / 2
    for auc_val, label in [(0.0, 'AUC=0'), (mid_val, f'AUC={mid_val}'), (auc_vmax, f'AUC={auc_vmax}')]:
        size = min_size + ((auc_val / auc_vmax) * (max_size - min_size))
        size_legend.append(plt.scatter([], [], c='gray', s=size, label=label, alpha=0.7))

    ax.legend(handles=size_legend, loc='lower right', fontsize=5,
              title='Node Size', title_fontsize=5, framealpha=0.9,
              labelspacing=1.0, handletextpad=1.5)

    ax.set_title(title, fontsize=7, fontweight='bold', pad=8)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Quick GRN visualization from CSV')
    parser.add_argument('--grn-dir', type=str, default='/Users/xinruqiu/Desktop/302025/grn',
                        help='Directory containing GRN CSV files')
    parser.add_argument('--output-dir', type=str, default='/Users/xinruqiu/Desktop/302025/grn/figures',
                        help='Output directory for figures')
    args = parser.parse_args()

    grn_dir = Path(args.grn_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process Eos networks
    for corr in ['03', '06']:
        nodes_file = grn_dir / f'eos_grn_corr{corr}_relax_nodes.csv'
        edges_file = grn_dir / f'eos_grn_corr{corr}_relax_edges.csv'

        if nodes_file.exists() and edges_file.exists():
            nodes_df = pd.read_csv(nodes_file)
            edges_df = pd.read_csv(edges_file)

            corr_val = 0.3 if corr == '03' else 0.6
            create_unified_network(
                nodes_df, edges_df,
                output_dir / f'eos_grn_unified_corr{corr}.pdf',
                cell_type='eos',
                title=f'Eos GRN (n={len(nodes_df)}, corr≥{corr_val})'
            )

    # Process Mac networks
    for corr in ['03', '06']:
        nodes_file = grn_dir / f'mac_grn_corr{corr}_relax_nodes.csv'
        edges_file = grn_dir / f'mac_grn_corr{corr}_relax_edges.csv'

        if nodes_file.exists() and edges_file.exists():
            nodes_df = pd.read_csv(nodes_file)
            edges_df = pd.read_csv(edges_file)

            corr_val = 0.3 if corr == '03' else 0.6
            create_unified_network(
                nodes_df, edges_df,
                output_dir / f'mac_grn_unified_corr{corr}.pdf',
                cell_type='mac',
                title=f'Mac GRN (n={len(nodes_df)}, corr≥{corr_val})'
            )

    print(f"\nOutput: {output_dir}")


if __name__ == "__main__":
    main()
