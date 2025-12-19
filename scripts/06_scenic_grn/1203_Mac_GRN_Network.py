#!/usr/bin/env python3
"""
GRN Network Visualization for Mac-Enriched Regulons

Creates a gene regulatory network visualization similar to Fig 3B in Frontiers paper.
- Nodes: Transcription factor regulons
- Edges: Co-regulation/co-activity relationships
- Node color: Mac subtype enrichment
- Node size: Regulon activity (mean AUC)

Input: Mac subtype-unique regulons from SCENIC analysis
Output: Network visualization (PDF) and edge list (CSV)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

# Mac subtype colors
MAC_COLORS = {
    'M1': '#E41A1C',           # Red
    'M2_Functional': '#377EB8', # Blue
    'M2_Surface': '#4DAF4A',    # Green
}

# Subsampling parameters
CELLS_PER_CLUSTER = 200
RANDOM_SEED = 42


def load_auc_matrices(scenic_dir):
    """Load AUC matrices for Mac subtypes"""

    mac_types = ['M1', 'M2_Functional', 'M2_Surface']
    mac_data = {}

    for mac_type in mac_types:
        auc_file = scenic_dir / "macrophages" / mac_type / "step3_auc_matrix.csv"
        if auc_file.exists():
            mac_data[mac_type] = pd.read_csv(auc_file, index_col=0)
            logging.info(f"Loaded {mac_type}: {mac_data[mac_type].shape[0]} cells, {mac_data[mac_type].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    loaded_mac = [m for m in mac_types if m in mac_data]
    return mac_data, loaded_mac


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_regulon_coactivity(mac_data, mac_types, target_regulons, cells_per_cluster=200):
    """
    Calculate co-activity correlation between regulons across all Mac cells.
    Returns correlation matrix for network edges.
    """
    logging.info("Calculating regulon co-activity correlations...")
    np.random.seed(RANDOM_SEED)

    # Build TF mapping
    tf_to_regulon = {}
    for mac_type in mac_types:
        if mac_type not in mac_data:
            continue
        for col in mac_data[mac_type].columns:
            tf_name = extract_tf_name(col)
            if tf_name not in tf_to_regulon:
                tf_to_regulon[tf_name] = {}
            tf_to_regulon[tf_name][mac_type] = col

    # Collect AUC values for target regulons across all cells
    all_cells_auc = []

    for mac_type in mac_types:
        if mac_type not in mac_data:
            continue
        mac_df = mac_data[mac_type]
        n_cells = mac_df.shape[0]

        # Subsample
        if n_cells >= cells_per_cluster:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=False)
        else:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=True)

        sampled_df = mac_df.iloc[sampled_idx]

        # Extract AUC for target regulons
        cell_data = {}
        for tf in target_regulons:
            if tf in tf_to_regulon and mac_type in tf_to_regulon[tf]:
                regulon_col = tf_to_regulon[tf][mac_type]
                cell_data[tf] = sampled_df[regulon_col].values
            else:
                cell_data[tf] = np.zeros(cells_per_cluster)

        mac_auc_df = pd.DataFrame(cell_data)
        mac_auc_df['subtype'] = mac_type
        all_cells_auc.append(mac_auc_df)

    # Combine all cells
    combined_auc = pd.concat(all_cells_auc, ignore_index=True)
    auc_matrix = combined_auc[target_regulons]

    # Calculate correlation matrix
    corr_matrix = auc_matrix.corr(method='pearson')

    logging.info(f"Calculated {len(target_regulons)}x{len(target_regulons)} correlation matrix")

    return corr_matrix, combined_auc


def build_grn_network(corr_matrix, regulon_info, correlation_threshold=0.3):
    """
    Build NetworkX graph from correlation matrix.

    Args:
        corr_matrix: Regulon-regulon correlation matrix
        regulon_info: DataFrame with regulon metadata (max_mac_subtype, max_mac_auc, etc.)
        correlation_threshold: Minimum correlation to create edge

    Returns:
        NetworkX graph object
    """
    logging.info(f"Building network with correlation threshold >= {correlation_threshold}")

    G = nx.Graph()

    # Add nodes with attributes
    for tf in corr_matrix.columns:
        if tf in regulon_info['tf'].values:
            info = regulon_info[regulon_info['tf'] == tf].iloc[0]
            G.add_node(tf,
                       subtype=info['max_mac_subtype'],
                       auc=info['max_mac_auc'],
                       uniqueness=info['uniqueness_score'],
                       color=MAC_COLORS.get(info['max_mac_subtype'], '#808080'))
        else:
            G.add_node(tf, subtype='Unknown', auc=0, uniqueness=0, color='#808080')

    # Add edges based on correlation
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


def detect_modules(G, resolution=1.0):
    """
    Detect network modules using Louvain community detection.
    """
    try:
        import community as community_louvain
        partition = community_louvain.best_partition(G, resolution=resolution, random_state=RANDOM_SEED)

        # Add module info to nodes
        for node, module in partition.items():
            G.nodes[node]['module'] = module

        n_modules = len(set(partition.values()))
        logging.info(f"Detected {n_modules} modules")

        return partition
    except ImportError:
        logging.warning("python-louvain not installed. Skipping module detection.")
        return None


def create_grn_visualization(G, output_path, title="Mac-Enriched Regulon Network",
                              node_size_scale=500, edge_width_scale=2):
    """
    Create network visualization with:
    - Node color by Mac subtype
    - Node size by AUC activity
    - Edge width by correlation strength
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Calculate layout
    pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Prepare node attributes
    node_colors = [G.nodes[n].get('color', '#808080') for n in G.nodes()]
    node_sizes = [G.nodes[n].get('auc', 0.1) * node_size_scale + 100 for n in G.nodes()]
    node_labels = {n: n for n in G.nodes()}

    # Prepare edge attributes
    edge_weights = [G.edges[e].get('abs_weight', 0.5) * edge_width_scale for e in G.edges()]
    edge_colors = ['#2166ac' if G.edges[e].get('weight', 0) > 0 else '#b2182b' for e in G.edges()]

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 12))

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_weights, alpha=0.5,
                           edge_color=edge_colors)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9, edgecolors='black', linewidths=1)

    # Draw labels
    nx.draw_networkx_labels(G, pos, ax=ax, labels=node_labels, font_size=8, font_weight='bold')

    # Create legend for Mac subtypes
    legend_elements = []
    for subtype, color in MAC_COLORS.items():
        count = sum(1 for n in G.nodes() if G.nodes[n].get('subtype') == subtype)
        if count > 0:
            legend_elements.append(plt.scatter([], [], c=color, s=100, label=f'{subtype} (n={count})'))

    ax.legend(handles=legend_elements, loc='upper left', title='Mac Subtype', fontsize=10)

    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_module_visualization(G, partition, output_path, title="Mac Regulon Modules"):
    """
    Create network visualization colored by detected modules.
    """
    if G.number_of_nodes() == 0 or partition is None:
        logging.warning("No nodes or modules. Skipping module visualization.")
        return

    # Generate module colors
    n_modules = len(set(partition.values()))
    module_colors = plt.cm.Set3(np.linspace(0, 1, max(n_modules, 1)))

    pos = nx.spring_layout(G, k=2, iterations=100, seed=RANDOM_SEED)

    # Color nodes by module
    node_colors = [module_colors[partition.get(n, 0)] for n in G.nodes()]
    node_sizes = [G.nodes[n].get('auc', 0.1) * 500 + 100 for n in G.nodes()]

    edge_weights = [G.edges[e].get('abs_weight', 0.5) * 2 for e in G.edges()]

    fig, ax = plt.subplots(figsize=(14, 12))

    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_weights, alpha=0.3, edge_color='gray')
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9, edgecolors='black', linewidths=1)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight='bold')

    # Legend for modules
    legend_elements = []
    for module_id in sorted(set(partition.values())):
        nodes_in_module = [n for n, m in partition.items() if m == module_id]
        # Get dominant subtype in this module
        subtypes = [G.nodes[n].get('subtype', 'Unknown') for n in nodes_in_module]
        dominant_subtype = max(set(subtypes), key=subtypes.count) if subtypes else 'Unknown'
        legend_elements.append(plt.scatter([], [], c=[module_colors[module_id]], s=100,
                                           label=f'Module {module_id} (n={len(nodes_in_module)}, {dominant_subtype})'))

    ax.legend(handles=legend_elements, loc='upper left', title='Modules', fontsize=9)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def analyze_module_enrichment(G, partition, mac_types):
    """
    Analyze which Mac subtypes are enriched in each module.
    """
    if partition is None:
        return None

    results = []

    for module_id in sorted(set(partition.values())):
        nodes_in_module = [n for n, m in partition.items() if m == module_id]

        subtype_counts = {s: 0 for s in mac_types}
        for node in nodes_in_module:
            subtype = G.nodes[node].get('subtype', 'Unknown')
            if subtype in subtype_counts:
                subtype_counts[subtype] += 1

        total = len(nodes_in_module)
        dominant_subtype = max(subtype_counts, key=subtype_counts.get)

        results.append({
            'module': module_id,
            'n_regulons': total,
            'dominant_subtype': dominant_subtype,
            'dominant_count': subtype_counts[dominant_subtype],
            'dominant_pct': subtype_counts[dominant_subtype] / total * 100 if total > 0 else 0,
            **{f'{s}_count': subtype_counts[s] for s in mac_types},
            'regulons': ', '.join(sorted(nodes_in_module))
        })

    return pd.DataFrame(results)


def save_network_data(G, corr_matrix, output_dir, prefix="mac_grn"):
    """
    Save network data for downstream analysis.
    """
    # Save edge list
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

    # Save node attributes
    nodes = []
    for node, data in G.nodes(data=True):
        nodes.append({
            'regulon': node,
            'subtype': data.get('subtype', 'Unknown'),
            'auc': data.get('auc', 0),
            'uniqueness': data.get('uniqueness', 0),
            'module': data.get('module', -1),
            'degree': G.degree(node)
        })
    node_df = pd.DataFrame(nodes)
    node_df.to_csv(output_dir / f'{prefix}_nodes.csv', index=False)

    # Save correlation matrix
    corr_matrix.to_csv(output_dir / f'{prefix}_correlation_matrix.csv')

    logging.info(f"Saved network data: {len(edges)} edges, {len(nodes)} nodes")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='GRN Network Visualization for Mac-Enriched Regulons')
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--regulon-file', type=str, required=True,
                        help='Path to Mac-enriched regulon CSV (from 1202 script)')
    parser.add_argument('--output-dir', type=str, default='results_mac_grn_network',
                        help='Path to output directory')
    parser.add_argument('--correlation-threshold', type=float, default=0.3,
                        help='Minimum correlation for network edges')
    parser.add_argument('--cells-per-cluster', type=int, default=200,
                        help='Number of cells to subsample per cluster')
    args = parser.parse_args()

    scenic_dir = Path(args.scenic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("="*80)
    logging.info("MAC-ENRICHED REGULON GRN NETWORK")
    logging.info("="*80)

    # Load regulon info
    logging.info(f"\nLoading regulon info from: {args.regulon_file}")
    regulon_info = pd.read_csv(args.regulon_file)
    target_regulons = regulon_info['tf'].tolist()
    logging.info(f"Found {len(target_regulons)} Mac-enriched regulons")

    # Extract cutoff level from input filename (e.g., "mac_subtype_unique_relax.csv" -> "relax")
    input_filename = Path(args.regulon_file).stem
    cutoff_level = "unknown"
    for level in ['relax', 'moderate', 'stringent']:
        if level in input_filename:
            cutoff_level = level
            break
    logging.info(f"Detected cutoff level: {cutoff_level}")

    # Load AUC matrices
    logging.info("\n=== Loading AUC Matrices ===")
    mac_data, mac_types = load_auc_matrices(scenic_dir)

    # Calculate co-activity correlations
    logging.info("\n=== Calculating Co-activity Correlations ===")
    corr_matrix, combined_auc = calculate_regulon_coactivity(
        mac_data, mac_types, target_regulons, cells_per_cluster=args.cells_per_cluster
    )

    # Build network
    logging.info("\n=== Building GRN Network ===")
    G = build_grn_network(corr_matrix, regulon_info, correlation_threshold=args.correlation_threshold)

    # Detect modules
    logging.info("\n=== Detecting Network Modules ===")
    partition = detect_modules(G)

    # Create visualizations
    logging.info("\n=== Creating Visualizations ===")

    # Main network colored by subtype
    create_grn_visualization(
        G, output_dir / f'mac_grn_network_by_subtype_{cutoff_level}.pdf',
        title=f'Mac-Enriched Regulon Network ({cutoff_level}, n={G.number_of_nodes()}, corr≥{args.correlation_threshold})'
    )

    # Network colored by module
    if partition:
        create_module_visualization(
            G, partition, output_dir / f'mac_grn_network_by_module_{cutoff_level}.pdf',
            title=f'Mac Regulon Modules ({cutoff_level}, n={len(set(partition.values()))} modules)'
        )

        # Module enrichment analysis
        module_enrichment = analyze_module_enrichment(G, partition, mac_types)
        if module_enrichment is not None:
            module_enrichment.to_csv(output_dir / f'mac_module_enrichment_{cutoff_level}.csv', index=False)
            logging.info("\nModule Enrichment Summary:")
            print(module_enrichment[['module', 'n_regulons', 'dominant_subtype', 'dominant_pct']].to_string())

    # Save network data
    logging.info("\n=== Saving Network Data ===")
    save_network_data(G, corr_matrix, output_dir, prefix=f"mac_grn_{cutoff_level}")

    # Summary statistics
    logging.info("\n" + "="*80)
    logging.info("NETWORK SUMMARY")
    logging.info("="*80)
    print(f"Nodes (regulons): {G.number_of_nodes()}")
    print(f"Edges (correlations ≥ {args.correlation_threshold}): {G.number_of_edges()}")
    print(f"Network density: {nx.density(G):.4f}")
    if G.number_of_nodes() > 0:
        print(f"Average clustering coefficient: {nx.average_clustering(G):.4f}")
        degrees = [d for n, d in G.degree()]
        print(f"Average degree: {np.mean(degrees):.2f}")

        # Hub regulons (highest degree)
        hub_nodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]
        print("\nTop 10 Hub Regulons:")
        for node, degree in hub_nodes:
            subtype = G.nodes[node].get('subtype', 'Unknown')
            print(f"  {node}: degree={degree}, subtype={subtype}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)

    print(f"\nOutput directory: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
