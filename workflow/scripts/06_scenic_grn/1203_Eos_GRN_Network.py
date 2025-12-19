#!/usr/bin/env python3
"""
GRN Network Visualization for Eos-Enriched Regulons

Creates a gene regulatory network visualization similar to Fig 3B in Frontiers paper.
- Nodes: Transcription factor regulons
- Edges: Co-regulation/co-activity relationships
- Node color: Eos subcluster enrichment
- Node size: Regulon activity (mean AUC)

Input: Eos subcluster-unique regulons from SCENIC analysis
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

# Eos subcluster colors
EOS_COLORS = {
    'Eos_C0': '#E41A1C',  # Red
    'Eos_C1': '#377EB8',  # Blue
    'Eos_C2': '#4DAF4A',  # Green
    'Eos_C3': '#984EA3',  # Purple
    'Eos_C4': '#FF7F00',  # Orange
    'Eos_C5': '#FFFF33',  # Yellow
    'Eos_C6': '#A65628',  # Brown
}

# Subsampling parameters
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
            logging.info(f"Loaded {cluster}: {eos_data[cluster].shape[0]} cells, {eos_data[cluster].shape[1]} regulons")
        else:
            logging.warning(f"File not found: {auc_file}")

    loaded_eos = [c for c in eos_clusters if c in eos_data]
    return eos_data, loaded_eos


def extract_tf_name(regulon_name):
    """Extract TF name from regulon format"""
    if '(' in regulon_name:
        return regulon_name.split('(')[0]
    return regulon_name


def calculate_regulon_coactivity(eos_data, eos_clusters, target_regulons, cells_per_cluster=200):
    """
    Calculate co-activity correlation between regulons across all Eos cells.
    Returns correlation matrix for network edges.
    """
    logging.info("Calculating regulon co-activity correlations...")
    np.random.seed(RANDOM_SEED)

    # Build TF mapping
    tf_to_regulon = {}
    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        for col in eos_data[cluster].columns:
            tf_name = extract_tf_name(col)
            if tf_name not in tf_to_regulon:
                tf_to_regulon[tf_name] = {}
            tf_to_regulon[tf_name][cluster] = col

    # Collect AUC values for target regulons across all cells
    all_cells_auc = []

    for cluster in eos_clusters:
        if cluster not in eos_data:
            continue
        cluster_df = eos_data[cluster]
        n_cells = cluster_df.shape[0]

        # Subsample
        if n_cells >= cells_per_cluster:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=False)
        else:
            sampled_idx = np.random.choice(n_cells, cells_per_cluster, replace=True)

        sampled_df = cluster_df.iloc[sampled_idx]

        # Extract AUC for target regulons
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
        regulon_info: DataFrame with regulon metadata (max_eos_cluster, max_eos_auc, etc.)
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
                       cluster=info['max_eos_cluster'],
                       auc=info['max_eos_auc'],
                       uniqueness=info['uniqueness_score'],
                       color=EOS_COLORS.get(info['max_eos_cluster'], '#808080'))
        else:
            G.add_node(tf, cluster='Unknown', auc=0, uniqueness=0, color='#808080')

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


def create_grn_visualization(G, output_path, title="Eos-Enriched Regulon Network",
                              node_size_scale=500, edge_width_scale=2):
    """
    Create network visualization with:
    - Node color by Eos subcluster
    - Node size by AUC activity
    - Edge width by correlation strength
    """
    if G.number_of_nodes() == 0:
        logging.warning("No nodes in network. Skipping visualization.")
        return

    # Calculate layout
    # Use spring layout with seed for reproducibility
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

    # Create legend for Eos subclusters
    legend_elements = []
    for cluster, color in EOS_COLORS.items():
        # Count nodes in this cluster
        count = sum(1 for n in G.nodes() if G.nodes[n].get('cluster') == cluster)
        if count > 0:
            legend_elements.append(plt.scatter([], [], c=color, s=100, label=f'{cluster} (n={count})'))

    ax.legend(handles=legend_elements, loc='upper left', title='Eos Subcluster', fontsize=10)

    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(output_path).replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def create_module_visualization(G, partition, output_path, title="Eos Regulon Modules"):
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
        # Get dominant cluster in this module
        clusters = [G.nodes[n].get('cluster', 'Unknown') for n in nodes_in_module]
        dominant_cluster = max(set(clusters), key=clusters.count) if clusters else 'Unknown'
        legend_elements.append(plt.scatter([], [], c=[module_colors[module_id]], s=100,
                                           label=f'Module {module_id} (n={len(nodes_in_module)}, {dominant_cluster})'))

    ax.legend(handles=legend_elements, loc='upper left', title='Modules', fontsize=9)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logging.info(f"Saved: {output_path}")
    plt.close()


def analyze_module_enrichment(G, partition, eos_clusters):
    """
    Analyze which Eos subclusters are enriched in each module.
    """
    if partition is None:
        return None

    results = []

    for module_id in sorted(set(partition.values())):
        nodes_in_module = [n for n, m in partition.items() if m == module_id]

        cluster_counts = {c: 0 for c in eos_clusters}
        for node in nodes_in_module:
            cluster = G.nodes[node].get('cluster', 'Unknown')
            if cluster in cluster_counts:
                cluster_counts[cluster] += 1

        total = len(nodes_in_module)
        dominant_cluster = max(cluster_counts, key=cluster_counts.get)

        results.append({
            'module': module_id,
            'n_regulons': total,
            'dominant_cluster': dominant_cluster,
            'dominant_count': cluster_counts[dominant_cluster],
            'dominant_pct': cluster_counts[dominant_cluster] / total * 100 if total > 0 else 0,
            **{f'{c}_count': cluster_counts[c] for c in eos_clusters},
            'regulons': ', '.join(sorted(nodes_in_module))
        })

    return pd.DataFrame(results)


def save_network_data(G, corr_matrix, output_dir, prefix="eos_grn"):
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
            'cluster': data.get('cluster', 'Unknown'),
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

    parser = argparse.ArgumentParser(description='GRN Network Visualization for Eos-Enriched Regulons')
    parser.add_argument('--scenic-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302025/scenic_results/scenic_cleaned',
                        help='Path to SCENIC results directory')
    parser.add_argument('--regulon-file', type=str, required=True,
                        help='Path to Eos-enriched regulon CSV (from 1202 script)')
    parser.add_argument('--output-dir', type=str, default='results_eos_grn_network',
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
    logging.info("EOS-ENRICHED REGULON GRN NETWORK")
    logging.info("="*80)

    # Load regulon info
    logging.info(f"\nLoading regulon info from: {args.regulon_file}")
    regulon_info = pd.read_csv(args.regulon_file)
    target_regulons = regulon_info['tf'].tolist()
    logging.info(f"Found {len(target_regulons)} Eos-enriched regulons")

    # Extract cutoff level from input filename (e.g., "eos_subcluster_unique_relax.csv" -> "relax")
    input_filename = Path(args.regulon_file).stem
    cutoff_level = "unknown"
    for level in ['relax', 'moderate', 'stringent']:
        if level in input_filename:
            cutoff_level = level
            break
    logging.info(f"Detected cutoff level: {cutoff_level}")

    # Load AUC matrices
    logging.info("\n=== Loading AUC Matrices ===")
    eos_data, eos_clusters = load_auc_matrices(scenic_dir)

    # Calculate co-activity correlations
    logging.info("\n=== Calculating Co-activity Correlations ===")
    corr_matrix, combined_auc = calculate_regulon_coactivity(
        eos_data, eos_clusters, target_regulons, cells_per_cluster=args.cells_per_cluster
    )

    # Build network
    logging.info("\n=== Building GRN Network ===")
    G = build_grn_network(corr_matrix, regulon_info, correlation_threshold=args.correlation_threshold)

    # Detect modules
    logging.info("\n=== Detecting Network Modules ===")
    partition = detect_modules(G)

    # Create visualizations
    logging.info("\n=== Creating Visualizations ===")

    # Main network colored by subcluster
    create_grn_visualization(
        G, output_dir / f'eos_grn_network_by_cluster_{cutoff_level}.pdf',
        title=f'Eos-Enriched Regulon Network ({cutoff_level}, n={G.number_of_nodes()}, corr≥{args.correlation_threshold})'
    )

    # Network colored by module
    if partition:
        create_module_visualization(
            G, partition, output_dir / f'eos_grn_network_by_module_{cutoff_level}.pdf',
            title=f'Eos Regulon Modules ({cutoff_level}, n={len(set(partition.values()))} modules)'
        )

        # Module enrichment analysis
        module_enrichment = analyze_module_enrichment(G, partition, eos_clusters)
        if module_enrichment is not None:
            module_enrichment.to_csv(output_dir / f'eos_module_enrichment_{cutoff_level}.csv', index=False)
            logging.info("\nModule Enrichment Summary:")
            print(module_enrichment[['module', 'n_regulons', 'dominant_cluster', 'dominant_pct']].to_string())

    # Save network data
    logging.info("\n=== Saving Network Data ===")
    save_network_data(G, corr_matrix, output_dir, prefix=f"eos_grn_{cutoff_level}")

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
            cluster = G.nodes[node].get('cluster', 'Unknown')
            print(f"  {node}: degree={degree}, cluster={cluster}")

    logging.info("\n" + "="*80)
    logging.info("ANALYSIS COMPLETE")
    logging.info("="*80)

    print(f"\nOutput directory: {output_dir}")
    print(f'scp -r "xinruqiu@cluster.hpcc.ucr.edu:{output_dir.absolute()}" .')


if __name__ == "__main__":
    main()
