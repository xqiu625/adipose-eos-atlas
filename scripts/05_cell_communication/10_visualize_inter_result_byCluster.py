import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']

def process_cluster_results(cluster):
    # Load significant means
    sig_means = pd.read_csv(f"./cellphonedb_results_{cluster}_vs_Mac/statistical_analysis_significant_means.txt", sep="\t")
    
    # Identify cluster-macrophage interactions
    cluster_mac_columns = [f"{cluster}|Macrophage", f"Macrophage|{cluster}"]
    
    # Filter and process data
    cluster_mac_interactions = sig_means[sig_means[cluster_mac_columns].notna().any(axis=1)].copy()
    cluster_mac_interactions['max_value'] = cluster_mac_interactions[cluster_mac_columns].max(axis=1)
    cluster_mac_interactions = cluster_mac_interactions.sort_values('max_value', ascending=False)
    
    return cluster_mac_interactions.head(20)  # Return top 20 interactions

# Process results for each cluster
cluster_results = {cluster: process_cluster_results(cluster) for cluster in clusters}

# Create a comparison heatmap
plt.figure(figsize=(80, 15))
for i, (cluster, results) in enumerate(cluster_results.items()):
    cluster_mac_columns = [f"{cluster}|Macrophage", f"Macrophage|{cluster}"]
    
    plt.subplot(1, len(cluster_results), i+1)
    sns.heatmap(results[cluster_mac_columns].T, annot=True, fmt='.2f', cmap='YlOrRd')
    plt.title(f"Top 20 Interactions - {cluster}")
    plt.xlabel("Interacting Pairs")
    plt.ylabel("Cell Type Interactions")
    plt.xticks(range(len(results)), results['interacting_pair'], rotation=90, ha='right')

plt.tight_layout()
plt.savefig('cluster_comparison_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# Print top interactions for each cluster
for cluster, results in cluster_results.items():
    print(f"\nTop 5 interactions for {cluster}:")
    print(results[['interacting_pair', 'partner_a', 'partner_b', 'max_value']].head())

# Create separate heatmaps for each cluster
for cluster, results in cluster_results.items():
    cluster_mac_columns = [f"{cluster}|Macrophage", f"Macrophage|{cluster}"]
    
    plt.figure(figsize=(20, 15))  # Adjust figure size as needed
    sns.heatmap(results[cluster_mac_columns].T, annot=True, fmt='.2f', cmap='YlOrRd')
    plt.title(f"Top 20 Interactions - {cluster}")
    plt.xlabel("Interacting Pairs")
    plt.ylabel("Cell Type Interactions")
    plt.xticks(range(len(results)), results['interacting_pair'], rotation=90, ha='right')
    
    plt.tight_layout()
    plt.savefig(f'heatmap_{cluster}.png', dpi=300, bbox_inches='tight')
    plt.close()

print("Heatmaps saved for each cluster.")
