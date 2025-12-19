import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

clusters = ['Eos_C0', 'Eos_C1', 'Eos_C2', 'Eos_C3', 'Eos_C4', 'Eos_C5', 'Eos_C6']

def analyze_cluster_interactions(cluster):
    # Load significant means
    sig_means = pd.read_csv(f"./cellphonedb_results_{cluster}_vs_Mac/statistical_analysis_significant_means.txt", sep="\t")
    
    # Identify cluster-macrophage interactions
    cluster_mac_columns = [f"{cluster}|Macrophage", f"Macrophage|{cluster}"]
    
    # Filter for significant interactions
    significant_interactions = sig_means[sig_means[cluster_mac_columns].notna().any(axis=1)]
    
    # Calculate metrics
    num_interactions = len(significant_interactions)
    total_strength = significant_interactions[cluster_mac_columns].sum().sum()
    max_strength = significant_interactions[cluster_mac_columns].max().max()
    
    return {
        'cluster': cluster,
        'num_interactions': num_interactions,
        'total_strength': total_strength,
        'max_strength': max_strength,
        'avg_strength': total_strength / (num_interactions * 2) if num_interactions > 0 else 0
    }

# Analyze results for each cluster
cluster_analyses = [analyze_cluster_interactions(cluster) for cluster in clusters]
analysis_df = pd.DataFrame(cluster_analyses)

# Print summary
print(analysis_df)


# 1. Scatter plot of number of interactions vs average strength
plt.figure(figsize=(10, 8))
sns.scatterplot(x='num_interactions', y='avg_strength', hue='cluster', data=analysis_df, s=100)
plt.title('Number of Interactions vs Average Strength')
plt.xlabel('Number of Interactions')
plt.ylabel('Average Interaction Strength')
for i, row in analysis_df.iterrows():
    plt.annotate(row['cluster'], (row['num_interactions'], row['avg_strength']))

plt.tight_layout()
plt.savefig('interactions_vs_strength_scatter_by_cluster.png')
plt.close()

# 2. Heatmap of top interactions across clusters
top_interactions = pd.DataFrame({cluster: results['interacting_pair'].head(20).tolist() 
                                 for cluster, results in cluster_results.items()})
all_top_interactions = pd.Series(top_interactions.values.ravel()).unique()

heatmap_data = pd.DataFrame(index=all_top_interactions, columns=clusters)
for cluster, results in cluster_results.items():
    heatmap_data[cluster] = results.set_index('interacting_pair')['max_value']

heatmap_data = heatmap_data.fillna(0)

plt.figure(figsize=(20, 30))
sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd')
plt.title('Top Eosinophil-Macrophage Interactions Across clusters')
plt.xlabel('cluster')
plt.ylabel('Interacting Pair')
plt.tight_layout()
plt.savefig('top_interactions_heatmap_by_cluster.png')
plt.close()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

conditions = ['WT_Male_HFD', 'WT_Female_HFD', 'KO_Male_HFD', 'KO_Female_HFD', 'WT_Female_HFDInfected', 'KO_Female_HFDInfected']

def process_condition_results(condition):
    # Load significant means
    sig_means = pd.read_csv(f"./cellphonedb_results_{condition}/statistical_analysis_significant_means.txt", sep="\t")
    
    # Identify eosinophil-macrophage interactions
    eos_mac_columns = [col for col in sig_means.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]
    
    # Filter and process data
    eos_mac_interactions = sig_means[sig_means[eos_mac_columns].notna().any(axis=1)].copy()
    eos_mac_interactions['max_value'] = eos_mac_interactions[eos_mac_columns].max(axis=1)
    eos_mac_interactions = eos_mac_interactions.sort_values('max_value', ascending=False)
    
    return eos_mac_interactions

# Process results for each condition
condition_results = {condition: process_condition_results(condition) for condition in conditions}

# Analyze interaction strengths
def analyze_interaction_strengths(condition_results):
    analysis = {}
    for condition, results in condition_results.items():
        eos_mac_columns = [col for col in results.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]
        analysis[condition] = {
            'num_interactions': len(results),
            'avg_strength': results[eos_mac_columns].mean().mean(),
            'max_strength': results['max_value'].max(),
            'total_strength': results[eos_mac_columns].sum().sum()
        }
    return pd.DataFrame(analysis).T

interaction_analysis = analyze_interaction_strengths(condition_results)
print(interaction_analysis)

# Visualizations

# 1. Number of interactions by condition
plt.figure(figsize=(12, 6))
sns.barplot(x=interaction_analysis.index, y='num_interactions', data=interaction_analysis)
plt.title('Number of Eosinophil-Macrophage Interactions by Condition')
plt.xlabel('Condition')
plt.ylabel('Number of Interactions')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('num_interactions_by_condition.png')
plt.close()

# 2. Average interaction strength by condition
plt.figure(figsize=(12, 6))
sns.barplot(x=interaction_analysis.index, y='avg_strength', data=interaction_analysis)
plt.title('Average Eosinophil-Macrophage Interaction Strength by Condition')
plt.xlabel('Condition')
plt.ylabel('Average Interaction Strength')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('avg_strength_by_condition.png')
plt.close()

# 3. Heatmap of top interactions across conditions
top_interactions = pd.DataFrame({condition: results['interacting_pair'].head(20).tolist() 
                                 for condition, results in condition_results.items()})
all_top_interactions = pd.Series(top_interactions.values.ravel()).unique()

heatmap_data = pd.DataFrame(index=all_top_interactions, columns=conditions)
for condition, results in condition_results.items():
    heatmap_data[condition] = results.set_index('interacting_pair')['max_value']

heatmap_data = heatmap_data.fillna(0)

plt.figure(figsize=(20, 30))
sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd')
plt.title('Top Eosinophil-Macrophage Interactions Across Conditions')
plt.xlabel('Condition')
plt.ylabel('Interacting Pair')
plt.tight_layout()
plt.savefig('top_interactions_heatmap.png')
plt.close()

# 4. Scatter plot of number of interactions vs average strength
plt.figure(figsize=(12, 8))
sns.scatterplot(x='num_interactions', y='avg_strength', hue=interaction_analysis.index, 
                data=interaction_analysis.reset_index(), s=100)
plt.title('Number of Interactions vs Average Strength by Condition')
plt.xlabel('Number of Interactions')
plt.ylabel('Average Interaction Strength')
for i, row in interaction_analysis.reset_index().iterrows():
    plt.annotate(row['index'], (row['num_interactions'], row['avg_strength']))

plt.tight_layout()
plt.savefig('interactions_vs_strength_scatter_by_condition.png')
plt.close()

print("Analysis complete. Check the generated PNG files for visualizations.")
