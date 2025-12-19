import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
sig_means = pd.read_csv('cellphonedb_results/statistical_analysis_significant_means_08_17_2024_014813.txt', sep='\t')


# Identify columns for eosinophil-macrophage interactions
eos_mac_columns = [col for col in sig_means.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]

# Create a new DataFrame with only the relevant columns
columns_to_keep = ['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b'] + eos_mac_columns
eos_mac_interactions = sig_means[columns_to_keep].copy()

# Filter for rows where there's a non-NaN value in any of the eosinophil-macrophage columns
eos_mac_interactions = eos_mac_interactions[eos_mac_interactions[eos_mac_columns].notna().any(axis=1)]

# Calculate the maximum value across eosinophil-macrophage columns
eos_mac_interactions['max_value'] = eos_mac_interactions[eos_mac_columns].max(axis=1)

# Sort by the maximum value
eos_mac_interactions = eos_mac_interactions.sort_values('max_value', ascending=False)


# Display top interactions
print(eos_mac_interactions[['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'max_value']].head(20))

# Visualize top interactions
top_20 = eos_mac_interactions.head(20)

plt.figure(figsize=(15, 12))
sns.heatmap(top_20[eos_mac_columns].T, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Interaction Score'})

# Customize x-axis labels
plt.xticks(range(len(top_20)), top_20['interacting_pair'], rotation=90, ha='right')
plt.xlabel('Interacting Pairs')
plt.ylabel('Cell Type Interactions')

plt.title('Top 20 Eosinophil-Macrophage Interactions', fontsize=16)
plt.tight_layout()
plt.savefig('top_20_eos_mac_interactions.png', dpi=300, bbox_inches='tight')
plt.close()

