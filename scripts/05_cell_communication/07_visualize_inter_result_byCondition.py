import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

conditions = ['WT_Male_HFD', 'WT_Female_HFD', 'KO_Male_HFD', 'KO_Female_HFD','WT_Female_HFDInfected', 'KO_Female_HFDInfected'] 

def process_condition_results(condition):
    # Load significant means
    sig_means = pd.read_csv(f"./cellphonedb_results_{condition}/statistical_analysis_significant_means.txt", sep="\t")
    
    # Identify eosinophil-macrophage interactions
    eos_mac_columns = [col for col in sig_means.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]
    
    # Filter and process data
    eos_mac_interactions = sig_means[sig_means[eos_mac_columns].notna().any(axis=1)].copy()
    eos_mac_interactions['max_value'] = eos_mac_interactions[eos_mac_columns].max(axis=1)
    eos_mac_interactions = eos_mac_interactions.sort_values('max_value', ascending=False)
    
    return eos_mac_interactions.head(20)  # Return top 20 interactions

# Process results for each condition
condition_results = {condition: process_condition_results(condition) for condition in conditions}

# Function to identify eosinophil-macrophage columns
def get_eos_mac_columns(df):
    return [col for col in df.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]

# Create a comparison heatmap
plt.figure(figsize=(80, 15))
for i, (condition, results) in enumerate(condition_results.items()):
    eos_mac_columns = get_eos_mac_columns(results)
    
    plt.subplot(1, len(condition_results), i+1)
    sns.heatmap(results[eos_mac_columns].T, annot=True, fmt='.2f', cmap='YlOrRd')
    plt.title(f"Top 20 Interactions - {condition}")
    plt.xlabel("Interacting Pairs")
    plt.ylabel("Cell Type Interactions")
    plt.xticks(range(len(results)), results['interacting_pair'], rotation=90, ha='right')

plt.tight_layout()
plt.savefig('condition_comparison_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

# Print top interactions for each condition
for condition, results in condition_results.items():
    print(f"\nTop 5 interactions for {condition}:")
    print(results[['interacting_pair', 'partner_a', 'partner_b', 'max_value']].head())


#  save separate heatmaps for each condition
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Function to identify eosinophil-macrophage columns
def get_eos_mac_columns(df):
    return [col for col in df.columns if 'Eosinophils|Macrophages' in col or 'Macrophages|Eosinophils' in col]

# Create separate heatmaps for each condition
for condition, results in condition_results.items():
    eos_mac_columns = get_eos_mac_columns(results)
    
    plt.figure(figsize=(20, 15))  # Adjust figure size as needed
    sns.heatmap(results[eos_mac_columns].T, annot=True, fmt='.2f', cmap='YlOrRd')
    plt.title(f"Top 20 Interactions - {condition}")
    plt.xlabel("Interacting Pairs")
    plt.ylabel("Cell Type Interactions")
    plt.xticks(range(len(results)), results['interacting_pair'], rotation=90, ha='right')
    
    plt.tight_layout()
    plt.savefig(f'heatmap_{condition}.png', dpi=300, bbox_inches='tight')
    plt.close()

print("Heatmaps saved for each condition.")
