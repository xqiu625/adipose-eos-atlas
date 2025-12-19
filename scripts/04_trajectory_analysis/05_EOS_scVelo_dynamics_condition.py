#!/usr/bin/env python3
import os
import scvelo as scv
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === 设置路径 === 
input_dir = "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/"
output_dir = "/bigdata/nairlab/shared/XQiu/EOS_project/R08_20250716/EOS_Trajectory/"
os.makedirs(output_dir, exist_ok=True)

# === EOS 文件列表 ===
eos_files = [
    "EOS_scVelo_adata_male_ko_naive.h5ad",
    "EOS_scVelo_adata_female_ko_naive.h5ad", 
    "EOS_scVelo_adata_female_ko_infection.h5ad",
    "EOS_scVelo_adata_male_wt_naive.h5ad",
    "EOS_scVelo_adata_female_wt_naive.h5ad",
    "EOS_scVelo_adata_female_wt_infection.h5ad",
    "EOS_scVelo_adata_infected.h5ad",
    "EOS_scVelo_adata_naive.h5ad", 
    "EOS_scVelo_adata_male.h5ad",
    "EOS_scVelo_adata_female.h5ad",
    "EOS_scVelo_adata_wt.h5ad",
    "EOS_scVelo_adata_ko.h5ad"
]

# === Cluster 配色 ===
eos_palette = {
    # Direct cluster names (C0, C1, etc.)
    "C0": "#6CA6CD",
    "C1": "#E84A5F", 
    "C2": "#FF6F1F",
    "C3": "#CBAACB",
    "C4": "#B03A48",
    "C5": "#5BA4A4",
    "C6": "#2F3C7E",
    # Also keep Eos_ prefixed versions as backup
    "Eos_C0": "#6CA6CD",
    "Eos_C1": "#E84A5F", 
    "Eos_C2": "#FF6F1F",
    "Eos_C3": "#CBAACB",
    "Eos_C4": "#B03A48",
    "Eos_C5": "#5BA4A4",
    "Eos_C6": "#2F3C7E"
}

# === scVelo 设置 ===
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# === 初始化总细胞数统计 ===
total_eos_cells = 0

# === 逐个处理文件 ===
for filename in eos_files:
    print(f"\n🧪 Processing {filename}")
    path = os.path.join(input_dir, filename)
    
    try:
        adata = scv.read(path)
        print(f"  - Shape: {adata.shape}")
        total_eos_cells += adata.n_obs
        
        # 确保有 velocity 信息
        if 'velocity' not in adata.layers:
            print("  🔄 Running velocity analysis...")
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            scv.tl.recover_dynamics(adata)
            scv.tl.velocity(adata, mode="dynamical")
            scv.tl.velocity_graph(adata)
        
        # 确保有 UMAP
        if 'X_umap' not in adata.obsm:
            print("  🔄 Computing UMAP...")
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
            sc.tl.umap(adata)
        
        # 处理 Cluster 信息
        cluster_col = None
        for col in ['Clusters', 'clusters', 'cluster', 'seurat_clusters', 'leiden']:
            if col in adata.obs.columns:
                cluster_col = col
                break
        
        if cluster_col is None:
            # 尝试从 obsm 中获取
            if 'X_cluster' in adata.obsm:
                adata.obs['Clusters'] = pd.Series(
                    adata.obsm['X_cluster'].flatten(), 
                    index=adata.obs_names
                ).astype(str).astype("category")
                cluster_col = 'Clusters'
            else:
                print("  ⚠️ No cluster info found, performing leiden clustering...")
                sc.tl.leiden(adata, resolution=0.5)
                cluster_col = 'leiden'
        
        # 确保 cluster 列是 categorical
        if not pd.api.types.is_categorical_dtype(adata.obs[cluster_col]):
            adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')
        
        # 设置颜色 - 关键修复
        categories = adata.obs[cluster_col].cat.categories
        color_key = f"{cluster_col}_colors"
        
        # 创建颜色列表
        colors = []
        for cat in categories:
            if str(cat) in eos_palette:
                colors.append(eos_palette[str(cat)])
            elif f"Eos_{cat}" in eos_palette:
                colors.append(eos_palette[f"Eos_{cat}"])
            else:
                # 为未定义的cluster生成颜色
                colors.append(plt.cm.tab10(len(colors) % 10))
        
        adata.uns[color_key] = colors
        
        print(f"  - Using cluster column: {cluster_col}")
        print(f"  - Categories found: {list(categories)}")
        print(f"  - Colors assigned: {len(colors)}")
        
        # 创建图形 - 分开绘制以确保颜色正确
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
        # 左图: UMAP with clusters
        scv.pl.scatter(
            adata, 
            basis='umap', 
            color=cluster_col,
            legend_loc='right margin',
            size=60,
            alpha=0.8,
            ax=ax1,
            show=False,
            title=f"{filename.replace('.h5ad', '')} - Clusters"
        )
        
        # 右图: Velocity stream plot
        scv.pl.velocity_embedding_stream(
            adata, 
            basis='umap', 
            color=cluster_col,
            legend_loc='none',  # 避免重复图例
            density=2,
            ax=ax2,
            show=False,
            title=f"{filename.replace('.h5ad', '')} - RNA Velocity"
        )
        
        plt.tight_layout()
        
        # 保存图片
        out_file = os.path.join(output_dir, f"{filename.replace('.h5ad', '')}_velocity.png")
        plt.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  ✅ Saved: {out_file}")
        
        # 额外保存单独的 velocity 图
        plt.figure(figsize=(12, 10))
        scv.pl.velocity_embedding_stream(
            adata, 
            basis='umap', 
            color=cluster_col,
            legend_loc='right margin',
            density=2,
            show=False,
            title=f"{filename.replace('.h5ad', '')} - RNA Velocity"
        )
        
        velocity_file = os.path.join(output_dir, f"{filename.replace('.h5ad', '')}_velocity_only.png")
        plt.savefig(velocity_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # 打印 cluster 分布
        print("  📊 Cluster counts:")
        print(adata.obs[cluster_col].value_counts())
        
        # 保存处理后的数据
        processed_file = os.path.join(output_dir, f"{filename.replace('.h5ad', '')}_processed.h5ad")
        adata.write(processed_file)
        
    except Exception as e:
        print(f"  ❌ Failed to process {filename}: {e}")
        import traceback
        traceback.print_exc()

# === 总结 ===
print("\n" + "="*50)
print(f"🎯 Total EOS cells across all files: {total_eos_cells}")
print(f"📂 Plots saved to: {output_dir}")
print("="*50)