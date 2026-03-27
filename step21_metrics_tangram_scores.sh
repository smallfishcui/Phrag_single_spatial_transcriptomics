cd /scratch/spatial
# Plot spatial clusters with fixed orientation and adjusted spot size
import tangram as tg
import anndata
import scanpy as sc
import scvelo as scv
import torch
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  
from scipy.stats import pearsonr, spearmanr
import pandas as pd

sc_data = sc.read_h5ad('./seurat_obj_NAint61_harmony_clusters_0.8.h5ad')
spatial_data =sc.read_h5ad('./Spatial_obj_NAint61_spatial_whole.h5ad')
tangram_map = sc.read_h5ad('./NAint61_tangram_map_cell_level_L4_allgenetraning.h5ad')
tg.pp_adatas(sc_data, spatial_data)
ad_ge = tg.project_genes(adata_map=tangram_map, adata_sc=sc_data)

#len(markers) #1415 genes
df = pd.read_csv('marker_genes_noblank2_lower.txt')  
markers = df.iloc[:, 0].tolist()  
print(len(markers))
print(markers[:5]) 
spatial_data.var_names = [v.lower() for v in spatial_data.var_names]
ad_ge.var_names = [v.lower() for v in ad_ge.var_names]
markers = [m.lower() for m in markers]

results_data = []

exec('''
for marker in markers:
    if marker not in spatial_data.var_names:
        print(f"Warning: {marker} not found in spatial_data, skipping")
        continue
    if marker not in ad_ge.var_names:
        print(f"Warning: {marker} not found in ad_ge, skipping")
        continue
    if marker not in tangram_map.uns['train_genes_df'].index:
        print(f"Warning: {marker} not found in train_genes_df, skipping")
        continue

    measured_expr = spatial_data[:, marker].X.toarray().flatten()
    predicted_expr = ad_ge[:, marker].X.toarray().flatten()
    pearson_corr = pearsonr(measured_expr, predicted_expr)[0]
    spearman_corr = spearmanr(measured_expr, predicted_expr)[0]

    mapping_score = tangram_map.uns['train_genes_df'].loc[marker, 'train_score']

    results_data.append({
        'Gene': marker,
        'Mapping score': mapping_score,
        'Pearson correlation': pearson_corr,
        'Spearman correlation': spearman_corr
    })
''')

results_df = pd.DataFrame(results_data)
print(results_df)

results_df.to_csv('marker_correlation_results.csv', index=False)


