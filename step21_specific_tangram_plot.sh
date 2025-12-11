
module load  python-data/3.10-24.04
module load cuda
export PYTHONUSERBASE=/projappl/my-python-env
export PATH=$PATH:/projappl/my-python-env
# Set environment variable for OpenMP
export OMP_NUM_THREADS=32
sc_data = sc.read_h5ad('./seurat_obj_NAint61_harmony_clusters_0.8.h5ad')
spatial_data =sc.read_h5ad('./Spatial_obj_NAint61_spatial_whole_L4.h5ad')
tangram_map = sc.read_h5ad('./NAint61_tangram_map_L4_level.h5ad')
ad_map = sc.read_h5ad('./NAint61_tangram_map_cluster_splotLevel4_1000hvg.h5ad')
# Project 'seurat_clusters' from `ad_map.obs` to `spatial_data.var`
tg.project_cell_annotations(ad_map, spatial_data, annotation="harmony_clusters_0.8")
# Extract spatial coordinates and cluster labels
# Extract spatial coordinates and cluster labels
imagecol = spatial_data.obs['imagecol']
imagerow = spatial_data.obs['imagerow']
# Extract the most probable cluster for each spot
# Assign the most probable cluster to spatial_data.obs
spatial_data.obs['harmony_clusters_0.8'] = spatial_data.obsm['tangram_ct_pred'].idxmax(axis=1)
spatial_data.obsm['spatial'] = spatial_data.obs[[ 'imagecol','imagerow']].to_numpy()

sc.pl.spatial(
    spatial_data,
    color='harmony_clusters_0.8',  # Use the projected cluster labels
    spot_size=4,             # Adjust for clarity
    show=False 
)
# Save the figure with high resolution
plt.savefig("NAint61_tangram_map_cluster_splotLevel4_1000hvg_plot.png", dpi=300, bbox_inches='tight')
plt.show()

# Plot only cluster 1
# Create a new column for cluster 1 flag
spatial_data.obs['is_cluster_10'] = (spatial_data.obs['harmony_clusters_0.8'] == '10').astype(str)
# Plot where color reflects if a spot is in cluster 1
sc.pl.spatial(
    spatial_data,
    color='is_cluster_10',  # New column indicating cluster 1
    spot_size=4,           # Adjust spot size
    alpha_img=1.0,         # Make spatial image fully visible
    show=False
)

# Save the plot
plt.savefig("NAint61_cluster10_cellmode_flagged_spatial_plot.png", dpi=300, bbox_inches='tight')
plt.show()

ad_ge = tg.project_genes(adata_map=tangram_map, adata_sc=sc_data)
ad_ge
genes = ['evm.tu.pachr12.52'] #XCP
genes = ['evm.tu.pachr1.636'] #CIPK19
genes = ['evm.tu.pachr1.2491'] #H3.1
genes = ['evm.tu.pachr1.2835'] 
genes = ['evm.tu.pachr9.421'] 
genes = ['evm.tu.pachr5.780', 'evm.tu.pachr6.2072', 'evm.tu.pachr9.421', 'evm.tu.pachr9.1485', 'evm.tu.pachr14.961']
ad_map.uns['train_genes_df'].loc[genes]
tg.pp_adatas(spatial_data, sc_data, genes=genes)
tg.plot_genes_sc(genes, spatial_data, ad_ge,spot_size=5, scale_factor=0.1, perc=0.02)
# Save the plot as a file (e.g., PNG)
plt.savefig("marker_gene_pachr9.421_plots.png", dpi=300, bbox_inches='tight')
plt.close()  # Close the plot to free resources