library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
##载入脚本
source("CreateBmkObject.R")
##创建对象
spatial<- CreateS1000Object(matrix_path="./L4_heAuto", png_path="./he_roi_small.png", spot_radius =NULL, min.features = 50) #important, default min.cells = 5,min.features = 100 
spatial@meta.data$imagerow <- spatial@images$sample1@coordinates$imagerow
spatial@meta.data$imagecol <- spatial@images$sample1@coordinates$imagecol
# Extract the raw count matrix from the Seurat v5 object
counts <- t(GetAssayData(spatial@assays$Spatial, layer = "counts"))
metadata <- spatial@meta.data
gene_names <- rownames(spatial@assays$Spatial)
dim(counts) # Should show (6868, 28606) if there are 6868 cells and 28606 genes
dim(metadata) # Should also show (6868, n) where n is the number of metadata columns
length(gene_names) #
rownames(metadata) <- as.character(rownames(metadata))
# Ensure gene names are strings
gene_names <- as.character(gene_names)
# Create the AnnData object with explicitly set observation and variable indices
anndata <- reticulate::import("anndata")
adata <- anndata$AnnData(
  X = counts,
  obs = metadata,
  var = data.frame(gene_names = gene_names, row.names = gene_names)
)
# Explicitly set var_names in the AnnData object
adata$var_names = gene_names

# Save the AnnData object as a .h5ad file
adata$write("Spatial_obj_NAint61_spatial_whole.h5ad")
module load  python-data/3.10-24.04
module load cuda
export PYTHONUSERBASE=/projappl/project_2009273/my-python-env
export PATH=$PATH:/projappl/project_2009273/my-python-env
# Set environment variable for OpenMP
export OMP_NUM_THREADS=32
# Set PyTorch memory management config
srun --gres=gpu:1 --exclusive --cpu-bind=none python3 sc_sp_tangram.py
#python script
import tangram as tg
import anndata
import scanpy as sc
import scvelo as scv
import torch
torch.cuda.empty_cache()
# Load single-cell data
sc_data = sc.read_h5ad('./seurat_obj_EU620_harmony_clusters_0.8.h5ad')
# Identify highly variable genes
import numpy as np
import scipy.sparse as sp
# Load spatial data
spatial_data = scv.read('./Spatial_obj_NAint61_spatial_whole.h5ad')
spatial_data.obsm['spatial'] = spatial_data.obs[['imagerow', 'imagecol']].to_numpy()
# Align genes between datasets
tg.pp_adatas(sc_data, spatial_data)
# Train Tangram
tangram_map = tg.map_cells_to_space(
    adata_sc=sc_data,
    adata_sp=spatial_data,
    device='cuda'  # Use 'cpu' if no GPU is availabl
)
# Save the cell-level mapping result
tangram_map.write('./tangram_map_cell_level.h5ad')
# Project single-cell data,cluster level
ad_map = tg.map_cells_to_space(
                   sc_data, 
                   spatial_data,         
                   mode='clusters',
                   cluster_label='harmony_clusters_0.8')
# Save the cluster-level mapping result
ad_map.write('./tangram_map_cluster_level.h5ad')

###############################plot

# Plot spatial clusters with fixed orientation and adjusted spot size
import tangram as tg
import anndata
import scanpy as sc
import scvelo as scv
import torch
import matplotlib.pyplot as plt

sc_data = sc.read_h5ad('./seurat_obj_EU620_harmony_clusters_0.8.h5ad')
spatial_data =sc.read_h5ad('./Spatial_obj_NAint61_spatial_whole.h5ad')
tangram_map = sc.read_h5ad('./tangram_map_cell_level.h5ad')
ad_map = sc.read_h5ad('./tangram_map_cluster_level.h5ad')
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
plt.savefig("spatial_clusters_plot.png", dpi=300, bbox_inches='tight')
plt.show()

ad_ge = tg.project_genes(adata_map=tangram_map, adata_sc=sc_data)
ad_ge
genes = ['evm.tu.pachr12.52', 'evm.tu.pachr1.636', 'evm.tu.pachr1.2491', 'evm.tu.pachr1.2835', 'evm.tu.pachr5.780', 'evm.tu.pachr6.2072', 'evm.tu.pachr9.421', 'evm.tu.pachr9.1485', 'evm.tu.pachr14.961']
ad_map.uns['train_genes_df'].loc[genes]
tg.pp_adatas(spatial_data, sc_data, genes=genes)
tg.plot_genes_sc(genes, spatial_data, ad_ge,spot_size=5, scale_factor=0.1, perc=0.02)
# Save the plot as a file (e.g., PNG)
plt.savefig("marker_gene_plots.png", dpi=300, bbox_inches='tight')
plt.close()  # Close the plot to free resources


# Plot only cluster 1
# Create a new column for cluster 1 flag
spatial_data.obs['is_cluster_4'] = (spatial_data.obs['harmony_clusters_0.8'] == '4').astype(str)
# Plot where color reflects if a spot is in cluster 1
sc.pl.spatial(
    spatial_data,
    color='is_cluster_4',  # New column indicating cluster 1
    spot_size=4,           # Adjust spot size
    alpha_img=1.0,         # Make spatial image fully visible
    show=False
)

# Save the plot
plt.savefig("cluster4_flagged_spatial_plot.png", dpi=300, bbox_inches='tight')
plt.show()