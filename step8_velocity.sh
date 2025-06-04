velocyto run --bcfile /scratch/Cui2/singlecellanalysis/EU620/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/E620 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/EU620/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 
velocyto run --bcfile /scratch/Cui2/singlecellanalysis/EU60/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/EU60 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/EU60/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 
velocyto run --bcfile /scratch/Cui2/singlecellanalysis/EU78/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/EU78 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/EU78/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 
velocyto run --bcfile /scratch/Cui2/singlecellanalysis/NAint61/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/NAint61 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/NAint61/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 
velocyto run --bcfile /scratch/Cui2/singlecellanalysis/NAint113/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/NAint113 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/NAint113/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 
velocyto run --bcfile /scratch/Cui2/singlecellanalysis/NAint191/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --outputfolder /scratch/Cui2/singlecellanalysis/velocyte/NAint191 --samtools-threads 30 --samtools-memory 8000 -vvv /scratch/Cui2/singlecellanalysis/NAint191/outs/possorted_genome_bam.bam PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf 

python
import scvelo as scv
import scanpy as sc
import numpy as np
import sklearn
import scipy
import pandas as pd
import seaborn as sns
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import igraph

adata = sc.read_h5ad("/scratch/project_2009273/Cui2/singlecellanalysis/seurat_obj_EU620_harmony_clusters_0.8.h5ad") #must use the full path
ldata = scv.read("/scratch/project_2009273/Cui2/singlecellanalysis/velocyte/E620/EU620.loom")

#rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
barcodes = ['EU620_' + bc[0:len(bc)-1] for bc in barcodes]
ldata.obs.index = barcodes
ldata.var.index=ldata.var.Accession

barcodes = [bc.split('-')[0] for bc in adata.obs.index.tolist()]
substring_to_replace1 = '1_'
replacement_substring1 = 'EU620_'
modified_list = [bc.replace(substring_to_replace1, replacement_substring1) for bc in barcodes]
adata.obs.index=modified_list

# set the corrsponding index for the variable names between ladata and adata
adata.var.index=adata.var.index.tolist()
modified_var = [bc.replace('-', '_') for bc in adata.var.index]
adata.var.index=modified_var

ldata.var_names_make_unique()

scv.utils.clean_obs_names(adata) 
scv.utils.clean_obs_names(ldata)

# Check how many barcodes match between ldata and adata
common_barcodes = set(ldata.obs.index).intersection(set(adata.obs.index))
print(f"Number of common barcodes: {len(common_barcodes)}")
print(f"Total barcodes in ldata: {len(ldata.obs.index)}")
print(f"Total barcodes in adata: {len(adata.obs.index)}")

# Filter ldata to include only barcodes present in adata
ldata = ldata[ldata.obs.index.isin(common_barcodes)]
adata = scv.utils.merge(adata, ldata)
print(f"Number of cells in merged adata: {adata.n_obs}")
adata.write('Allcelltype_dynamicModel.h5ad', compression = 'gzip')

adata= scv.read('Allcelltype_dynamicModel.h5ad')
# make Umap to check data
#sc.pp.pca(adata, n_comps=35)
#sc.pp.neighbors(adata, n_neighbors=50, use_rep='X_pca')
#sc.tl.umap(adata)

# Load the UMAP coordinates from the CSV file
umap_coords = pd.read_csv("umap_coordinates_harmony_clusters_0.8.csv", index_col=0)

# Assign the UMAP coordinates to the AnnData object
adata.obsm["X_umap"] = umap_coords.values
sc.pl.umap(adata, color='harmony_clusters_0.8', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')
#sc.pl.umap(adata, color='CellType', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

#get the proportion
#scv.pl.proportions(adata, groupby='CellType', save='_celltypes.pdf')
scv.pl.proportions(adata, groupby='harmony_clusters_0.8', save='_celltypes.pdf')
scv.pp.filter_and_normalize(adata)
# pre-process, deepvelo
scv.pp.moments(adata, n_pcs=30, n_neighbors=40)
# Extract dynamical velocity vectors
scv.tl.recover_dynamics(adata, n_jobs=40)
scv.tl.velocity(adata, mode='dynamical', use_raw=True)
scv.tl.velocity_graph(adata, n_jobs=40)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='harmony_clusters_0.8', save='velocity_stream_EU620.svg')
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='CellType_cluster', save='velocity_stream_EU620.svg')

# Latent time inference
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='latent_time_EU620.svg')

#https://colab.research.google.com/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb#scrollTo=GHjDOr1IBzP0

#PAGA
#scv.tl.rank_velocity_genes(adata, groupby='CellType', min_corr=.3)
scv.tl.rank_velocity_genes(adata, groupby='harmony_clusters_0.8', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

kwargs = dict(frameon=False, size=10, linewidth=1.5, add_outline='procambium, shoot apical meristem, intercalary meristem')
scv.pl.scatter(adata, df['shoot apical meristem'][:5], ylabel='shoot apical meristem', frameon=False, color='CellType', size=10, linewidth=1.5, save="shoot apical meristem_scatter.pdf")
#scv.pl.scatter(adata, df['Shoot system epidermis'][:5], ylabel='Protoxylem', frameon=False, color='CellType', size=10, linewidth=1.5, save='Protoxylem_scatter.pdf')
#scv.pl.scatter(adata, df['Shoot system epidermis'][:5], ylabel='Xylem', frameon=False, color='CellType', size=10, linewidth=1.5, save='Xylem_scatter.pdf')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save='velocity_confidence.pdf')

#scv.pl.velocity_graph(adata, threshold=.1, color='CellType',save='velocity_graph.pdf')
scv.pl.velocity_graph(adata, threshold=.1, color='harmony_clusters_0.8',save='velocity_graph.pdf')

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax,save='velocity_graph.pdf')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='velocity_pseudotime.pdf')

#scv.tl.paga(adata, groups='CellType')
scv.tl.paga(adata, groups='harmony_clusters_0.8')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save='velocity_paga.pdf')

#examine the low latent cells
# Define the threshold for low latent time
low_latent_threshold = adata.obs['latent_time'].quantile(0.05)  # bottom 5%

# Filter cells with low latent time
low_latent_cells = adata[adata.obs['latent_time'] <= low_latent_threshold]

# Print the number of cells with low latent time
print(f"Number of cells with low latent time: {low_latent_cells.n_obs}")

# Optional: Inspect the CellType distribution of low latent time cells
print("Cell type distribution in low latent time cells:")
print(low_latent_cells.obs['CellType'].value_counts())

# Visualize these cells on UMAP
scv.pl.scatter(low_latent_cells, color='CellType', basis='umap', size=80, 
               title='Low Latent Time Cells on UMAP', frameon=False, save='low_latent_cells.pdf')
# Save filtered AnnData object to H5AD format
low_latent_cells.write("low_latent_time_cells.h5ad")