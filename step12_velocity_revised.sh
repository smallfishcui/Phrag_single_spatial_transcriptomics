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
import loompy

adata= sc.read('seurat_obj_EU620_harmony_clusters_0.8.h5ad')
ldata = sc.read_loom("/scratch/cuiwang/nodelete5/singlecellanalysis/velocyte/E620/EU620.loom")
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

############revise for other individuals#####################
adata= sc.read('seurat_obj_NAint113_harmony_clusters_0.8.h5ad')
ldata = sc.read_loom("/scratch/cuiwang/nodelete5/singlecellanalysis/velocyte/NAint113/NAint113.loom")
#rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
barcodes = ['NAint113_' + bc[0:len(bc)-1] for bc in barcodes]
ldata.obs.index = barcodes
ldata.var.index=ldata.var.Accession

barcodes = [bc.split('-')[0] for bc in adata.obs.index.tolist()]
substring_to_replace1 = '1_'
replacement_substring1 = 'NAint113_'
modified_list = [bc.replace(substring_to_replace1, replacement_substring1) for bc in barcodes]
adata.obs.index=modified_list

# set the corrsponding index for the variable names between ladata and adata
adata.var.index=adata.var.index.tolist()
modified_var = [bc.replace('-', '_') for bc in adata.var.index]
adata.var.index=modified_var
ldata.var_names_make_unique()

# Check how many barcodes match between ldata and adata
common_barcodes = set(ldata.obs.index).intersection(set(adata.obs.index))
print(f"Number of common barcodes: {len(common_barcodes)}")
print(f"Total barcodes in ldata: {len(ldata.obs.index)}")
print(f"Total barcodes in adata: {len(adata.obs.index)}")

# Filter ldata to include only barcodes present in adata
ldata = ldata[ldata.obs.index.isin(common_barcodes)]
adata = scv.utils.merge(adata, ldata)
print(f"Number of cells in merged adata: {adata.n_obs}")

adata.write('NAint113_Allcelltype_dynamicModel.h5ad', compression = 'gzip')
adata= sc.read('NAint113_Allcelltype_dynamicModel.h5ad')
# Load the UMAP coordinates from the CSV file
umap_coords = pd.read_csv("umap_coordinates_harmony_clusters_0.8_NAint113.csv", index_col=0)
umap_coords.index = umap_coords.index.str.replace('-1$', '', regex=True)
print(sum(umap_coords.index.isin(adata.obs_names)), "/", len(umap_coords))
adata.obsm["X_umap"] = umap_coords.loc[adata.obs_names].values

########################for NAint61,barcode in loom file is in different format, 胞条形码示例: Index(['possorted_genome_bam_RKSXH:AAAGTCCTCTGTCTCGx',需要特殊处理
##############NAint191#######3# 去掉末尾的 "-1"
adata.obs.index = [bc.rstrip('-1') if bc.endswith('-1') else bc for bc in adata.obs.index]

# 或者使用 split
adata.obs.index = [bc.split('-')[0] for bc in adata.obs.index]

# 检查前 10 个
print(adata.obs.index[:10])