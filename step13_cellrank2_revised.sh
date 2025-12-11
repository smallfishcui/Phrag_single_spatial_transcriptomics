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
import cellrank as cr

cr.settings.verbosity = 2
sc.settings.set_figure_params(
    frameon=False,
    dpi=300,
)
import warnings
warnings.simplefilter("ignore", category=UserWarning)
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

adata= sc.read('Allcelltype_dynamicModel.h5ad')
# Load the UMAP coordinates from the CSV file
umap_coords = pd.read_csv("umap_coordinates_harmony_clusters_0.8.csv", index_col=0)
umap_coords.index = umap_coords.index.str.replace('-1$', '', regex=True)
print(sum(umap_coords.index.isin(adata.obs_names)), "/", len(umap_coords))
adata.obsm["X_umap"] = umap_coords.loc[adata.obs_names].values

#follow https://www.sc-best-practices.org/trajectories/rna_velocity.html
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=22)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#scv.tl.velocity(adata, mode="deterministic")
#scv.tl.velocity_graph(adata, n_jobs=40)
#scv.pl.velocity_embedding_stream(adata, basis="umap", color=["harmony_clusters_0.8"])
#plt.savefig("deterministic_EU620_projection_umap_plot.png", dpi=300, bbox_inches="tight")
scv.tl.recover_dynamics(adata, n_jobs=40)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=40)
scv.pl.velocity_embedding_stream(adata, basis="umap", color=["harmony_clusters_0.8"])
plt.savefig("Dynamical_EU620_projection_umap_plot.png", dpi=300, bbox_inches="tight")
vk = cr.kernels.VelocityKernel(
    adata,
    mode="dynamical"        # 使用 dynamical 模型
)
#vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix(model='monte_carlo')
#ck = cr.kernels.ConnectivityKernel(adata)
#ck.compute_transition_matrix()
#combined_kernel = 0.8 * vk + 0.2 * ck
vk.plot_projection(color="harmony_clusters_0.8")
# Save the plot to a file
plt.savefig("Dynamical_EU620_projection_plot.png", dpi=300, bbox_inches="tight")
vk.write_to_adata()
adata.write(
    "./EU620_day1.11.2025_velocity_kernel.h5ad", compression="gzip"
)

############################333333velocity doesnt work, need to use pseudotime kernel
###########https://zhuanlan.zhihu.com/p/605287055
sc.tl.diffmap(adata)

# 选择根
SAM = np.isin(adata.obs['harmony_clusters_0.8'], ['14'])
min_stem_id = np.argmin(adata.obsm['X_diffmap'][SAM, 1])
root_id = np.arange(len(SAM))[SAM][min_stem_id]
adata.uns['iroot'] = root_id
sc.tl.dpt(adata)
sc.pl.umap(adata, color="dpt_pseudotime")
plt.savefig("DPT_EU620_projection_plot.png", dpi=300, bbox_inches="tight")
#adata.obsm['X_diffmap_'] = adata.obsm['X_diffmap'][:,1:]

from cellrank.kernels import PseudotimeKernel
ptk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
ptk.compute_transition_matrix()
ptk.plot_random_walks(
    seed=0,
    n_sims=100,
    start_ixs={"harmony_clusters_0.8": "14"},
    basis="umap",
    legend_loc="right",
    dpi=300,
)
plt.savefig("random_walks_cluster14.png", dpi=300, bbox_inches="tight")
import matplotlib.colors as mcolors
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime"],
    color_map="inferno",
    legend_loc="on data"
)
plt.savefig("pseudotime_umap_plotEU620.png", dpi=300, bbox_inches="tight")

ptk.plot_projection(basis="umap", recompute=True)
plt.savefig("stream_umap_plot_EU620.png", dpi=300, bbox_inches="tight")

from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.4 * vk + 0.1 * ck + 0.5 * ptk
combined_kernel
pk.write_to_adata()
#automatically calculate the number of states between 5-22
# 创建 GPCCA 对象并拟合模型
from cellrank import estimators as cr_estimators  # 导入 estimators 模块
g = cr_estimators.GPCCA(combined_kernel, random_state=42)  # 使用联合核
g.fit(n_states=[5, 22], cluster_key="harmony_clusters_0.8")
#g.fit(cluster_key="harmony_clusters_0.8")
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100)
plt.savefig("macrostates.png", dpi=300, bbox_inches="tight")
g.predict_terminal_states(method="top_n", n_states=10)
g.plot_macrostates(which="terminal")
plt.savefig("macrostates_terminal.png", dpi=300, bbox_inches="tight")
g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="right",same_plot=False)
plt.savefig("fate_probabilities.png", dpi=300, bbox_inches="tight")
cr.pl.circular_projection(adata, keys="harmony_clusters_0.8", legend_loc="right")
plt.savefig("circular_projection.png", dpi=300, bbox_inches="tight")
mono_drivers = g.compute_lineage_drivers(lineages="18")
mono_drivers.head(10)
model = cr.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["evm.TU.PaChr3.2674","evm.TU.PaChr9.74","evm.TU.PaChr12.68","evm.TU.PaChr21.930"],
    same_plot=True,
    ncols=2,
    time_key="dpt_pseudotime",
    hide_cells=True,
)
# Remove grid from all subplots
for ax in plt.gcf().get_axes():
    ax.grid(False)

plt.show()

plt.savefig("gene_trends_cluster18.png", dpi=300, bbox_inches="tight")



from cellrank import estimators as cr_estimators  # 导入 estimators 模块
g = cr_estimators.GPCCA(combined_kernel)  # 使用联合核
g.fit(cluster_key="harmony_clusters_0.8", n_states=[4, 8])

#g.fit(cluster_key="harmony_clusters_0.8")
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100)
plt.savefig("EU620_macrostates.png", dpi=300, bbox_inches="tight")



