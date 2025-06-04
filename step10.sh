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
import cellrank as cr

cr.settings.verbosity = 2
sc.settings.set_figure_params(
    frameon=False,
    dpi=300,
)
import warnings
warnings.simplefilter("ignore", category=UserWarning)
#adata= scv.read('Allcelltype_dynamicModel.h5ad')
#adata= scv.read('NAint61_Allcelltype_dynamicModel.h5ad')
adata= sc.read('NAint61_Allcelltype_dynamicModel.h5ad')
# Load the UMAP coordinates from the CSV file
umap_coords = pd.read_csv("umap_coordinates_harmony_clusters_0.8_NAint61.csv", index_col=0)
# Assign the UMAP coordinates to the AnnData object
adata.obsm["X_umap"] = umap_coords.values
scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(adata, n_jobs=40)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
sc.pl.embedding(adata, basis="umap", color=["harmony_clusters_0.8"])
plt.savefig("projection_umap_plot.png", dpi=300, bbox_inches="tight")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
vk.plot_projection(color="harmony_clusters_0.8")
# Save the plot to a file
plt.savefig("projection_plot.png", dpi=300, bbox_inches="tight")
##########slingshot in R

########
sc.tl.diffmap(adata)
# Step 2: Get all cells in 'Cluster_14'
cluster_cells = adata.obs.index[adata.obs['harmony_clusters_0.8'] == '14']
# Step 2: Find the integer indices corresponding to the cell IDs
# We need to find the positions of `cluster_cells` within `adata.obs.index`, which will give us integer indices
cluster_cell_indices = np.where(np.isin(adata.obs.index, cluster_cells))[0]
# Step 3: Calculate the centroid (average position) of the cluster in the diffusion space
diffmap_coords = adata.obsm['X_diffmap'][cluster_cell_indices]
centroid = np.mean(diffmap_coords, axis=0)
# Step 5: Find the closest cell to the centroid
distances = np.linalg.norm(diffmap_coords - centroid, axis=1)
closest_cell_index = cluster_cell_indices[np.argmin(distances)]
root_cell_index = closest_cell_index
# Step 1: Set the root cell index in adata.uns['iroot']
adata.uns['iroot'] = root_cell_index
# Step 2: Run pseudotime calculation
scv.pl.scatter(
    adata,
    basis="diffmap",
    c=["harmony_clusters_0.8", root_cell_index],
    legend_loc="right",
    components=["2, 3"]
)
plt.savefig("scatter_diffmap_plot_2_3.png", dpi=300, bbox_inches="tight")
sc.tl.dpt(adata)
root_cell_id = adata.obs.index[root_cell_index]
import palantir
#dm_res = palantir.utils.run_diffusion_maps(adata, n_components=30)
#ms_data = palantir.utils.determine_multiscale_space(adata)
#imputed_X = palantir.utils.run_magic_imputation(adata)
#diffmap_coords = adata.obsm['X_diffmap']
# 按照第三列（[:, 2]）从小到大排序
#sorted_indices = diffmap_coords[:, 3].argsort()
# 获取排序后**第三个**（索引为2的）细胞的索引
##third_leftmost_cell_index = sorted_indices[3]  # 取第3个最左侧细胞
##root_ixs = root_cell_index
#adata.uns["iroot"] = root_ixs
#sc.tl.dpt(adata)
#early_cell_name = adata.obs_names[root_cell_index]
#early_cells = palantir.utils.early_cell(adata, "14", "harmony_clusters_0.8")
#early_cell_name = adata.obs_names[root_ixs]
#pr_res = palantir.core.run_palantir(adata, early_cell=early_cell_name, num_waypoints=1000)
#adata.obs['palantir_pseudotime'] = pr_res.pseudotime

#sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime","palantir_pseudotime"],
    color_map="gnuplot2",
)
#plt.savefig("dpt_palantir_pseudotime_plot.png", dpi=300, bbox_inches="tight")
#adata.obs["palantir_pseudotime"] = pr_res.pseudotime
#palantir.compute_pseudotime(adata, root_cell=root_cell_id)

mono_trajectory = ["14", "15", "2",  "12", "18"]
ery_trajectory = ["14", "15", "11"]
mask = np.in1d(adata.obs["harmony_clusters_0.8"], mono_trajectory)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime", "palantir_pseudotime"],
    groupby="harmony_clusters_0.8",
    rotation=-90,
    order=mono_trajectory,
)
plt.savefig("violin_plot_mono_trajectory.png", dpi=300, bbox_inches="tight")
mask = np.in1d(adata.obs["harmony_clusters_0.8"], ery_trajectory)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime", "palantir_pseudotime"],
    groupby="harmony_clusters_0.8",
    rotation=-90,
    order=ery_trajectory,
)
plt.savefig("violin_plot_ery_trajectory.png", dpi=300, bbox_inches="tight")

from cellrank.kernels import PseudotimeKernel
pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()
pk.plot_random_walks(
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
plt.savefig("pseudotime_umap_plot.png", dpi=300, bbox_inches="tight")

pk.plot_projection(basis="umap", recompute=True)
plt.savefig("stream_umap_plot.png", dpi=300, bbox_inches="tight")

from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * pk + 0.2 * ck
combined_kernel
pk.write_to_adata()
#automatically calculate the number of states between 5-22
# 创建 GPCCA 对象并拟合模型
from cellrank import estimators as cr_estimators  # 导入 estimators 模块
g = cr_estimators.GPCCA(combined_kernel)  # 使用联合核
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