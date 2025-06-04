library(reticulate)
reticulate::py_install("anndata", pip = TRUE)
# 导入anndata包
anndata <- import("anndata")
#########
# Extract the raw count matrix from the Seurat v5 object
counts <- t(GetAssayData(seurat_obj_all.list[[3]]@assays$RNA, layer = "counts"))
metadata <- seurat_obj_all.list[[3]]@meta.data
gene_names <- rownames(seurat_obj_all.list[[3]]@assays$RNA)
# Check that the number of rows in counts matches the number of rows in metadata
dim(counts) # Should show (6456, 28606) if there are 6868 cells and 28606 genes
dim(metadata) # Should also show (6456, n) where n is the number of metadata columns
length(gene_names) # Should be 28606
# Now create the AnnData object
# Now create the AnnData object with gene names
# Ensure the row names of metadata are strings
rownames(metadata) <- as.character(rownames(metadata))
# Ensure gene names are strings
gene_names <- as.character(gene_names)
# Assuming your Seurat object is named `seurat_obj`
umap_coordinates <- Embeddings(seurat_obj_all.list[[3]], reduction = "harmony_clusters_0.8")
#umap_coordinates <- Embeddings(seurat_obj_all.list[[3]], reduction = "umap.unintegrated")
write.csv(umap_coordinates, "umap_coordinates_harmony_clusters_0.8.csv")
# Create the AnnData object with explicitly set observation and variable indices
adata <- anndata$AnnData(
  X = counts,
  obs = metadata,
  var = data.frame(gene_names = gene_names, row.names = gene_names)
)
# Explicitly set var_names in the AnnData object
adata$var_names = gene_names
################optional, be careful to use, remove all columns which are problematic with chracters, use when there is error with saving############
for (col in colnames(adata$obs)) {
    tryCatch({
        adata$obs[[col]] <- as.character(adata$obs[[col]])
        adata$obs[[col]][is.na(adata$obs[[col]])] <- "NA"
    }, error = function(e) {
        cat("Removing problematic column:", col, "\n")
        adata$obs[[col]] <- NULL
    })
}

# Save the AnnData object as a .h5ad file
adata$write("seurat_obj_EU620_harmony_clusters_0.8.h5ad")