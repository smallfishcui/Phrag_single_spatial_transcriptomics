seurat_obj_all<-readRDS("seurat_obj_all_harmony_res0.2_0.8.rds")
EU <- c("EU60","EU78","EU620")
seurat_obj_all$group <- sapply(seurat_obj_all$orig.ident, function(ita) ifelse(ita %in% EU,"EU","invasive"))
options(future.globals.maxSize = 10000 * 1024^2)  #
#seurat_obj_all <- SCTransform(seurat_obj_all, vst.flavor = "v2", verbose = TRUE)

#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
# Load libraries
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(seurat_obj_all, layer = "counts")
metadata <- seurat_obj_all@meta.data
# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(Idents(seurat_obj_all))

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names
# Total number of clusters
length(cluster_names)
# Extract unique names of samples (= levels of sample_id factor variable)
# Ensure orig.ident is a character variable
colData(sce)$orig.ident <- as.character(colData(sce)$orig.ident)
# Extract unique sample names
sample_names <- unique(colData(sce)$orig.ident)
# Print unique sample names
sample_names
# Total number of samples
length(sample_names)
groups <- colData(sce)[, c("cluster_id", "orig.ident")]
head(groups)
# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]
# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]
## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
# Reminder: explore structure of metadata
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(group, orig.ident)
# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
rownames(metadata) <- metadata$orig.ident
head(metadata)
# Number of cells per sample and cluster
t <- table(colData(sce)$orig.ident,
           colData(sce)$cluster_id)
t[1:6, 1:6]

# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names
# Loop over all cell types to extract corresponding counts, and store information in a list
## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {

  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]

}

# Explore the different components of the list
str(counts_ls)

# Creating metadata list
## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$cluster_id))
    cell_counts <- t[, idx]
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample_id, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    #df <- plyr::join(df, metadata, 
    #                 by = intersect(names(df), names(metadata)))
    df <- merge(df, metadata, by.x = "sample_id", by.y = "orig.ident", all.x = TRUE)
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$cluster_sample_id
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)

}

# Explore the different components of the list
str(metadata_ls)

idx <- which(names(counts_ls) == "13")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]
cluster_metadata$group <- as.factor(cluster_metadata$group)
# Verify it's now a factor
str(cluster_metadata)
table(cluster_metadata$group, useNA = "always")  # Should not
# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# Plot PCA
pdf("clusterall_pca_deseq2.pdf")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")
dev.off()
pdf("cluster18_pca_cellcount_deseq2.pdf")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
dev.off()
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
# Plot heatmap
pdf("cluster13_heatmap.pdf")
pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop=F])
dev.off()

#PCAfor all clusters
# 提取所有 cluster 的 metadata 数据
all_metadata <- do.call(rbind, lapply(names(metadata_ls), function(cluster) {
  metadata_ls[[cluster]]
}))
# 提取所有 cluster 的 counts 数据
all_counts <- do.call(cbind, lapply(names(counts_ls), function(cluster) {
  counts_ls[[cluster]]
}))
# 按样本聚合 metadata 数据
all_metadata_aggregated <- all_metadata[!duplicated(all_metadata$sample_id), ]
rownames(all_metadata_aggregated) <- all_metadata_aggregated$sample_id
# 检查列名和行名是否一致
identical(colnames(all_counts_aggregated), rownames(all_metadata_aggregated))
# 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = all_counts_aggregated,
  colData = all_metadata_aggregated,
  design = ~ group
)
# 运行 DESeq2
dds <- DESeq(dds)
# 提取结果
res <- results(dds)
summary(res)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# Plot PCA
pdf("clusterall_pca_deseq2.pdf")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")
dev.off()
#
# Run the script on all clusters comparing stimulated
get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
  print(clustx) # useful for debugging
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  # Ensure group is a factor with correct levels
  cluster_metadata$group <- factor(cluster_metadata$group, levels = c(A, B))
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group)
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  # Generate QC plots
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "sample_id")
  if (!dir.exists("results")) { dir.create("results") }
  ggsave(paste0("results_allsamples/", clustx, "_specific_PCAplot.png"))
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  ## Plot and save heatmap
  png(paste0("results_allsamples/", clustx, "_specific_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
    pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop = FALSE])
  dev.off()
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  ## Plot dispersion estimates
  png(paste0("results_allsamples/", clustx, "_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
    plotDispEsts(dds)
  dev.off()
  ## Output and shrink results of Wald test for contrast A vs B
  contrast <- paste(c("group", A, "vs", B), collapse = "_")
  # resultsNames(dds)
  res <- results(dds, name = resultsNames(dds)[2], alpha = 0.05)
  res <- lfcShrink(dds, coef = resultsNames(dds)[2], res = res)
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  #write.csv(res_tbl,
  #          paste0("results/", clustx, "_", resultsNames(dds)[2], "_all_genes.csv"),
   #         quote = FALSE, 
   #         row.names = FALSE)
## Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff & abs(log2FoldChange) > 2) %>%
    dplyr::arrange(padj)
  write.csv(sig_res,
            paste0("results_allsamples/", clustx, "_", resultsNames(dds)[2], "_signif_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  # Generate results visualization plots
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  ## Extract top 20 DEG from resLFC (make sure to order by padj)
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)
  ## Extract matching normalized count values from matrix
  top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]