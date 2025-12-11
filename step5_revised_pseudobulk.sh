seurat_obj_all<-readRDS("seurat_obj_all_harmony_res0.2_0.8.rds")
EU <- c("EU60","EU78","EU620")
seurat_obj_all$group <- sapply(seurat_obj_all$orig.ident, function(ita) ifelse(ita %in% EU,"EU","invasive"))
options(future.globals.maxSize = 100000 * 1024^2)  #
#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
# Load libraries
.libPaths(c("/projappl/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
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
library(muscat)
library(Seurat)
library(DESeq2)
library(tibble)
library(dplyr)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(seurat_obj_all, layer = "counts")
metadata <- seurat_obj_all@meta.data
Mean_reads_per_cell <- c(
  "EU60" = 30786,
  "EU78" = 29401,
  "EU620" = 59994,
  "NAint61" = 42478,
  "NAint113" = 24167,
  "NAint191" = 33900
)
metadata$mean_reads_per_cell <- Mean_reads_per_cell[metadata$orig.ident]
# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(Idents(seurat_obj_all))
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)
sce <- prepSCE(
  sce,
  kid = "cluster_id",   # cluster column name
  sid = "orig.ident",    # sample column name
  gid = "group",        # group column name
  drop = FALSE
)
pb <- aggregateData(
  sce,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)
metadata(pb)$experiment_info

metadata_df <- metadata(pb)$experiment_info
# 添加协变量
Mean_reads_per_cell <- c(
  "EU60" = 30786,
  "EU78" = 29401,
  "EU620" = 59994,
  "NAint61" = 42478,
  "NAint113" = 24167,
  "NAint191" = 33900
)
metadata_df$mean_reads_per_cell <- Mean_reads_per_cell[metadata_df$sample_id]

# 创建输出文件夹
output_dir <- "DEG_results_per_cluster"
if (!dir.exists(output_dir)) dir.create(output_dir)

cluster_ids <- names(assays(pb)) # muscat aggregateData生成的每个cluster的名字
p_list <- list()

for (cl in cluster_ids) {
  counts_matrix <- assay(pb, cl)
  
  # 对列名做匹配
  colnames(counts_matrix) <- metadata_df$sample_id
  rownames(metadata_df) <- metadata_df$sample_id
  metadata_df <- metadata_df[match(colnames(counts_matrix), metadata_df$sample_id), ]
  
  # 确保协变量是数值并标准化
  metadata_df$mean_reads_per_cell <- scale(as.numeric(metadata_df$mean_reads_per_cell))
  
  # 创建 DESeq2 对象
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = metadata_df,
    design = ~ mean_reads_per_cell + group_id
  )
##########################################################################################################3
rld <- rlog(dds, blind = FALSE)
pca_data <- plotPCA(rld, intgroup = c("group_id"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# 使用色盲友好的调色板 (Okabe-Ito)
color_palette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7"  # reddish purple
)
p <- ggplot(pca_data, aes(PC1, PC2, color = group_id, label = name)) +
    geom_point(size = 5) +
    geom_text_repel(size = 6, show.legend = FALSE) +
    scale_color_manual(values = color_palette) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    ggtitle(paste0("Cluster ", cl)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.8),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14)
    )
  
  p_list[[cl]] <- p
}
combined_plot <- wrap_plots(p_list, ncol = 2, guides = "collect") & theme(legend.position = "bottom")

# 保存大图
ggsave(file.path(output_dir, "All_clusters_PCA.png"),
       combined_plot,
       width = 16, height = 28, dpi = 300)  # 比 A4 长，16x28 英寸，300dpi

 #####################################################################################################3 
  # 过滤低表达基因
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  # 运行 DESeq2
  dds <- DESeq(dds)
  # lfcShrink
  res_shrunk <- lfcShrink(dds, coef = "group_id_invasive_vs_EU", type = "apeglm")
  
  # 转 tibble
  res_tbl <- res_shrunk %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  # 5. 筛选 padj < 0.05 且 |log2FoldChange| > 2 的基因
  #sig_genes <- res_tbl %>%
  #filter(padj < 0.05 & abs(log2FoldChange) > 1)
  #write.csv(sig_genes, file.path(output_dir, paste0("Cluster", cl, "_DEG_total_invasive_vs_EU_logFC1.csv")), row.names = FALSE, quote=FALSE)
  

  # 输出显著上调基因
  deg_up <- res_tbl %>%
    filter(padj < 0.05 & log2FoldChange > 2) %>%
    arrange(padj)
  write.csv(deg_up, file.path(output_dir, paste0("Cluster", cl, "_DEG_up_invasive_vs_EU_logFC1.csv")), row.names = FALSE, quote=FALSE)
  
  # 输出显著下调基因
  deg_down <- res_tbl %>%
    filter(padj < 0.05 & log2FoldChange < -2) %>%
    arrange(padj)
  write.csv(deg_down, file.path(output_dir, paste0("Cluster", cl, "_DEG_down_invasive_vs_EU_logFC1.csv")), row.names = FALSE, quote=FALSE)
  
  cat("Finished cluster:", cl, "\n")
}

#group into larger clusters
cluster_group_map <- list(
  parenchyma = c(2, 5, 9, 16),
  epidermis = c(6, 7, 8, 11),
  vascular = c(4, 10, 12, 13, 15, 18)
)
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(seurat_obj_all, layer = "counts")
metadata <- seurat_obj_all@meta.data
Mean_reads_per_cell <- c(
  "EU60" = 30786, "EU78" = 29401, "EU620" = 59994,
  "NAint61" = 42478, "NAint113" = 24167, "NAint191" = 33900
)
metadata$mean_reads_per_cell <- Mean_reads_per_cell[metadata$orig.ident]
metadata$cluster_id <- factor(Idents(seurat_obj_all))
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)
sce <- prepSCE(
  sce,
  kid = "cluster_id",   # cluster column name
  sid = "orig.ident",    # sample column name
  gid = "group",        # group column name
  drop = FALSE
)

sce$super_cluster <- sapply(sce$cluster_id, function(cl) {
  for (grp in names(cluster_group_map)) {
    if (as.numeric(cl) %in% cluster_group_map[[grp]]) return(grp)
  }
  return(NA)
})
sce$super_cluster <- factor(sce$super_cluster)

# --------------------------
# 4️⃣ 聚合到 super_cluster + sample
# --------------------------
pb <- aggregateData(
  sce,
  assay = "counts",
  fun = "sum",
  by = c("super_cluster", "sample_id")
)

metadata_df <- metadata(pb)$experiment_info
metadata_df$mean_reads_per_cell <- Mean_reads_per_cell[metadata_df$sample_id]
output_dir <- "DEG_results_per_supercluster"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --------------------------
# 6️⃣ 循环处理每个 super_cluster
# --------------------------
super_clusters <- names(assays(pb))  # parenchyma, epidermis, vascular
p_list <- list()

for (cl in super_clusters) {
  
  counts_matrix <- assay(pb, cl)
  
  # 匹配 colData
  colnames(counts_matrix) <- metadata_df$sample_id
  rownames(metadata_df) <- metadata_df$sample_id
  metadata_df <- metadata_df[match(colnames(counts_matrix), metadata_df$sample_id), ]
  
  # 标准化协变量
  metadata_df$mean_reads_per_cell <- scale(as.numeric(metadata_df$mean_reads_per_cell))
  
  # 创建 DESeq2 对象
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = metadata_df,
    design = ~ mean_reads_per_cell + group_id
  )
  
  # 过滤低表达基因
  dds <- dds[rowSums(counts(dds)) > 10, ]

  # rlog 变换
  rld <- rlog(dds, blind = FALSE)

  # PCA 数据
  pca_data <- plotPCA(rld, intgroup = c("group_id"), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  # 色盲友好调色板（Okabe–Ito）
  color_palette <- c(
    "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )

  # PCA 图
  p <- ggplot(pca_data, aes(PC1, PC2, color = group_id, label = name)) +
      geom_point(size = 5) +
      geom_text_repel(size = 6, show.legend = FALSE) +
      scale_color_manual(values = color_palette) +
      xlab(paste0("PC1: ", percentVar[1], "%")) +
      ylab(paste0("PC2: ", percentVar[2], "%")) +
      ggtitle(paste0("Super-cluster: ", cl)) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14)
      )

  p_list[[cl]] <- p

  cat("Finished PCA for super_cluster:", cl, "\n")
}

# 合并所有 super-cluster 的 PCA 图
combined_plot <- wrap_plots(p_list, ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# 保存
ggsave(
  file.path(output_dir, "All_superclusters_PCA.png"),
  combined_plot,
  width = 16, height = 28, dpi = 300
)

  # 运行 DESeq2
  dds <- DESeq(dds)
  
  # lfcShrink
  res_shrunk <- lfcShrink(dds, coef = "group_id_invasive_vs_EU", type = "apeglm")
  
  # 转 tibble
  res_tbl <- res_shrunk %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    as_tibble()
  
  # 输出显著上调基因
  deg_up <- res_tbl %>% filter(padj < 0.05 & log2FoldChange > 1) %>% arrange(padj)
  write.csv(
    deg_up,
    file.path(output_dir, paste0(cl, "_DEG_up_invasive_vs_EU_logFC1.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  # 输出显著下调基因
  deg_down <- res_tbl %>% filter(padj < 0.05 & log2FoldChange < -1) %>% arrange(padj)
  write.csv(
    deg_down,
    file.path(output_dir, paste0(cl, "_DEG_down_invasive_vs_EU_logFC1.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("Finished super_cluster:", cl, "\n")
}
#########################################################plot

