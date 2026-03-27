cd /scratch/Cui3/singlecellanalysis
module load r-env
start-r
.libPaths(c("/projappl/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(Seurat)

seurat_obj_all<-readRDS("seurat_obj_all_harmony_res0.2_0.8.rds")
Idents(seurat_obj_all) <- seurat_obj_all@meta.data$harmony_clusters_0.8
EU <- c("EU60","EU78","EU620")
seurat_obj_all$group <- sapply(seurat_obj_all$orig.ident, function(ita) ifelse(ita %in% EU,"EU","invasive"))
options(future.globals.maxSize = 100000 * 1024^2)  #
library(Seurat)
library(ggplot2)

# 提取 harmony_clusters_0.8 的二维坐标
coord <- Embeddings(seurat_obj_all, "harmony_clusters_0.8")[, 1:2]

# 组装数据框
df <- data.frame(
  UMAP1 = coord[, 1],
  UMAP2 = coord[, 2],
  expr = FetchData(seurat_obj_all, vars = "evm.TU.PaChr24.43")[, 1],
  group = seurat_obj_all@meta.data$group
)

# 绘图
p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = expr)) +
  geom_point(size = 0.1) +
  facet_wrap(~group) +
  scale_color_viridis_c(option = "magma") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "PaChr24B.43"   # 修改图例标题
  )

# 显示图
print(p)

# 保存为 JPG 文件
ggsave(
  filename = "PaChr24B_43_by_group_harmony_clusters_0.8.jpg",
  plot = p,
  width = 10,
  height = 4,
  dpi = 300
)
