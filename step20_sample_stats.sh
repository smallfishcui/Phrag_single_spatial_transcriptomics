cd /scratch/Cui3/spatial/spatial
# 
sample_obj_all<-readRDS("/scratch/Cui3/singlecellanalysis/seurat_obj_all_harmony_res0.2_0.8.rds")
sample_obj_list <- SplitObject(sample_obj_all, split.by = "orig.ident")
sample_obj <- sample_obj_list[[4]]
# 提取 seurat_obj_all.list 中的第 4 个样本
sample_name <- names(sample_obj_list)[4]  # 第 4 个样本的名称

sample_log10_genes <- log10(sample_obj@meta.data$nFeature_RNA + 1)  # log10(nFeature_RNA + 1)
sample_data <- data.frame(Sample = sample_name, Log10_Gene_Number = sample_log10_genes)
#spatial
spatial<-readRDS("spatial_res0.5_L4.rds")
spatial_log10_genes <- log10(spatial@meta.data$nFeature_Spatial + 1)  # log10(nFeature_Spatial + 1)
spatial_data <- data.frame(Sample = "Spatial", Log10_Gene_Number = spatial_log10_genes)
log10_gene_data <- bind_rows(sample_data, spatial_data)
head(log10_gene_data)
vln_plot <- ggplot(log10_gene_data, aes(x = Sample, y = Log10_Gene_Number, fill = Sample)) +
  geom_violin(scale = "width", trim = TRUE) + 
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +  
  scale_fill_viridis_d() + 
  theme_minimal() +  
  labs(title = "Log10(Feature Counts) per Sample",
       x = "Sample",
       y = "Log10(Number of Genes + 1)") +
  theme(legend.position = "none",  
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black") )  

pdf("log10_feature_counts_violinplot_sample4_spatial_L4.pdf", width = 8, height = 6)
print(vln_plot)
dev.off()


sample_cell_number <- ncol(sample_obj) 
sample_data <- data.frame(Sample = sample_name, Cell_Number = sample_cell_number)

spatial_cell_number <- ncol(spatial)  
spatial_data <- data.frame(Sample = "Spatial", Cell_Number = spatial_cell_number)

cell_number_data <- bind_rows(sample_data, spatial_data)

print(cell_number_data)

bar_plot <- ggplot(cell_number_data, aes(x = Sample, y = Cell_Number, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6) +  
  scale_fill_viridis_d() +  
  theme_minimal() +  
  labs(title = "Cell Number per Sample",
       x = "Sample",
       y = "Number of Cells") +
  theme(
    legend.position = "none",  
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black")  
  )

print(bar_plot)
pdf("cell_number_barplot_sample4_spatial.pdf", width = 6, height = 6)
print(bar_plot)
dev.off()
