##############use SCT v2, and remember when subsetting the cells, use recorrect_umi=FALSE
setwd("/scratch/Cui3/singlecellanalysis")
seurat_obj_all<-readRDS("seurat_obj_all_harmony_res0.2_0.8.rds")
EU <- c("EU60","EU78","EU620")
seurat_obj_all$group <- sapply(seurat_obj_all$orig.ident, function(ita) ifelse(ita %in% EU,"EU","invasive"))
options(future.globals.maxSize = 10000 * 1024^2)  #

Idents(seurat_obj_all) <- seurat_obj_all$harmony_clusters_0.8
clusters <- sort(unique(Idents(seurat_obj_all)))
#dir.create("cluster_umap_harmony0.8", showWarnings = FALSE)
Idents(seurat_obj_all) <- seurat_obj_all$harmony_clusters_0.8
clusters <- sort(unique(Idents(seurat_obj_all)))
#dir.create("cluster_umap_harmony0.8", showWarnings = FALSE)

for (cl in clusters) {
  cat("Plotting cluster:", cl, "\n")
  
  p <- DimPlot(
    seurat_obj_all,
    reduction = "harmony_clusters_0.8",   
    cells.highlight = WhichCells(seurat_obj_all, idents = cl),
    cols.highlight = "red",
    cols = "gray80",
    pt.size = 0.8
  ) +
    ggtitle(paste("Cluster", cl)) +        
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),  
      panel.grid = element_blank(),                       
      legend.position = "none",                           
      axis.title = element_text(size = 16),              
      axis.text = element_text(size = 16)               
    )
  
  ggsave(
    filename = paste0("cluster_umap_harmony0.8/Cluster_", cl, "_UMAP.jpg"),  
    plot = p,
    width = 6,
    height = 5
  )
}
