library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
##载入脚本
source("CreateBmkObject.R")
spatial<- CreateS1000Object(matrix_path="./cell_split", png_path="./he_cell_split_small.png", spot_radius =NULL, min.cells =5, min.features = 100)
spatial <- SCTransform(spatial, assay = "Spatial", new.assay.name = "SCT", verbose = FALSE, vst.flavor = "v2")
#spatial<-readRDS("object.RDS")
spatial@meta.data$imagerow <- spatial@images$sample1@coordinates$imagerow
spatial@meta.data$imagecol <- spatial@images$sample1@coordinates$imagecol
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:35)
spatial <- FindClusters(spatial, resolution = 0.7, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:35)
saveRDS(spatial, "spatial_cellsplit_res0.7.rds")


Large <- subset(spatial, imagecol > 9845 )

pdf("visual.pdf")
SpatialDimPlot(spatial, crop = FALSE, label = TRUE, label.size = 3, pt.size.factor = 7)
dev.off()