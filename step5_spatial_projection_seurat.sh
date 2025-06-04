source("CreateBmkObject.R")
##创建对象
spatial<- CreateS1000Object(matrix_path="./L7_heAuto", png_path="./he_roi_small.png", spot_radius =NULL, min.cells =5, min.features = 100)
spatial <- SCTransform(spatial, assay = "Spatial", new.assay.name = "SCT", verbose = FALSE)
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:50)
spatial <- FindClusters(spatial, resolution = 0.5, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:50)
saveRDS(spatial,"spatial_res0.5.rds")
seurat_obj_all.list <- SplitObject(seurat_obj_all, split.by = "orig.ident")
###memory consuming, need 200g
seurat_obj_all.list <- lapply(X = seurat_obj_all.list, FUN = SCTransform, method = "glmGamPoi") # vars.to.regress not in example
anchors <- FindTransferAnchors(reference = spatial, query = seurat_obj_all.list[[1]], normalization.method = "SCT" )
summary(as.data.frame(anchors@anchors)$score)
anchors_data <- as.data.frame(anchors@anchors)
anchors_filtered <- anchors_data[anchors_data$score > 0.5, ]#filter score
anchors@anchors <- anchors_filtered
predictions.assay <- TransferData(anchorset = anchors, refdata = spatial$seurat_clusters, prediction.assay = TRUE,
    weight.reduction = seurat_obj_all.list[[1]][["pca"]], dims = 1:30)

spatial[["predictions"]] <- predictions.assay
# unique(seurat_obj_all_combined$seurat_clusters)
pdf("integrated_clusters_allsix_spatial_.pdf")
SpatialFeaturePlot(spatial, features = c("0", "1","2","3","4","5","6","7","8"), pt.size.factor = 7, ncol = 2, crop = TRUE)
dev.off()

anchors <- FindTransferAnchors(reference = seurat_obj_all.list[[1]], query = spatial, normalization.method = "SCT" )
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat_obj_all_combined$seurat_clusters, prediction.assay = TRUE,
    weight.reduction = spatial[["pca"]], dims = 1:30)
spatial[["predictions"]] <- predictions.assay
# unique(seurat_obj_all_combined$seurat_clusters)
pdf("integrated_clusters_allsix_spatial_.pdf")
SpatialFeaturePlot(spatial, features = c("0", "1","2","3","4","5","6","7","8"), pt.size.factor = 7, ncol = 2, crop = TRUE)
dev.off()

features <- SelectIntegrationFeatures(object.list = list(seurat_obj_all.list[[6]], spatial), nfeatures = 3000)
anchors <- FindTransferAnchors(reference = seurat_obj_all.list[[6]], query = spatial, normalization.method = "SCT" )
refdata <- seurat_obj_all.list[[6]]@meta.data$harmony_clusters_0.8
refdata <- as.vector(refdata)
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = refdata,
  prediction.assay = TRUE,
  weight.reduction = spatial[["pca"]],
  dims = 1:30
)
spatial[["predictions"]] <- predictions.assay
# unique(seurat_obj_all_combined$seurat_clusters)
pdf("integrated_clusters_allsix_spatial_nodoublets_0.8_6_part1.pdf")
SpatialFeaturePlot(spatial, features = c("0", "1","2","3"), pt.size.factor = 7, ncol = 2, crop = TRUE)+ theme(
  legend.text = element_text(size = 6), # Adjust text size
  legend.title = element_text(size = 8) # Adjust title size (if needed)
)
dev.off()