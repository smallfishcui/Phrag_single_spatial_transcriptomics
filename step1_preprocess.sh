library(Seurat)
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(dplyr)
library(stringr)
library(clustree)
library(SeuratDisk)
library(SeuratData)
library(harmony)
library(remotes)
library(DoubletFinder)
library(sctransform)

bp.params <- MulticoreParam(workers = 8)
EU60.path <- "EU60/outs/filtered_feature_bc_matrix/"
EU78.path <- "EU78/outs/filtered_feature_bc_matrix/"
EU620.path <- "EU620/outs/filtered_feature_bc_matrix/"
NAint61.path <- "NAint61/outs/filtered_feature_bc_matrix/"
NAint113.path <- "NAint113/outs/filtered_feature_bc_matrix/"
NAint191.path <- "NAint191/outs/filtered_feature_bc_matrix/"
sce.EU60.sing <- read10xCounts(EU60.path, col.names=TRUE, BPPARAM = bp.params)
sce.EU78.sing <- read10xCounts(EU78.path, col.names=TRUE, BPPARAM = bp.params)
sce.EU620.sing <- read10xCounts(EU620.path, col.names=TRUE, BPPARAM = bp.params)
sce.NAint61.sing <- read10xCounts(NAint61.path, col.names=TRUE, BPPARAM = bp.params)
sce.NAint113.sing <- read10xCounts(NAint113.path, col.names=TRUE, BPPARAM = bp.params)
sce.NAint191.sing <- read10xCounts(NAint191.path, col.names=TRUE, BPPARAM = bp.params)
genesPerCellEU60 <- colSums(counts(sce.EU60.sing) > 0)
genesPerCellEU78 <- colSums(counts(sce.EU78.sing) > 0)
genesPerCellEU620 <- colSums(counts(sce.EU620.sing) > 0)
genesPerCellNAint61 <- colSums(counts(sce.NAint61.sing) > 0)
genesPerCellNAint113 <- colSums(counts(sce.NAint113.sing) > 0)
genesPerCellNAint191 <- colSums(counts(sce.NAint191.sing) > 0)
detected_genes_EU60 <- rowSums(counts(sce.EU60.sing)) > 0
detected_genes_EU78 <- rowSums(counts(sce.EU78.sing)) > 0
detected_genes_EU620 <- rowSums(counts(sce.EU620.sing)) > 0
detected_genes_NAint61 <- rowSums(counts(sce.NAint61.sing)) > 0
detected_genes_NAint113 <- rowSums(counts(sce.NAint113.sing)) > 0
detected_genes_NAint191 <- rowSums(counts(sce.NAint191.sing)) > 0
sce.EU60.sing <- sce.EU60.sing[detected_genes_EU60,]
sce.EU78.sing <- sce.EU78.sing[detected_genes_EU78,]
sce.EU620.sing <- sce.EU620.sing[detected_genes_EU620,]
sce.NAint61.sing <- sce.NAint61.sing[detected_genes_NAint61,]
sce.NAint113.sing <- sce.NAint113.sing[detected_genes_NAint113,]
sce.NAint191.sing <- sce.NAint191.sing[detected_genes_NAint191,]
rowData(sce.EU60.sing)
rowData(sce.EU78.sing)
rowData(sce.EU620.sing)
rowData(sce.NAint61.sing)
rowData(sce.NAint113.sing)
rowData(sce.NAint191.sing)
rowData(sce.EU60.sing)$Chromosome <- sapply(strsplit(rownames(sce.EU60.sing), "\\."), `[`, 3)
rowData(sce.EU78.sing)$Chromosome <- sapply(strsplit(rownames(sce.EU78.sing), "\\."), `[`, 3)
rowData(sce.EU620.sing)$Chromosome <- sapply(strsplit(rownames(sce.EU620.sing), "\\."), `[`, 3)
rowData(sce.NAint61.sing)$Chromosome <- sapply(strsplit(rownames(sce.NAint61.sing), "\\."), `[`, 3)
rowData(sce.NAint113.sing)$Chromosome <- sapply(strsplit(rownames(sce.NAint113.sing), "\\."), `[`, 3)
rowData(sce.NAint191.sing)$Chromosome <- sapply(strsplit(rownames(sce.NAint191.sing), "\\."), `[`, 3)
scaffolds_chloro <- c("HiC_scaffold_411", "HiC_scaffold_293", "HiC_scaffold_220","HiC_scaffold_217","HiC_scaffold_183","HiC_scaffold_622","HiC_scaffold_1055","HiC_scaffold_197","HiC_scaffold_395","HiC_scaffold_382","HiC_scaffold_549","HiC_scaffold_441","HiC_scaffold_284","HiC_scaffold_919","HiC_scaffold_138")
rowData(sce.EU60.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU60.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
rowData(sce.EU78.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU78.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
rowData(sce.EU620.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU620.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
rowData(sce.NAint61.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint61.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
rowData(sce.NAint113.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint113.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
rowData(sce.NAint191.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint191.sing)$Chromosome, 
  setNames(rep("PaChloro", length(scaffolds_chloro)), scaffolds_chloro)
)
scaffolds_mito <- c("HiC_scaffold_27","HiC_scaffold_1247","HiC_scaffold_277","HiC_scaffold_597","HiC_scaffold_507","HiC_scaffold_30","HiC_scaffold_29","HiC_scaffold_1312","HiC_scaffold_135","HiC_scaffold_511","HiC_scaffold_50")
rowData(sce.EU60.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU60.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
rowData(sce.EU78.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU78.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
rowData(sce.EU620.sing)$Chromosome <- str_replace_all(
  rowData(sce.EU620.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
rowData(sce.NAint61.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint61.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
rowData(sce.NAint113.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint113.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
rowData(sce.NAint191.sing)$Chromosome <- str_replace_all(
  rowData(sce.NAint191.sing)$Chromosome, 
  setNames(rep("PaMito", length(scaffolds_mito)), scaffolds_mito)
)
is.mito_EU60 <- which(rowData(sce.EU60.sing)$Chromosome=="PaMito")
is.mito_EU78 <- which(rowData(sce.EU78.sing)$Chromosome=="PaMito")
is.mito_EU620 <- which(rowData(sce.EU620.sing)$Chromosome=="PaMito")
is.mito_NAint61 <- which(rowData(sce.NAint61.sing)$Chromosome=="PaMito")
is.mito_NAint113 <- which(rowData(sce.NAint113.sing)$Chromosome=="PaMito")
is.mito_NAint191 <- which(rowData(sce.NAint191.sing)$Chromosome=="PaMito")
is.chloro_EU60 <-which(rowData(sce.EU60.sing)$Chromosome=="PaChloro")
is.chloro_EU78 <-which(rowData(sce.EU78.sing)$Chromosome=="PaChloro")
is.chloro_EU620 <-which(rowData(sce.EU620.sing)$Chromosome=="PaChloro")
is.chloro_NAint61 <-which(rowData(sce.NAint61.sing)$Chromosome=="PaChloro")
is.chloro_NAint113 <-which(rowData(sce.NAint113.sing)$Chromosome=="PaChloro")
is.chloro_NAint191 <-which(rowData(sce.NAint191.sing)$Chromosome=="PaChloro")
sce.EU60.sing <- addPerCellQC(sce.EU60.sing, subsets=list(Mito=is.mito_EU60),BPPARAM = bp.params)
sce.EU78.sing <- addPerCellQC(sce.EU78.sing, subsets=list(Mito=is.mito_EU78),BPPARAM = bp.params)
sce.EU620.sing <- addPerCellQC(sce.EU620.sing, subsets=list(Mito=is.mito_EU620),BPPARAM = bp.params)
sce.NAint61.sing <- addPerCellQC(sce.NAint61.sing, subsets=list(Mito=is.mito_NAint61),BPPARAM = bp.params)
sce.NAint113.sing <- addPerCellQC(sce.NAint113.sing, subsets=list(Mito=is.mito_NAint113),BPPARAM = bp.params)
sce.NAint191.sing <- addPerCellQC(sce.NAint191.sing, subsets=list(Mito=is.mito_NAint191),BPPARAM = bp.params)
#colData(sce.EU60.sing)
sce.EU60.sing <- addPerCellQC(sce.EU60.sing, subsets=list(Chloro=is.chloro_EU60),BPPARAM = bp.params)
sce.EU78.sing <- addPerCellQC(sce.EU78.sing, subsets=list(Chloro=is.chloro_EU78),BPPARAM = bp.params)
sce.EU620.sing <- addPerCellQC(sce.EU620.sing, subsets=list(Chloro=is.chloro_EU620),BPPARAM = bp.params)
sce.NAint61.sing <- addPerCellQC(sce.NAint61.sing, subsets=list(Chloro=is.chloro_NAint61),BPPARAM = bp.params)
sce.NAint113.sing <- addPerCellQC(sce.NAint113.sing, subsets=list(Chloro=is.chloro_NAint113),BPPARAM = bp.params)
sce.NAint191.sing <- addPerCellQC(sce.NAint191.sing, subsets=list(Chloro=is.chloro_NAint191),BPPARAM = bp.params)
colData(sce.EU60.sing)
colData(sce.EU78.sing)
colData(sce.EU620.sing)
colData(sce.NAint61.sing)
colData(sce.NAint113.sing)
colData(sce.NAint191.sing)
low_lib_size_EU60 <- isOutlier(sce.EU60.sing$sum, log=TRUE, type="lower")
low_lib_size_EU78 <- isOutlier(sce.EU78.sing$sum, log=TRUE, type="lower")
low_lib_size_EU620 <- isOutlier(sce.EU620.sing$sum, log=TRUE, type="lower")
low_lib_size_NAint61 <- isOutlier(sce.NAint61.sing$sum, log=TRUE, type="lower")
low_lib_size_NAint113 <- isOutlier(sce.NAint113.sing$sum, log=TRUE, type="lower")
low_lib_size_NAint191 <- isOutlier(sce.NAint191.sing$sum, log=TRUE, type="lower")
table(low_lib_size_EU60)
#low_lib_size_EU60
#FALSE  TRUE 
#11427   278 
table(low_lib_size_EU78)
#low_lib_size_EU78
#FALSE  TRUE 
#11854   371 
table(low_lib_size_EU620)
#low_lib_size_EU620
#FALSE 
#7154 
table(low_lib_size_NAint61)
#low_lib_size_NAint61
#FALSE 
# 8645 
table(low_lib_size_NAint113)
#low_lib_size_NAint113
#FALSE 
#15297 
table(low_lib_size_NAint191)
#low_lib_size_NAint191
#FALSE 
#12168 
attr(low_lib_size_EU60, "thresholds")
#  lower   higher 
#635.2608      Inf 
attr(low_lib_size_EU78, "thresholds")
#  lower   higher 
#733.2962      Inf 
attr(low_lib_size_EU620, "thresholds")
#  lower  higher 
#302.955     Inf 
attr(low_lib_size_NAint61, "thresholds")
#   lower   higher 
#67.01136      Inf 
attr(low_lib_size_NAint113, "thresholds")
#   lower   higher 
#407.0389      Inf 
attr(low_lib_size_NAint191, "thresholds")
#   lower   higher 
#379.1528      Inf 
low_n_features_EU60 <- isOutlier(sce.EU60.sing$detected, log=TRUE, type="lower")
low_n_features_EU78 <- isOutlier(sce.EU78.sing$detected, log=TRUE, type="lower")
low_n_features_EU620 <- isOutlier(sce.EU620.sing$detected, log=TRUE, type="lower")
low_n_features_NAint61 <- isOutlier(sce.NAint61.sing$detected, log=TRUE, type="lower")
low_n_features_NAint113 <- isOutlier(sce.NAint113.sing$detected, log=TRUE, type="lower")
low_n_features_NAint191 <- isOutlier(sce.NAint191.sing$detected, log=TRUE, type="lower")
table(low_n_features_EU60)
#low_n_features_EU60
#FALSE  TRUE 
#11292   413 
table(low_n_features_EU60)
#FALSE  TRUE 
#11292   413 
table(low_n_features_EU78)
#FALSE  TRUE 
#11787   438 
table(low_n_features_EU620)
#FALSE  TRUE 
# 7137    17 
table(low_n_features_NAint61)
#FALSE 
# 8645 
table(low_n_features_NAint113)
#FALSE  TRUE 
#15296     1 
table(low_n_features_NAint191)
#FALSE  TRUE 
#12064   104 
attr(low_n_features_EU60, "thresholds")[1]
attr(low_n_features_EU78, "thresholds")[1]
attr(low_n_features_EU620, "thresholds")[1]
attr(low_n_features_NAint61, "thresholds")[1]
attr(low_n_features_NAint113, "thresholds")[1]
attr(low_n_features_NAint191, "thresholds")[1]
high_Mito_percent_EU60 <- isOutlier(sce.EU60.sing$subsets_Mito_percent, type="higher")
high_Mito_percent_EU78 <- isOutlier(sce.EU78.sing$subsets_Mito_percent, type="higher")
high_Mito_percent_EU620 <- isOutlier(sce.EU620.sing$subsets_Mito_percent, type="higher")
high_Mito_percent_NAint61 <- isOutlier(sce.NAint61.sing$subsets_Mito_percent, type="higher")
high_Mito_percent_NAint113 <- isOutlier(sce.NAint113.sing$subsets_Mito_percent, type="higher")
high_Mito_percent_NAint191 <- isOutlier(sce.NAint191.sing$subsets_Mito_percent, type="higher")
table(high_Mito_percent_EU60)
#FALSE  TRUE 
#10531  1174 
table(high_Mito_percent_EU78)
#FALSE  TRUE 
#10969  1256 
table(high_Mito_percent_EU620)
#FALSE  TRUE 
# 6246   908 
table(high_Mito_percent_NAint61)
#FALSE  TRUE 
# 7949   696 
table(high_Mito_percent_NAint113)
#FALSE  TRUE 
#13993  1304
table(high_Mito_percent_NAint191)
#FALSE  TRUE 
#11003  1165 
attr(high_Mito_percent_EU60, "thresholds")[2]
#   higher 
#0.8331792 
attr(high_Mito_percent_EU78, "thresholds")[2]
#  higher 
#0.8711276 
attr(high_Mito_percent_EU620, "thresholds")[2]
# higher 
#1.602231 
attr(high_Mito_percent_NAint61, "thresholds")[2]
#  higher 
#1.608463 
attr(high_Mito_percent_NAint113, "thresholds")[2]
#  higher 
#1.143472 
attr(high_Mito_percent_NAint191, "thresholds")[2]
#  higher 
#1.392472 
high_Chloro_percent_EU60 <- isOutlier(sce.EU60.sing$subsets_Chloro_percent, type="higher")
high_Chloro_percent_EU78 <- isOutlier(sce.EU78.sing$subsets_Chloro_percent, type="higher")
high_Chloro_percent_EU620 <- isOutlier(sce.EU620.sing$subsets_Chloro_percent, type="higher")
high_Chloro_percent_NAint61 <- isOutlier(sce.NAint61.sing$subsets_Chloro_percent, type="higher")
high_Chloro_percent_NAint113 <- isOutlier(sce.NAint113.sing$subsets_Chloro_percent, type="higher")
high_Chloro_percent_NAint191 <- isOutlier(sce.NAint191.sing$subsets_Chloro_percent, type="higher")
table(high_Chloro_percent_EU60)
#FALSE  TRUE 
# 9309  2396 
table(high_Chloro_percent_EU78)
#FALSE  TRUE 
# 9716  2509 
table(high_Chloro_percent_EU620)
#FALSE  TRUE 
# 5714  1440 
table(high_Chloro_percent_NAint61)
#FALSE  TRUE 
# 7158  1487 
table(high_Chloro_percent_NAint113)
#FALSE  TRUE 
#12964  2333 
table(high_Chloro_percent_NAint191)
#FALSE  TRUE 
#10331  1837
attr(high_Chloro_percent_EU60, "thresholds")[2]
#   higher 
#0.2472901 
attr(high_Chloro_percent_EU78, "thresholds")[2]
#   higher 
#0.1986436 
attr(high_Chloro_percent_EU620, "thresholds")[2]
#   higher 
#0.7996841 
attr(high_Chloro_percent_NAint61, "thresholds")[2]
#   higher 
#0.5717528 
attr(high_Chloro_percent_NAint113, "thresholds")[2]
#   higher 
#0.2510121 
attr(high_Chloro_percent_NAint191, "thresholds")[2]
#  higher 
#1.109069 
sce.EU60.sing <- addPerFeatureQC(sce.EU60.sing, BPPARAM = bp.params)
sce.EU78.sing <- addPerFeatureQC(sce.EU78.sing, BPPARAM = bp.params)
sce.EU620.sing <- addPerFeatureQC(sce.EU620.sing, BPPARAM = bp.params)
sce.NAint61.sing <- addPerFeatureQC(sce.NAint61.sing, BPPARAM = bp.params)
sce.NAint113.sing <- addPerFeatureQC(sce.NAint113.sing, BPPARAM = bp.params)
sce.NAint191.sing <- addPerFeatureQC(sce.NAint191.sing, BPPARAM = bp.params)
rowData(sce.EU60.sing)
rowData(sce.EU78.sing)
rowData(sce.EU620.sing)
rowData(sce.NAint61.sing)
rowData(sce.NAint113.sing)
rowData(sce.NAint191.sing)
colData(sce.EU60.sing)$cell_sparsity <- 1 - (colData(sce.EU60.sing)$detected / nrow(sce.EU60.sing))
colData(sce.EU78.sing)$cell_sparsity <- 1 - (colData(sce.EU78.sing)$detected / nrow(sce.EU78.sing))
colData(sce.EU620.sing)$cell_sparsity <- 1 - (colData(sce.EU620.sing)$detected / nrow(sce.EU620.sing))
colData(sce.NAint61.sing)$cell_sparsity <- 1 - (colData(sce.NAint61.sing)$detected / nrow(sce.NAint61.sing))
colData(sce.NAint113.sing)$cell_sparsity <- 1 - (colData(sce.NAint113.sing)$detected / nrow(sce.NAint113.sing))
colData(sce.NAint191.sing)$cell_sparsity <- 1 - (colData(sce.NAint191.sing)$detected / nrow(sce.NAint191.sing))
rowData(sce.EU60.sing)$gene_sparsity <- (100 - rowData(sce.EU60.sing)$detected) / 100
rowData(sce.EU78.sing)$gene_sparsity <- (100 - rowData(sce.EU78.sing)$detected) / 100
rowData(sce.EU620.sing)$gene_sparsity <- (100 - rowData(sce.EU620.sing)$detected) / 100
rowData(sce.NAint61.sing)$gene_sparsity <- (100 - rowData(sce.NAint61.sing)$detected) / 100
rowData(sce.NAint113.sing)$gene_sparsity <- (100 - rowData(sce.NAint113.sing)$detected) / 100
rowData(sce.NAint191.sing)$gene_sparsity <- (100 - rowData(sce.NAint191.sing)$detected) / 100