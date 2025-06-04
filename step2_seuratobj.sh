seurat_obj_EU60 <- CreateSeuratObject(
  counts = counts(sce.EU60.sing),
  meta.data = as.data.frame(colData(sce.EU60.sing)), project = "EU60", min.cells = 20, min.features = 500
)
seurat_obj_EU78 <- CreateSeuratObject(
  counts = counts(sce.EU78.sing),
  meta.data = as.data.frame(colData(sce.EU78.sing)), project = "EU78", min.cells = 20, min.features = 500
)
seurat_obj_EU620 <- CreateSeuratObject(
  counts = counts(sce.EU620.sing),
  meta.data = as.data.frame(colData(sce.EU620.sing)), project = "EU620", min.cells = 20, min.features = 500
)
seurat_obj_NAint61 <- CreateSeuratObject(
  counts = counts(sce.NAint61.sing),
  meta.data = as.data.frame(colData(sce.NAint61.sing)), project = "NAint61", min.cells = 20, min.features = 500
)
seurat_obj_NAint113 <- CreateSeuratObject(
  counts = counts(sce.NAint113.sing),
  meta.data = as.data.frame(colData(sce.NAint113.sing)), project = "NAint113", min.cells = 20, min.features = 500
)
seurat_obj_NAint191 <- CreateSeuratObject(
  counts = counts(sce.NAint191.sing),
  meta.data = as.data.frame(colData(sce.NAint191.sing)), project = "NAint191", min.cells = 20, min.features = 500
)
#filter data, chloroplast and mitochondira less than 30%, 
seurat_obj_EU60 <- subset(seurat_obj_EU60, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
seurat_obj_EU78 <- subset(seurat_obj_EU78, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
seurat_obj_EU620 <- subset(seurat_obj_EU620, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
seurat_obj_NAint61 <- subset(seurat_obj_NAint61, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
seurat_obj_NAint113 <- subset(seurat_obj_NAint113, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
seurat_obj_NAint191 <- subset(seurat_obj_NAint191, subset = nFeature_RNA > 10 & nFeature_RNA < 15000 & subsets_Mito_percent < 10 & subsets_Chloro_percent < 30)
#seurat_obj_EU60<- SCTransform( seurat_obj_EU60, assay = "RNA", new.assay.name = "SCT")
seurat_obj_EU60 <- SCTransform(seurat_obj_EU60, vst.flavor = "v2")
seurat_obj_EU78 <- SCTransform(seurat_obj_EU78, vst.flavor = "v2")
seurat_obj_EU620 <- SCTransform(seurat_obj_EU620, vst.flavor = "v2")
seurat_obj_NAint61 <- SCTransform(seurat_obj_NAint61, vst.flavor = "v2")
seurat_obj_NAint113 <- SCTransform(seurat_obj_NAint113, vst.flavor = "v2")
seurat_obj_NAint191 <- SCTransform(seurat_obj_NAint191, vst.flavor = "v2")
seurat_obj_EU60 <- ScaleData(seurat_obj_EU60, verbose = FALSE)
seurat_obj_EU78 <- ScaleData(seurat_obj_EU78, verbose = FALSE)
seurat_obj_EU620 <- ScaleData(seurat_obj_EU620, verbose = FALSE)
seurat_obj_NAint61 <- ScaleData(seurat_obj_NAint61, verbose = FALSE)
seurat_obj_NAint113 <- ScaleData(seurat_obj_NAint113, verbose = FALSE)
seurat_obj_NAint191 <- ScaleData(seurat_obj_NAint191, verbose = FALSE)
seurat_obj_EU60<-RunPCA(seurat_obj_EU60)
seurat_obj_EU78<-RunPCA(seurat_obj_EU78)
seurat_obj_EU620<-RunPCA(seurat_obj_EU620)
seurat_obj_NAint61<-RunPCA(seurat_obj_NAint61)
seurat_obj_NAint113<-RunPCA(seurat_obj_NAint113)
seurat_obj_NAint191<-RunPCA(seurat_obj_NAint191)

sweep.res.list <- paramSweep(seurat_obj_EU60, PCs = 1:25, sct = TRUE)# pN = 0.25, pK = 0.05)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_EU60.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.04  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_EU60))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_EU60 <- doubletFinder(seurat_obj_EU60, PCs = 1:25, pN = 0.25,         # Proportion of artificial doublets
  pK = 0.04,  nExp = nExp,      
  reuse.pANN = FALSE,
  sct = TRUE
)
seurat_obj_EU60 <- subset(seurat_obj_EU60, subset = DF.classifications_0.25_0.04_456 == "Singlet")
seurat_obj_EU60 <- FindNeighbors(seurat_obj_EU60, reduction = "pca", dims = 1:35)
seurat_obj_EU60 <- FindClusters(seurat_obj_EU60, resolution = 0.3, algorithm = 2, random.seed = 100)
seurat_obj_EU60 <- RunUMAP(seurat_obj_EU60, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.EU60_no_doublets.pdf")
DimPlot(seurat_obj_EU60, reduction = "umap", group.by = "seurat_clusters")
dev.off()

sweep.res.list <- paramSweep(seurat_obj_EU78, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_EU78.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.02  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_EU78))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_EU78 <- doubletFinder(
  seurat_obj_EU78,
  PCs = 1:25,        # Number of PCs used in clustering
  pN = 0.25,         # Proportion of artificial doublets
  pK = 0.04,   # Optimal pK found
  nExp = nExp,       # Expected number of doublets
  reuse.pANN = FALSE,
  sct = TRUE
)

# Step 5: Filter out doublets
seurat_obj_EU78 <- subset(seurat_obj_EU78, subset = DF.classifications_0.25_0.04_239 == "Singlet")
seurat_obj_EU78 <- FindNeighbors(seurat_obj_EU78, reduction = "pca", dims = 1:35)
seurat_obj_EU78 <- FindClusters(seurat_obj_EU78, resolution = 0.2, algorithm = 2, random.seed = 100)
seurat_obj_EU78 <- RunUMAP(seurat_obj_EU78, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.EU78_no_doublets.pdf")
DimPlot(seurat_obj_EU78, reduction = "umap", group.by = "seurat_clusters")
dev.off()

sweep.res.list <- paramSweep(seurat_obj_EU620, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_EU620.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.06  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_EU620))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_EU620 <- doubletFinder(
  seurat_obj_EU620,
  PCs = 1:25,        # Number of PCs used in clustering
  pN = 0.25,         # Proportion of artificial doublets
  pK = 0.06,   # Optimal pK found
  nExp = nExp,       # Expected number of doublets
  reuse.pANN = FALSE,
  sct = TRUE
)

# Step 5: Filter out doublets
seurat_obj_EU620 <- subset(seurat_obj_EU620, subset = DF.classifications_0.25_0.06_412 == "Singlet")
seurat_obj_EU620 <- FindNeighbors(seurat_obj_EU620, reduction = "pca", dims = 1:35)
seurat_obj_EU620 <- FindClusters(seurat_obj_EU620, resolution = 0.2, algorithm = 2, random.seed = 100)
seurat_obj_EU620 <- RunUMAP(seurat_obj_EU620, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.EU620_no_doublets.pdf")
DimPlot(seurat_obj_EU620, reduction = "umap", group.by = "seurat_clusters")
dev.off()

sweep.res.list <- paramSweep(seurat_obj_NAint61, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_NAint61.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.14  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_NAint61))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_NAint61 <- doubletFinder(
  seurat_obj_NAint61,
  PCs = 1:25,        # Number of PCs used in clustering
  pN = 0.25,         # Proportion of artificial doublets
  pK = 0.14,   # Optimal pK found
  nExp = nExp,       # Expected number of doublets
  reuse.pANN = FALSE,
  sct = TRUE
)

# Step 5: Filter out doublets
seurat_obj_NAint61 <- subset(seurat_obj_NAint61, subset = DF.classifications_0.25_0.14_1107 == "Singlet")
seurat_obj_NAint61 <- FindNeighbors(seurat_obj_NAint61, reduction = "pca", dims = 1:35)
seurat_obj_NAint61 <- FindClusters(seurat_obj_NAint61, resolution = 0.3, algorithm = 2, random.seed = 100)
seurat_obj_NAint61 <- RunUMAP(seurat_obj_NAint61, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.NAint61_no_doublets.pdf")
DimPlot(seurat_obj_NAint61, reduction = "umap", group.by = "seurat_clusters")
dev.off()

sweep.res.list <- paramSweep(seurat_obj_NAint113, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_NAint113.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.11  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_NAint113))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_NAint113 <- doubletFinder(
  seurat_obj_NAint113,
  PCs = 1:25,        # Number of PCs used in clustering
  pN = 0.25,         # Proportion of artificial doublets
  pK = 0.11,   # Optimal pK found
  nExp = nExp,       # Expected number of doublets
  reuse.pANN = FALSE,
  sct = TRUE
)

# Step 5: Filter out doublets
seurat_obj_NAint113 <- subset(seurat_obj_NAint113, subset = DF.classifications_0.25_0.11_1554 == "Singlet")
seurat_obj_NAint113 <- FindNeighbors(seurat_obj_NAint113, reduction = "pca", dims = 1:35)
seurat_obj_NAint113 <- FindClusters(seurat_obj_NAint113, resolution = 0.3, algorithm = 2, random.seed = 100)
seurat_obj_NAint113 <- RunUMAP(seurat_obj_NAint113, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.NAint113_no_doublets.pdf")
DimPlot(seurat_obj_NAint113, reduction = "umap", group.by = "seurat_clusters")
dev.off()

sweep.res.list <- paramSweep(seurat_obj_NAint191, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
pdf("find_pK_plot_NAint191.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()  ##highest pk=0.03
# Step 2: Choose optimal pK
pK_optimal <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% pull(pK)
# Step 3: Set expected doublet rate and calculate nExp
doublet_rate <- 0.25  # Example rate of 3%
nExp <- round(doublet_rate * ncol(seurat_obj_NAint191))
# Step 4: Run DoubletFinder with chosen parameters
pK <- as.numeric(as.character(pK_optimal))
seurat_obj_NAint191 <- doubletFinder(
  seurat_obj_NAint191,
  PCs = 1:25,        # Number of PCs used in clustering
  pN = 0.25,         # Proportion of artificial doublets
  pK = 0.25,   # Optimal pK found
  nExp = nExp,       # Expected number of doublets
  reuse.pANN = FALSE,
  sct = TRUE
)

# Step 5: Filter out doublets
seurat_obj_NAint191 <- subset(seurat_obj_NAint191, subset = DF.classifications_0.25_0.25_2905 == "Singlet")
seurat_obj_NAint191 <- FindNeighbors(seurat_obj_NAint191, reduction = "pca", dims = 1:35)
seurat_obj_NAint191 <- FindClusters(seurat_obj_NAint191, resolution = 0.2, algorithm = 2, random.seed = 100)
seurat_obj_NAint191 <- RunUMAP(seurat_obj_NAint191, dims = 1:35, reduction = "pca", min.dist = 0.1, n.neighbors = 50)
pdf("umap.NAint191_no_doublets.pdf")
DimPlot(seurat_obj_NAint191, reduction = "umap", group.by = "seurat_clusters")
dev.off()

seurat_obj_all <- merge(seurat_obj_EU60, y = c(seurat_obj_EU78, seurat_obj_EU620, seurat_obj_NAint61, seurat_obj_NAint113, seurat_obj_NAint191), add.cell.ids = c("EU60", "EU78", "EU620", "NAint61", "NAint113", "NAint191"), project = "Phrag")
#head(colnames(seurat_obj_all))
#tail(colnames(seurat_obj_all))
#unique(sapply(X = strsplit(colnames(seurat_obj_all), split = "_"), FUN = "[", 1))
#table(seurat_obj_all$orig.ident)
seurat_obj_all$batch <- seurat_obj_all$orig.ident
seurat_obj_all <- SCTransform(seurat_obj_all, vst.flavor = "v2")
#2.使用 ScaleData
seurat_obj_all <- ScaleData(seurat_obj_all, verbose = FALSE)
#seurat_obj_all[["RNA"]] <- split(seurat_obj_all[["RNA"]], f = seurat_obj_all$orig.ident)
seurat_obj_all<-RunPCA(seurat_obj_all)

seurat_obj_all<- RunHarmony(seurat_obj_all, group.by.vars = "batch", verbose = TRUE)
seurat_obj_all <- FindNeighbors(seurat_obj_all, reduction = "harmony", dims = 1:35)
#seurat_obj_all <- IntegrateLayers(object = seurat_obj_all, method = HarmonyIntegration, orig.reduction = "pca",
#  new.reduction = 'harmony', assay = "SCT", verbose = FALSE)
# re-join layers after integration
#seurat_obj_all[["RNA"]] <- JoinLayers(seurat_obj_all[["RNA"]])
#seurat_obj_all <- FindNeighbors(seurat_obj_all, reduction = "harmony", dims = 1:35)
#seurat_obj_all <- FindNeighbors(seurat_obj_all, dims = 1:50, reduction = "pca")
# Perform clustering at a series of resolutions from 0 to 1
for (i in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)){
  seurat_obj_all <- FindClusters(seurat_obj_all, resolution = i, algorithm = 2, random.seed = 100)}

# Generate clustering tree
treeplot<-clustree(seurat_obj_all, prefix = "SCT_snn_res.")

cluster_sc3s<-as.data.frame(treeplot$data %>% group_by(SCT_snn_res.) %>% summarise(mean_sc3=mean(sc3_stability)))
most_stable_res<-subset(cluster_sc3s,mean_sc3==max(cluster_sc3s$mean_sc3))[,1]
print(paste("Most stable resolution determined to be:",most_stable_res))

# Make a plot of SC3 index values
sc3_plot<-ggplot(cluster_sc3s,aes(x=SCT_snn_res.,y=mean_sc3,group = 1))+
  theme_bw()+
  geom_line()+
  geom_point()+
  labs(x="Clustering resolution",y="SC3 stability index")
pdf("six_samples_harmony_sc3.pdf")
sc3_plot
dev.off()
#choose 0.2 as the best resolution, and then try 0.5
seurat_obj_all <- FindClusters(seurat_obj_all, resolution = 0.2, cluster.name = "harmony_clusters_0.2")
#seurat_obj_all <- FindClusters(seurat_obj_all, resolution = 0.8, cluster.name = "harmony_clusters_0.8")
seurat_obj_all <- RunUMAP(seurat_obj_all, dims = 1:35, reduction = "harmony", min.dist = 0.1, n.neighbors = 50, reduction.name = "harmony_clusters_0.2")
#seurat_obj_all <- RunUMAP(seurat_obj_all, dims = 1:35, reduction = "harmony", min.dist = 0.1, n.neighbors = 50, reduction.name = "harmony_clusters_0.8")
pdf("umap.harmony_sixsample_nodoublets_resolution0.2.pdf")
#DimPlot(seurat_obj_all, reduction = "harmony_clusters_0.2", group.by = "orig.ident")
DimPlot(seurat_obj_all, reduction = "harmony_clusters_0.2", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + scale_color_brewer(palette = "Set2") + labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")
dev.off()
pdf("umap.harmony_sixsample_nodoublets_resolution0.2.pdf")
DimPlot(seurat_obj_all, reduction = "harmony_clusters_0.2", group.by = "seurat_clusters") +
  scale_color_manual(values = c(
  "#E63946", "#F4A261", "#2A9D8F", "#264653",
  "#1D3557", "#A8DADC", "#457B9D", "#E76F51" # Add a 9th color
)) + # Adjust colors as needed
  theme(
    panel.background = element_rect(fill = "black", color = "black"),  # Black background
    plot.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "white"),  # Keep x and y axis lines in white
    axis.text = element_text(color = "white"),  # White axis text
    axis.title = element_text(color = "white"),  # White axis labels
    legend.background = element_rect(fill = "black", color = "black"),  # Black legend background
    legend.text = element_text(color = "white"),  # White legend text
    legend.title = element_text(color = "white")  # White legend title
  ) +
  labs(
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    title = "UMAP of Harmonized Clusters with Black Background"
  )
dev.off()


saveRDS(seurat_obj_all, file = "./seurat_obj_all_harmony_res0.2_0.8.rds")