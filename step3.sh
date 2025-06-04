Idents(seurat_obj_all) <- seurat_obj_all@meta.data$harmony_clusters_0.8
seurat_obj_all <- PrepSCTFindMarkers(seurat_obj_all)
all_markers <- FindAllMarkers(seurat_obj_all, logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE, test.use = "roc")
# Step 2: Rank genes within each cluster
ranked_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(power)) %>%
  ungroup()

# Step 3: Filter top genes based on avg_log2FC and then min.pct (pct.1)
top_markers <- ranked_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%  # Sort by avg_log2FC first, then pct.1
  slice_head(n = 2) %>%                      # Keep top 10 markers per cluster
  ungroup() 

# Extract unique genes from top_markers
unique_genes <- unique(top_markers$gene)

# Check for duplicates
print(length(top_markers$gene))  # Original length
print(length(unique_genes))
# Save the bubble plot to PDF
pdf("marker_bubble_plot.pdf", width = 12, height = 8)