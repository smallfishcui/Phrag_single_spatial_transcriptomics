cd /scratch/Cui3/spatial/spatial
####!becareful! Remove the old file before running this script. otherwise results will apprehend
for i in {0..16}; do
  cat "spatial_markers_res0.5_cluster_arabdop${i}.gene" >> merged.spatial_markers_res0.5_cluster_arabdopL4.gene
done
#convert to symbol, spatial
library(org.At.tair.db)
library(AnnotationDbi)
genes_orig <- read.csv("merged.spatial_markers_res0.5_cluster_arabdopL4.gene", sep=" ", header=F)
keytypes(org.At.tair.db)
genes <- as.character(genes_orig$V2)
genes <- gsub("\\.\\d+$", "", genes)
mapped_symbol<-mapIds(org.At.tair.db, keys = genes,
  column = c('SYMBOL'), keytype = 'TAIR')
# Create a data frame with original V1, V2, and mapped symbols
# Create a data frame combining original data and mapped symbols
output <- data.frame(
  Original_ID = genes_orig$V1,  # Original gene ID (V1)
  TAIR_ID = genes_orig$V2,      # Unprocessed TAIR ID (V2)
  Processed_TAIR_ID = genes,    # Processed TAIR ID
  SYMBOL = mapped_symbol[genes]  # Mapped SYMBOLs
)

# Write the result to a CSV file
write.csv(output, file = "merged.spatial_markers_res0.5_cluster_arabdop.gene_with_symbolsL4.csv", row.names = FALSE, quote = FALSE)

# Load necessary libraries
library(dplyr)
# Read the first file
file1 <- read.csv("spatial_markers_res0.5_L4.csv", stringsAsFactors = FALSE)
# 获取第一列的列名
first_col <- colnames(file1)[1]

# 把第一列中所有包含 "TU" 的基因 ID 中的 "TU" 替换为 "model"
file1[[first_col]] <- gsub("TU", "model", file1[[first_col]])

# Read the second file
file2 <- read.csv("merged.spatial_markers_res0.5_cluster_arabdop.gene_with_symbolsL4.csv", stringsAsFactors = FALSE)
# Rename the first column in file2 to match the first column in file1
colnames(file2)[1] <- colnames(file1)[1]
# Merge the two files based on the first column
merged_data <- merge(file1, file2, by = colnames(file1)[1], all.x = TRUE)
# Save the merged data to a new CSV file
write.csv(merged_data, file = "result_spatial_markers_res0.5_cluster_arabdop.gene_with_symbolL4.csv", row.names = FALSE, quote = FALSE)
#################GO
