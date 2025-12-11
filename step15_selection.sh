export C_INCLUDE_PATH=/scratch/Cui3/Salt/VariantWGS_1/Selection/RAiSD-AI/RAiSD-AI-master/libpng-1.6.40/install/include
export LIBRARY_PATH=/scratch/Cui3/Salt/VariantWGS_1/Selection/RAiSD-AI/RAiSD-AI-master/libpng-1.6.40/install/lib
#test
ln -s /scratch/cuiwang/nodelete5/WGS/EUlineage/dv_9samples.biallelic.recode.vcf.gz .
/scratch/Cui3/Salt/VariantWGS_1/Selection/RAiSD-AI/RAiSD-AI-master/bin/release/RAiSD-AI -n dv_9samples.biallelic -I dv_9samples.biallelic.recode.vcf -w 100000
#sweepFinder2
#install via tykky
#conda-containerize new --mamba --prefix SweepFinder2 env.yml
cd /scratch/Cui3/Salt/VariantWGS_1/Selection/
export PATH="/projappl/SweepFinder2/bin:$PATH"

##this is non-sense, be careful when using it, especially the snp sites, modified by chr number
for i in {1..25}; do
    sed '1d' ./PaChr${i}.sf2.sites.txt | sed "s/^/$i/"
done > ./AllChr.sf2.sites.txt

####instead this
(echo -e "position\tx\tn\tfolded"; awk 'FNR>1 { $1=""; sub(/^[ \t]+/, ""); print ++n "\t" $0 }' OFS='\t' ./PaChr{1..25}.sf2.sites.folded.txt) > ./AllChr.sf2.sites.folded.txt
###for swapped Derived allele
(echo -e "position\tx\tn\tfolded"; awk 'FNR>1 { $1=""; sub(/^[ \t]+/, ""); print ++n "\t" $0 }' OFS='\t' ./PaChr{1..25}.sf2.sites.swapped.txt) > ./AllChr.sf2.sites.swapped.txt
#for EU population
(echo -e "position\tx\tn\tfolded"; awk 'FNR>1 { $1=""; sub(/^[ \t]+/, ""); print ++n "\t" $0 }' OFS='\t' ./EU_sourcePaChr{1..25}.sf2.sites.swapped.txt ) > ./EU_sourceAllChr.sf2.sites.swapped.txt


./SF2/SweepFinder2 -f ./AllChr.sf2.sites.folded.txt ./AllChr.sf2.sites.folded.spect
#for swapped
./SF2/SweepFinder2 -f ./AllChr.sf2.sites.swapped.txt ./AllChr.sf2.sites.swapped.spect
#for EU swapped
./SF2/SweepFinder2 -f ./EUsource/EU_sourceAllChr.sf2.sites.swapped.txt ./EU_sourceAllChr.sf2.sites.swapped.spect
#angsd -bam EU_invasive.bamlist -anc ancestral.fasta.fa.gz -out EU_invasive_monor -dosaf 1 -GL 1 -doMajorMinor 1 -doMaf 1
#SweepFinder2 -f EU_invasive.sfs EU_invasive.sf2.sites.txt ./SweepFinder2/EU_invasive_output.txt
#sed -i '1ichrom\tposition\tx\tn' EU_invasive.sf2.sites.txt
#for i in {2..25}; do
 #   ./SweepFinder2 -f ../PaChr$i.sf2.sites.txt ../PaChr$i.sf2.spect
#done
for i in {1..25}; do
    ./SF2/SweepFinder2 -lg 5000 ./PaChr$i.sf2.sites.folded.txt ./AllChr.sf2.sites.folded.spect ./PaChr$i.sf2.olded_output.txt
done
#for swapped
for i in {1..25}; do
    ./SF2/SweepFinder2 -lg 5000 ./PaChr$i.sf2.sites.swapped.txt ./AllChr.sf2.sites.swapped.spect ./PaChr$i.sf2.swapped_output.txt
done
#For EU population
for i in {21..25}; do
    ../SF2/SweepFinder2 -lg 5000 ./EU_sourcePaChr$i.sf2.sites.swapped.txt ./EU_sourceAllChr.sf2.sites.swapped.spect ./EU_sourcePaChr$i.sf2.swapped_output.txt
done

module load r-env/451
##this is non-sense, be careful when using it, especially the snp sites, modified by chr number
for i in {1..25}; do
    sed '1d' ./EU_sourcePaChr${i}.sf2.sites.txt | sed "s/^/$i/"
done > ./EU_source_AllChr.sf2.sites.txt
#nano add header, 
./SweepFinder2 -f /scratch/Cui3/Salt/VariantWGS_1/Selection/EUsource/EU_source_AllChr.sf2.sites.txt /scratch/Cui3/Salt/VariantWGS_1/Selection/EUsource/EU_source_AllChr.sf2.sites.spect

for i in {1..25}; do
    ./SweepFinder2 -lg 5000 ../EUsource/EU_sourcePaChr$i.sf2.sites.txt ../EUsource/EU_source_AllChr.sf2.sites.spect ../EUsource/EU_sourcePaChr$i.sf2_output.txt
done

awk 'NR>3 {print $5 "\t" $6 "\t" $7 "\t" $10 "\t" $11}' Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta.out | grep 'PaChr' > repeat.bed

################extra high, LR>10
> high_LR.bed  # 先清空文件
for i in {1..25}; do
  awk -v i="$i" '$2>10 { 
      start=int($1-2500); 
      if(start>=0) {
          end=int($1+2500); 
          print "PaChr"i "\t" start "\t" end
      }
  }' PaChr$i.sf2.swapped_output.txt >> high_LR.bed
done
> low_LR.bed  # 先清空文件
for i in {1..25}; do
  awk -v i="$i" '$2<10 { 
      start=int($1-2500); 
      if(start>=0) {
          end=int($1+2500); 
          print "PaChr"i "\t" start "\t" end "\t" $2
      }
  }' PaChr$i.sf2.swapped_output.txt >> low_LR.bed
done

################extra high, LR>50
> high_LRchr1_25.bed  # 先清空文件
for i in {1..25}; do
  awk -v i="$i" '$2>10 { 
      start=int($1-2500); 
      if(start>=0) {
          end=int($1+2500); 
          print "PaChr"i "\t" start "\t" end
      }
  }' EU_sourcePaChr$i.sf2.swapped_output.txt >> high_LRchr1_25.bed
done

#bedtools intersect -a high_LR.bed -b repeat.bed -wo > high_LR_in_repeat.bed
bedtools intersect -a high_LR.bed -b high_LRchr1_25.bed -wa |sort | uniq | wc -l #36
bedtools intersect -a high_LR.bed -b high_LRchr1_25.bed -wa > EU_inv_overlap.bed
bedtools intersect -a high_LRchr1_25.bed -b repeat.bed -wa
#bedtools intersect -a high_LR.bed -b repeat.bed -wb | sort | uniq| wc -l
bedtools intersect -a EU_inv_overlap.bed -b gene.tab.bed -wa | sort | uniq| wc -l
bedtools intersect -a high_LR.bed -b gene.tab.bed -wb | sort | uniq| wc -l #102
bedtools intersect -a high_LR.bed -b gene.tab.bed -wb | grep 'PaChr24' | sort | uniq| wc -l #19
bedtools intersect -a high_LRchr1_25.bed -b gene.tab.bed -wb | sort | uniq| wc -l #73

cp /scratch/Cui3/singlecellanalysis/PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 .
bedtools intersect -a high_LR.bed -b PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -wb | grep '\bgene\b' #130 genes
bedtools intersect -a high_LR.bed -b PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -wb | grep '\bgene\b' | grep 'PaChr24' #13 genes

grep 'PaChr24' PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 | grep -w -f HighExpBchr.list | grep '\bgene\b' | awk '{ print $1"\t"$4"\t"$5 }'> HighExpBchr.bed
#grep 'PaChr24' PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 | grep -w -f selected.list | grep '\bgene\b' | awk '{ print $1"\t"$4"\t"$5 }'> selected.bed
cp HighExpBchr.bed /scratch/Cui3/Salt/VariantWGS_1/Selection
#cp selected.bed /scratch/Cui3/Salt/VariantWGS_1/Selection
bedtools intersect -a selected.bed -b low_LR.bed -wb

#check coverage
bedtools intersect -a high_LR.bed -b 66cov1000bin.txt -wa > high_LR_coverage_66.txt
df1 <- read.csv("66cov1000bin.txt", header = FALSE, sep="\t")
df2 <- read.csv("high_LR_coverage_66.txt", header = FALSE, sep="\t")
cov1 <- df1[, 4]
cov2 <- df2[, 7]
test_result <- wilcox.test(cov1, cov2, alternative = "two.sided")
print(test_result)
#
        Wilcoxon rank sum test with continuity correction
data:  cov1 and cov2
W = 571392979, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0

library(readr)
library(ggplot2)

# 读取数据
df <- read_tsv("PaChr24.sf2.swapped_output.txt")

# 绘图
p<-ggplot(df, aes(x = location, y = LR)) +
  geom_point(color = "red", size = 0.8, alpha = 0.7) +
  labs(
    x = "Position on PaChr24",
    y = "LR (selective sweep)",
    title = "PaChr24 selective sweep signal"
  ) +
  theme_bw(base_size = 12) +                 # 纯白背景
  theme(
    panel.grid = element_blank(),            # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5)   # 标题居中
  )

# 保存成PNG
ggsave("PaChr24_LR_swapped_dot.png", plot = p, width = 8, height = 4, dpi = 300)

library(readr)
library(ggplot2)

df <- read_tsv("PaChr24.sf2.swapped_output.txt")

ggplot(df, aes(x = location/1e6, y = LR)) +
  geom_point(aes(color = LR > 10), size = 1.2, alpha = 0.9) +
  scale_color_manual(values = c("FALSE" = "#d88900",
                                "TRUE"  = "black")) +
  labs(
    x = "Position (Mb)",
    y = "CLR"
  ) +
  theme_bw(base_size = 22) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave("PaChr24_CLR_bigfont.png", width = 8, height = 4.5, dpi = 300)


library(readr)
library(ggplot2)
library(dplyr)

# 读入数据
df_lr <- read_tsv("PaChr24.sf2_output.txt")
df_tajima <- read_tsv("PaChr24.TajimaD.txt")

# 合并数据
df_lr$stat <- "LR"
df_tajima$stat <- "TajimaD"
df_lr <- df_lr %>% select(location, value = LR, stat)
df_tajima <- df_tajima %>% select(location, value = TajimaD, stat)
df <- bind_rows(df_lr, df_tajima)

# 绘图
p <- ggplot(df, aes(x = location, y = value, color = stat)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("LR" = "red", "TajimaD" = "blue")) +
  labs(
    x = "Position on PaChr24",
    y = "Statistic value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5)
  )

# 保存图片
ggsave("PaChr24_LR_TajimaD.pdf", plot = p, width = 8, height = 4, dpi = 300)













以下为chatgpt乱写的
# ----------------- 配置 -----------------
sf_file     <- "PaChr24.sf2_output.txt"      # SweepFinder2 输出
genes_bed   <- "gene.tab.bed"                # 基因 bed 文件
repeats_bed <- "repeat.bed"                  # RepeatMasker 输出转成 BED
chrom       <- "PaChr24"                     # 目标染色体
out_pdf     <- "PaChr24_LR_Gypsy_Copia.pdf" # 输出 PDF 文件
window_size <- 5000
# ---------------------------------------

library(GenomicRanges)
library(dplyr)
library(ggplot2)

# --- 1. 读取数据 ---
sf <- read.table(sf_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
genes_all <- read.table(genes_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(genes_all) <- c("chr","start","end")

repeats_all <- read.table(repeats_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if(ncol(repeats_all) >= 5){
  colnames(repeats_all)[1:5] <- c("chr","start","end","family","class")
} else {
  stop("repeat.bed 列数不足 5 列，请检查输入")
}

# --- 2. 过滤目标染色体 ---
genes <- genes_all %>% filter(chr == chrom)
repeats <- repeats_all %>% filter(chr == chrom)

# --- 3. TE 分类 ---
repeats <- repeats %>%
  mutate(
    class_lower = tolower(class),
    TE_type = case_when(
      grepl("ltr/gypsy", class_lower) ~ "Gypsy",
      grepl("ltr/copia", class_lower) ~ "Copia",
      grepl("^ltr/", class_lower) ~ "OtherLTR",  # 其他 LTR
      TRUE ~ "OtherTE"                           # 非 LTR
    )
  )

# --- 4. 构建 GenomicRanges ---
sf_gr <- GRanges(
    seqnames = chrom,
    ranges = IRanges(
        start = sf$location - window_size/2 + 1,  # +1 because IRanges is inclusive
        end   = sf$location + window_size/2
    )
)
sf_gr     <- GRanges(seqnames = chrom,
                     ranges = IRanges(start = floor(sf$location),
                                      end = ceiling(sf$location)))
gene_gr   <- GRanges(seqnames = genes$chr,
                     ranges = IRanges(start = genes$start, end = genes$end))
repeat_gr <- GRanges(seqnames = repeats$chr,
                     ranges = IRanges(start = repeats$start, end = repeats$end))

# --- 5. 判断是否落在 TE 和基因区 ---
sf$TE_hit <- FALSE
sf$TE_type <- "None"

ov_rep <- findOverlaps(sf_gr, repeat_gr)
if (length(ov_rep) > 0) {
  ov_df <- as.data.frame(ov_rep)
  ov_df$TE_type <- repeats$TE_type[ov_df$subjectHits]
  
  # 优先级 Gypsy > Copia > OtherLTR > OtherTE
  by_sf <- ov_df %>% group_by(queryHits) %>%
    summarize(TE_pick = case_when(
      any(TE_type == "Gypsy") ~ "Gypsy",
      any(TE_type == "Copia") ~ "Copia",
      any(TE_type == "OtherLTR") ~ "OtherLTR",
      TRUE ~ "OtherTE"
    ))
  
  sf$TE_hit[by_sf$queryHits] <- TRUE
  sf$TE_type[by_sf$queryHits] <- by_sf$TE_pick
}

sf$in_gene <- countOverlaps(sf_gr, gene_gr) > 0

# --- 6. 区域类型标记 ---
sf$region_type <- "intergenic"

# 基因内
sf$region_type[sf$in_gene & sf$TE_hit & sf$TE_type=="Gypsy"]     <- "gene+Gypsy"
sf$region_type[sf$in_gene & sf$TE_hit & sf$TE_type=="Copia"]     <- "gene+Copia"
sf$region_type[sf$in_gene & sf$TE_hit & sf$TE_type=="OtherLTR"]  <- "gene+OtherLTR"
sf$region_type[sf$in_gene & sf$TE_hit & sf$TE_type=="OtherTE"]   <- "gene+OtherTE"
sf$region_type[sf$in_gene & !sf$TE_hit]                           <- "gene"

# 基因外
sf$region_type[!sf$in_gene & sf$TE_hit & sf$TE_type=="Gypsy"]     <- "Gypsy"
sf$region_type[!sf$in_gene & sf$TE_hit & sf$TE_type=="Copia"]     <- "Copia"
sf$region_type[!sf$in_gene & sf$TE_hit & sf$TE_type=="OtherLTR"]  <- "OtherLTR"
sf$region_type[!sf$in_gene & sf$TE_hit & sf$TE_type=="OtherTE"]   <- "OtherTE"

# --- 7. 绘图 ---
pdf(out_pdf, width = 12, height = 6)
ggplot(sf, aes(x = location, y = LR)) +
  geom_line(color = "grey70") +
  geom_point(aes(color = region_type), size = 0.7, alpha = 0.8) +
  scale_color_manual(values = c(
    "gene+Gypsy"    = "purple",
    "gene+Copia"    = "magenta",
    "gene+OtherLTR" = "pink",
    "gene+OtherTE"  = "brown",
    "gene"          = "blue",
    "Gypsy"         = "red",
    "Copia"         = "orange",
    "OtherLTR"      = "darkgreen",
    "OtherTE"       = "grey50",
    "intergenic"    = "black"
  )) +
  theme_minimal() +
  labs(title = paste0(chrom, " LR with TE overlap"),
       x = "Position (bp)", y = "LR", color = "Region") +
  theme(legend.position = "right")
dev.off()

cat("Plot saved to:", out_pdf, "\n")
# 筛选 LR > 50
sf_high <- sf[sf$LR > 50, ]
sf_high_gr <- sf_gr[sf$LR > 50]

# 找 overlaps，获取基因“名称”（坐标）
ov <- findOverlaps(sf_high_gr, gene_gr)
overlap_genes <- gene_gr[subjectHits(ov)]
overlap_genes_df <- data.frame(
  chr   = as.character(seqnames(overlap_genes)),
  start = start(overlap_genes),
  end   = end(overlap_genes)
)
overlap_genes_df
#> overlap_genes_df
#      chr    start      end
#1 PaChr24   537069   537230 #evm.TU.PaChr24.18
#2 PaChr24  1042146  1044229 #evm.model.PaChr24.47
#3 PaChr24  9771006  9777202 #evm.TU.PaChr24.101
#4 PaChr24 10074926 10079339  #evm.model.PaChr24.110
#5 PaChr24 11865260 11870949 #evm.model.PaChr24.173

#if use windows +-2500bp, LRhigher than 50, with overlap at least 1000bp, then the result is:
  window_index   location        LR gene_start gene_end overlap_len
1            3   537084.9  90.79605     535615   536872        1258
2            5   577091.2  77.73239     578369   579816        1223
3            7  1042164.0  53.59228    1042147  1044229        2083
4           44  9773530.0  57.55439    9771007  9777202        5000
5           46 10078577.7 100.45597   10074927 10079339        3262
6           56 11553808.5 109.42569   11550453 11552313        1005
7           58 11863857.0  56.14621   11865261 11870949        1097
8           59 11868857.8 117.29013   11865261 11870949        4592
9           60 11893861.7  80.72489   11894606 11896295        1690
#if use windows +-5000bp, LRhigher than 50, with overlap at least 500bp, then the result is:
> sf_gene_overlap
   window_index   location        LR gene_start gene_end overlap_len
1             2   532084.2  99.97181     535615   536872        1258
2             3   537084.9  90.79605     535615   536872        1258
3             5   577091.2  77.73239     578369   579816        1448
4             7  1042164.0  53.59228    1042147  1044229        2083
5             8  1047164.7 151.47773    1042147  1044229        2066
6            11  1082170.2  59.64755    1084098  1084828         731
7            44  9773530.0  57.55439    9771007  9777202        6196
8            46 10078577.7 100.45597   10074927 10079339        4413
9            55 11548807.7 124.58793   11550453 11552313        1861
10           56 11553808.5 109.42569   11550453 11552313        1861
11           57 11573811.6  70.59861   11572987 11573755         769
12           58 11863857.0  56.14621   11865261 11870949        3597
13           59 11868857.8 117.29013   11865261 11870949        5689
14           59 11868857.8 117.29013   11871057 11873841        2785
15           60 11893861.7  80.72489   11894606 11896295        1690
16           61 11903863.3 218.18003   11900399 11901028         630
17           61 11903863.3 218.18003   11902692 11903208         517
awk '{ print "PaChr24\t"$5"\t"$6 }' sf_gene_overlap.txt | sed '/PaChr24[[:space:]]*$/d' > sf_gene_overlap.bed
rep 'mRNA' /scratch/project_2009273/Cui3/singlecellanalysis/PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 | grep 'PaChr24' | awk '{ print $1"\t"$4"\t"$5"\t"$9 }' > B_gene.gff
#cd /scratch/project_2009273/Cui3/singlecellanalysis/
#bedtools intersect -wb -a sf_gene_overlap.bed -b B_gene.gff 
bedtools intersect -wb -a sf_gene_overlap.bed -b B_gene.gff | awk '{ print $7 }' | uniq
#ID=evm.model.PaChr24.17;Parent=evm.TU.PaChr24.17;gene_id=evm.TU.PaChr24.17;gene_name=EVM%2520prediction%2520PaChr24.17;transcript_id=evm.model.PaChr24.17
evm.TU.PaChr24.19
evm.TU.PaChr24.47
evm.TU.PaChr24.48
evm.TU.PaChr24.101
evm.TU.PaChr24.110
evm.TU.PaChr24.160
evm.TU.PaChr24.161
evm.TU.PaChr24.173
evm.TU.PaChr24.174
evm.TU.PaChr24.175
evm.TU.PaChr24.176
evm.TU.PaChr24.177

# --- 8. 提取 LR > 50 且 intergenic 的窗口 ---
sf_intergenic_high <- sf %>%
  filter(LR > 50, region_type == "intergenic")

cat("Number of intergenic windows with LR > 50:", nrow(sf_intergenic_high), "\n")

# --- 9. 计算每个 intergenic window 距离最近基因的距离 ---
sf_intergenic_high_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = sf_intergenic_high$location - window_size / 2,
    end   = sf_intergenic_high$location + window_size / 2
  )
)

nearest_gene <- distanceToNearest(sf_intergenic_high_gr, gene_gr)

# 整理结果
sf_intergenic_high$nearest_gene_dist <- mcols(nearest_gene)$distance[match(seq_along(sf_intergenic_high_gr), queryHits(nearest_gene))]

# --- 10. 过滤出距离基因 <= 2000 bp 的窗口（靠近基因）---
sf_intergenic_near <- sf_intergenic_high %>%
  filter(nearest_gene_dist <= 5000)

cat("Number of intergenic LR>50 windows within 2kb of genes:", nrow(sf_intergenic_near), "\n")

# --- 11. 保存 BED 文件 ---
sf_intergenic_near_bed <- data.frame(
  chr   = chrom,
  start = as.integer(floor(sf_intergenic_high$location - window_size / 2)),
  end   = as.integer(ceiling(sf_intergenic_high$location + window_size / 2))
  #LR    = sf_intergenic_near$LR,
  #dist2gene = sf_intergenic_near$nearest_gene_dist
)
# 去除NA、start>=end的错误行
sf_bed_clean <- sf_intergenic_near_bed[!is.na(sf_intergenic_near_bed$start) &
                             !is.na(sf_intergenic_near_bed$end), ]

write.table(
  sf_bed_clean,
  file = paste0(chrom, "_intergenic_LR50_near5kb.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Near (<5kb) intergenic high-LR windows saved to:", paste0(chrom, "_intergenic_LR50_near5kb.bed"), "\n")
cp /scratch/project_2009273/Cui3/Salt/VariantWGS_1/Selection/PaChr24_intergenic_LR50_near5kb.bed .
#nano manulaly remove header  in PaChr24_intergenic_LR50_near5kb.bed
bedtools intersect -wb -a PaChr24_intergenic_LR50_near5kb.bed -b B_gene.gff | awk '{ print $7 }' | uniq
#evm.TU.PaChr24.19;
evm.TU.PaChr24.86;
evm.TU.PaChr24.160;
evm.TU.PaChr24.161;
evm.TU.PaChr24.175;