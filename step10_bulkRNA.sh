cd /scratch/project_2009273/Cui3/RNAseq_EU
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//g'> namelist
#because some study showed soft cliping of adapters and low quality sliding windows are not helpful, even harmful, so we only keep the paired reads, whichi is nearly all
trimmomatic PE -threads 40 ${name}"_R1.fastq.gz" ${name}"_R2.fastq.gz"  MINLEN:50 -baseout ${name}"_trimmed"
#mapping
ln -s /scratch/project_2009273/Cui3/singlecellanalysis/Pa.ChrandChloroMito.fasta .
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir ./ --genomeFastaFiles Pa.ChrandChloroMito.fasta --genomeSAindexNbases 13
STAR --runMode alignReads --runThreadN 40 --genomeDir ./ --readFilesIn ${name}"_trimmed_1P" ${name}"_trimmed_2P" --outFilterMismatchNmax 15 --outSAMstrandField intronMotif 
samtools sort -@ 40 Aligned.out.sam -o Aligned.out.bam.sort -O BAM
export PATH="/projappl/project_2009273/stringtie/bin:$PATH"
stringtie -o $name".out.gtf" -G ../PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -p 40 -A $name".abundance" -l $name -b $name marked_duplicates.bam
#stringtie --merge -G PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -o merged.gtf -p 28 merge.list
#stringtie -e -p 28 -G ../merged.gtf -o ${name}".merged.gtf" -A ${name}".abundance"  -l ${name} -b ${name} marked_duplicates.bam
stringtie -e -p 40 -o $name".out.e.gtf" -G ../PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -A $name".abundance" -l $name -b $name marked_duplicates.bam

export PROJAPPL=/projappl/project_2001259 
python prepDE.py3 -i sample_lst.txt

start-r
 .libPaths(c("/projappl/project_2001259/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(DESeq2)
library(ggplot2)
library(ggrepel)
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
condition<-factor(c('invasive','invasive','invasive','invasive','invasive','EU','EU','EU','EU'))
colData<-data.frame(row.names=colnames(countData),condition)
all(rownames(colData) %in% colnames(countData))#check the data
countData <- na.omit(countData)
dds_condition <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#PCA
rld <- rlog(dds_condition)
pdf("PCA_condition.pdf")

all(rownames(colData) %in% colnames(countData))#check the data
countData<- countData[, rownames(colData)]#for condition
all(rownames(colData) == colnames(countData))#check dataset
dds_condition <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds_condition <- dds_condition[rowSums(counts(dds_condition)) > 1, ]
#setup the reference to make comparisons
dds_condition$condition <- relevel( dds_condition$condition, ref ="EU" )
dds_dds_condition <- DESeq(dds_condition)
res_condition <- results(dds_dds_condition)
#to compare OL and OC
res_dds_condition<- results( dds_dds_condition, contrast = c("condition","invasive", "EU"))
res_dds_condition_Ordered <- res_dds_condition[order(res_dds_condition$pvalue),]
summary(res_dds_condition_Ordered)
#heatmap
rld_dds_condition_Ordered <- rlog(dds_dds_condition)
sum(res_dds_condition_Ordered$padj < 0.001, na.rm=TRUE)
library("pheatmap")
library("genefilter")
#topVarGenes <- head(order(-rowVars(assay(rld))),500)
topVarGenes <- head(order(-rowVars(assay(rld_dds_condition_Ordered))),100)
mat <- assay(rld_dds_condition_Ordered)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld_dds_condition_Ordered)[,"condition"])
rownames(df) <- colnames(mat)
pdf('heatmap100_invasive_highvaria.pdf',width = 6, height = 7)
pheatmap(mat, annotation_col=df,fontsize = 5)
dev.off()

DEG_condition <-subset(res_dds_condition_Ordered,padj < 0.05 & (log2FoldChange > 2))
DEG_condition_out <- cbind(row.names(DEG_condition),DEG_condition$baseMean,DEG_condition$log2FoldChange,
DEG_condition$lfcSE,DEG_condition$stat,DEG_condition$pvalue,DEG_condition$padj)
colnames(DEG_condition_out)<-c('gene',colnames(DEG_condition))
res_dds_condition_Ordereddata <-  merge(as.data.frame(res_dds_condition_Ordered),as.data.frame(counts(dds_dds_condition,normalize=TRUE)),by="row.names",sort=FALSE)
List<-NULL;
for(i in 1:dim(DEG_condition_out)[1]){
  gene<-DEG_condition_out[i,1]
  p<-res_dds_condition_Ordereddata[which(res_dds_condition_Ordereddata$Row.names==gene),]
  List<-rbind(List,p)
}
write.csv(List,file= "InvasivevsEU.difference-gene.p0.05-f2.upregulate.result.csv",row.names = F)#

#downregulation
DEG_condition <-subset(res_dds_condition_Ordered,padj < 0.05 & (log2FoldChange < -2))
DEG_condition_out <- cbind(row.names(DEG_condition),DEG_condition$baseMean,DEG_condition$log2FoldChange,
DEG_condition$lfcSE,DEG_condition$stat,DEG_condition$pvalue,DEG_condition$padj)
colnames(DEG_condition_out)<-c('gene',colnames(DEG_condition))
res_dds_condition_Ordereddata <-  merge(as.data.frame(res_dds_condition_Ordered),as.data.frame(counts(dds_dds_condition,normalize=TRUE)),by="row.names",sort=FALSE)
List<-NULL;
for(i in 1:dim(DEG_condition_out)[1]){
  gene<-DEG_condition_out[i,1]
  p<-res_dds_condition_Ordereddata[which(res_dds_condition_Ordereddata$Row.names==gene),]
  List<-rbind(List,p)
}
write.csv(List,file= "InvasivevsEU.difference-gene.p0.05-f2.downregulate.result.csv",row.names = F)

awk -F, '{ print $1 }' InvasivevsEU.difference-gene.p0.01-f2.upregulate.result.csv | sed 's/"//g' | awk -F'|' '{ print $1}' | sed 's/TU/model/g' | sed '1d' > InvasivevsEU.difference-gene.p0.01-f2.upregulate.gene
awk -F, '{ print $1 }' InvasivevsEU.difference-gene.p0.01-f2.downregulate.result.csv | sed 's/"//g' | awk -F'|' '{ print $1}' | sed 's/TU/model/g' | sed '1d' > InvasivevsEU.difference-gene.p0.01-f2.downregulate.gene
awk -F, '{ print $1 }' InvasivevsEU.difference-gene.p0.05-f2.upregulate.result.csv | sed 's/"//g' | awk -F'|' '{ print $1}' | sed 's/TU/model/g' | sed '1d' > InvasivevsEU.difference-gene.p0.05-f2.upregulate.gene
awk -F, '{ print $1 }' InvasivevsEU.difference-gene.p0.05-f2.downregulate.result.csv | sed 's/"//g' | awk -F'|' '{ print $1}' | sed 's/TU/model/g' | sed '1d' > InvasivevsEU.difference-gene.p0.05-f2.downregulate.gene
mv InvasivevsEU.difference-gene.p0.05-f2.upregulate.gene /scratch/project_2001259/cuiwang/nodelete5/CNWGS2/GO
mv InvasivevsEU.difference-gene.p0.05-f2.downregulate.gene /scratch/project_2001259/cuiwang/nodelete5/CNWGS2/GO
cd /scratch/project_2001259/cuiwang/nodelete5/CNWGS2/GO
python /scratch/project_2001259/cuiwang/nodelete5/CNWGS/ChrBenrichment/goatools/scripts/find_enrichment.py InvasivevsEU.difference-gene.p0.05-f2.upregulate.gene population.txt Allgene.association --pval=0.05 --outfile=InvasivevsEU.difference-gene.p0.05-f2.upregulate.gene"_enrichment.txt";
python /scratch/project_2001259/cuiwang/nodelete5/CNWGS/ChrBenrichment/goatools/scripts/find_enrichment.py InvasivevsEU.difference-gene.p0.05-f2.downregulate.gene population.txt Allgene.association --pval=0.05 --outfile=InvasivevsEU.difference-gene.p0.05-f2.downregulate.gene"_enrichment.txt";
