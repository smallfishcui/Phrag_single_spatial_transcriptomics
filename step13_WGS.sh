module load angsd
cd /scratch/Cui3/Salt/VariantWGS_1/Fst
angsd -b EU_source.bamlist  -anc ancestral.fasta.fa.gz -out EU_source -dosaf 1 -gl 1
angsd -b EU_invasive.bamlist  -anc ancestral.fasta.fa.gz -out EU_invasive -dosaf 1 -gl 1
#EU_source.bamlist
/scratch/Cui3/Salt/Ploidy/94.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/92.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/72.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/71.marked_duplicates.filter.bam
#EU_invasive.bamlist
/scratch/Cui3/Salt/Ploidy/66.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/103.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/68.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/70.marked_duplicates.filter.bam
/scratch/Cui3/Salt/Ploidy/62.marked_duplicates.filter.bam

realSFS EU_invasive.saf.idx
realSFS EU_source.saf.idx
realSFS print EU_invasive.saf.idx
realSFS EU_invasive.saf.idx -maxIter 100 -P 8 >EU_invasive.sfs

gunzip EU_invasive_monor.mafs.gz
#sweepFinder2
#install via tykky
#cd /projappl/project_2009273
#conda-containerize new --mamba --prefix SweepFinder2 env.yml
angsd -bam EU_invasive.bamlist -anc ancestral.fasta.fa.gz -out EU_invasive_monor -doCounts 1 -dosaf 1 -GL 1 -doMajorMinor 1 -doMaf 8
realSFS print EU_invasive_monor.saf.idx > EU_invasive_monor.saf.txt

angsd -bam EU_source.bamlist -anc ancestral.fasta.fa.gz -out EU_source_minor -doCounts 1 -dosaf 1 -GL 1 -doMajorMinor 1 -doMaf 8
realSFS print EU_source_minor.saf.idx > EU_source_minor.saf.txt

cd /scratch/Cui3/Salt/VariantWGS_1/Fst

cat EU_invasive_monor.mafs | \
awk '
NR>1 && $1 ~ /^PaChr/ {
    chrom=$1; pos=$2; phat=$6; nInd=$7;
    if(nInd==0) next;
    x=int(phat * 2 * nInd + 0.5);   # 衍生等位基因数
    if(x<=0) next;                   # 去掉 x=0 的位点
    n=2*nInd;                        # 总等位基因数
    if(x>n) x=n;
    # 输出到对应染色体文件，如果文件不存在先加 header
    file=chrom".sf2.sites.txt"
    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file]=1
    }
    print pos"\t"x"\t"n"\t0" >> file
}'

zcat EU_source_minor.mafs.gz | \
awk '
NR>1 && $1 ~ /^PaChr/ {
    chrom=$1
    pos=$2
    phat=$6
    nInd=$7
    if (nInd == 0) next
    # 衍生等位基因数 (四舍五入)
    x = int(phat * 2 * nInd + 0.5)
    if (x <= 0) next    # 去掉 x=0 的位点
    n = 2 * nInd        # 总等位基因数
    if (x > n) x = n
    # 输出文件名（加前缀 EU_source）
    file = "EU_source" chrom ".sf2.sites.txt"
    # 如果是第一次写这个文件，加表头
    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file] = 1
    }
    # 写数据
    print pos "\t" x "\t" n "\t0" >> file
}'

mv PaChr*.sf2.sites.txt /scratch/Cui3/Salt/VariantWGS_1/Selection/
mv /scratch/Cui3/Salt/VariantWGS_1/Selection/EUsource

#treat them as folded
cat EU_invasive_monor.mafs | \
awk '
NR>1 && $1 ~ /^PaChr/ {
    chrom=$1; pos=$2; phat=$6; nInd=$7;
    if(nInd==0) next;
    x=int(phat * 2 * nInd + 0.5);   # 衍生等位基因数
    if(x<=0) next;                   # 去掉 x=0 的位点
    n=2*nInd;                        # 总等位基因数
    if(x>n) x=n;
    # 输出到对应染色体文件，如果文件不存在先加 header
    file=chrom".sf2.sites.folded.txt"
    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file]=1
    }
    print pos"\t"x"\t"n"\t1" >> file
}'

zcat EU_source_minor.mafs.gz | \
awk '
NR>1 && $1 ~ /^PaChr/ {
    chrom=$1
    pos=$2
    phat=$6
    nInd=$7
    if (nInd == 0) next
    # 衍生等位基因数 (四舍五入)
    x = int(phat * 2 * nInd + 0.5)
    if (x <= 0) next    # 去掉 x=0 的位点
    n = 2 * nInd        # 总等位基因数
    if (x > n) x = n
    # 输出文件名（加前缀 EU_source）
    file = "EU_source" chrom ".sf2.sites.folded.txt"
    # 如果是第一次写这个文件，加表头
    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file] = 1
    }
    # 写数据
    print pos "\t" x "\t" n "\t1" >> file
}'

mv PaChr*.sf2.sites.folded.txt /scratch/Cui3/Salt/VariantWGS_1/Selection/
mv /scratch/Cui3/Salt/VariantWGS_1/Selection/EUsource

#separate folded and unfolded,convert all folded to unfolded
cat EU_invasive_monor.mafs | \
awk '
NR > 1 && $1 ~ /^PaChr/ {
    chrom = $1
    pos   = $2
    major = $3
    minor = $4
    anc   = $5
    phat  = $6
    nInd  = $7

    if (nInd == 0) next

    x = int(phat * 2 * nInd + 0.5)
    n = 2 * nInd

    if (x <= 0) next
    if (x > n) x = n

    # minor == anc → folded = 1
    folded = (minor == anc ? 1 : 0)

    # 如果 folded=1，折叠 x 然后将 folded 设为 0
    if (folded == 1) {
        x = (x < (n - x) ? x : (n - x))
        folded = 0
    }

    file = chrom ".sf2.sites.swapped.txt"

    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file] = 1
    }

    print pos "\t" x "\t" n "\t" folded >> file
}'

zcat EU_source_minor.mafs.gz | \
awk '
NR > 1 && $1 ~ /^PaChr/ {
    chrom = $1
    pos   = $2
    major = $3
    minor = $4
    anc   = $5
    phat  = $6
    nInd  = $7

    if (nInd == 0) next

    x = int(phat * 2 * nInd + 0.5)
    n = 2 * nInd

    if (x <= 0) next
    if (x > n) x = n

    # minor == anc → folded = 1
    folded = (minor == anc ? 1 : 0)

    # 如果 folded=1，折叠 x 然后将 folded 设为 0
    if (folded == 1) {
        x = (x < (n - x) ? x : (n - x))
        folded = 0
    }

    file = "EU_source" chrom ".sf2.sites.swapped.txt"

    if (!(file in seen)) {
        print "position\tx\tn\tfolded" > file
        seen[file] = 1
    }

    print pos "\t" x "\t" n "\t" folded >> file
}'
mv PaChr*.sf2.sites.swapped.txt /scratch/Cui3/Salt/VariantWGS_1/Selection/
mv EU_source*.sf2.sites.swapped.txt /scratch/Cui3/Salt/VariantWGS_1/Selection/EUsource
