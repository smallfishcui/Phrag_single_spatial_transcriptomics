cd /scratch/singlecellanalysis
samtools faidx -r chr.txt Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta > Pa_noorganelle.Chr.fasta
cat Chloroplast_path_sequence.fasta Pa_noorganelle.Chr.fasta > Pa.ChrandChloro.fasta
cat embplant_mt.complete.graph1.1.path_sequence.fasta Pa.ChrandChloro.fasta > Pa.ChrandChloroMito.fasta
conda activate liftoff
liftoff -g EVM.all.gff3 ./Pa.ChrandChloroMito.fasta Pa_noorganelle.asm.bp.p_ctg.wrapped.FINAL.sorted.gapclose.finalrenamed.fasta -o ./Pa.ChrandChloroMito.liftoff.gff
module load sqlite
ORIGINALTRANSCRIPTOME=$(ls *.clean | sed 's/.clean//g')
CLEANTRANSCRIPTOME=$(ls *.clean)
CPU=$SLURM_CPUS_PER_TASK
REFERENCE_GENOME=./Pa.ChrandChloroMito.fasta
singularity exec -B ${PWD}:${PWD} -B /scratch/genome/:/scratch/genome/ ../pasa/pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c alignment_config.txt -C -R -g $REFERENCE_GENOME --ALIGNERS gmap,blat,minimap2 -t $CLEANTRANSCRIPTOME -T -u $ORIGINALTRANSCRIPTOME --CPU 32 &> alignment_assembly.log
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta PhragRNA.sqlite.assemblies.fasta  --pasa_transcripts_gff3 PhragRNA.sqlite.pasa_assemblies.gff3 &> transdecoder.log
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg  /usr/local/src/PASApipeline/misc_utilities/pasa_gff3_validator.pl Pa.ChrandChloroMito.liftoff.gff
mv Pa.ChrandChloroMito.liftoff.gff Pa.ChrandChloroMito.liftoff.gff3
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg  /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
     -c alignment_config.txt -g Pa.ChrandChloroMito.fasta \
     -P Pa.ChrandChloroMito.liftoff.gff3
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg cp /usr/local/src/PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt .   
mv pasa.annotationCompare.Template.txt annotCompare.config #change database setting using nano
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg  /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g Pa.ChrandChloroMito.fasta \
        -t transcripts_all.fasta.clean --CPU 32

#another rounnd
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg  /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
     -c alignment_config.txt -g Pa.ChrandChloroMito.fasta \
     -P PhragRNA.sqlite.gene_structures_post_PASA_updates.2473363.gff3
singularity exec -B ${PWD}:${PWD} ../pasa/pasapipeline.v2.5.3.simg  /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g Pa.ChrandChloroMito.fasta \
        -t transcripts_all.fasta.clean --CPU 32
# cell ranger
cd /scratch/cuiwang/nodelete3/CNWGS2/annotation  #puhti
gffread -E PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gff3 -T -o PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf
cd /scratch/singlecellanalysis
cellranger mkref --nthreads=8 --genome=CN --fasta=Pa.ChrandChloroMito.fasta --genes=PhragRNA.sqlite.gene_structures_post_PASA_updates.142237.gtf --jobmode=local
cellranger count --id NAint61 --transcriptome=CN --fastqs /scratch/project_2009273/Cui/singlecell/NAint61 --include-introns true --create-bam true --jobmode=local --localcores 128
#replace name with other samples
cellranger count --id $name --transcriptome=CN --fastqs "/scratch/project_2009273/Cui/singlecell/"$name --include-introns true --create-bam true --jobmode=local --localcores 128
