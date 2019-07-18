#!/usr/bin/env bash
#
# Process a single read pair to VCFs in place.
#
# NOTE: relies on module load infrastructure
#

set -e
set -x

TARGET_FOLDER=$1
cd $TARGET_FOLDER

# set up target files
FQ_FILES=`ls *.fastq.gz`
FQ_TRIM_UNTRIM_FILES=`for f in $FQ_FILES; do filename="${f%%.*}"; echo ${filename}_trim_pair.fastq.gz; echo ${filename}_trim_unpair.fastq.gz; done | tr '\r\n' ' '`

# run trimming
adapter_file=/home/rd389/refseq/illumina/NexteraPE-PE.fa
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
    -threads 40 \
    $FQ_FILES \
    $FQ_TRIM_UNTRIM_FILES \
    ILLUMINACLIP:${adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:5

star_index=/home/rd389/refseq/reference/
sjdb_gtf_file=/home/rd389/refseq/reference/gencode.v26.GRCh38.ERCC.genes.gtf

FQ_TRIM_PAIRED_FILES=`for f in $FQ_FILES; do filename="${f%%.*}"; echo ${filename}_trim_pair.fastq.gz; done | tr '\r\n' ' '`

mkdir -p star

STAR --runThreadN 40 \
     --genomeDir $star_index \
     --readFilesIn $FQ_TRIM_PAIRED_FILES \
     --readFilesCommand zcat \
     --outFileNamePrefix star/ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs RemoveNoncanonical \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.1 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --quantMode TranscriptomeSAM GeneCounts \
     --sjdbGTFfile $sjdb_gtf_file \
     --twopassMode Basic

mv star/Aligned.sortedByCoord.out.bam ./Aligned.sortedByCoord.out.bam
mv star/Aligned.toTranscriptome.out.bam Aligned.toTranscriptome.out.bam
mv star/ReadsPerGene.out.tab ReadsPerGene.out.tab

BAMFILE=Aligned.sortedByCoord.out.bam
FILEBASE=${BAMFILE%.bam}
HUMAN_REFERENCE=/home/rd389/refseq/reference/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta

java -jar ~/bin/picard.jar AddOrReplaceReadGroups \
      I=$BAMFILE \
      O=rg_$BAMFILE \
      SO=coordinate \
      RGID=id \
      RGLB=library \
      RGPL=platform \
      RGPU=machine \
      RGSM=sample 

# Filter optical duplicates for genomic alignment:
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
     I=rg_$BAMFILE \
     O=nodup_$BAMFILE \
     M=star_genome_duplicates.log \
     REMOVE_DUPLICATES=TRUE \
     CREATE_INDEX=true \
     VALIDATION_STRINGENCY=SILENT

# Run GATK
java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T SplitNCigarReads \
      -R $HUMAN_REFERENCE \
      -I nodup_$BAMFILE \
      -o split_$BAMFILE \
      -rf ReassignOneMappingQuality \
      -RMQF 255 \
      -RMQT 60 \
      -U ALLOW_N_CIGAR_READS

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R $HUMAN_REFERENCE \
      -I split_$BAMFILE \
      -dontUseSoftClippedBases \
      -stand_call_conf 20.0 \
      -o $FILEBASE.raw.vcf

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R $HUMAN_REFERENCE \
      -V $FILEBASE.raw.vcf \
      -window 35 \
      -cluster 3 \
      -filterName FS -filter "FS > 30.0" \
      -filterName QD -filter "QD < 2.0" \
      -o $FILEBASE.trimmed.vcf 

