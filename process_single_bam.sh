#!/usr/bin/env bash
#
# GATK-process a single .bam file
# 
# NOTE: relies on module load infrastructure and pre-installed binaries
#       available on rd389's Farnam home directory.

. ~/genome_module_load.sh

set -e
set -x

# generate some name hooks

BAMFILE=$1
FILEBASE=${BAMFILE%.bam}

## Prepare RNAseq data

java -jar ~/bin/picard.jar AddOrReplaceReadGroups \
      I=$BAMFILE \
      O=rg_$BAMFILE \
      SO=coordinate \
      RGID=id \
      RGLB=library \
      RGPL=platform \
      RGPU=machine \
      RGSM=sample 

java -Xmx16g -jar ~/bin/picard.jar MarkDuplicates \
      I=rg_$BAMFILE \
      O=dedup_$BAMFILE  \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=SILENT \
      M=output.metrics 

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T SplitNCigarReads \
      -R ~/refseq/hg19/genome_ercc.fa \
      -I dedup_$BAMFILE \
      -o split_$BAMFILE \
      -rf ReassignOneMappingQuality \
      -RMQF 255 \
      -RMQT 60 \
      -U ALLOW_N_CIGAR_READS

# Execute GATK

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R ~/refseq/hg19/genome_ercc.fa \
      -I split_$BAMFILE \
      -dontUseSoftClippedBases \
      -stand_call_conf 20.0 \
      -o $FILEBASE.raw.vcf

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R ~/refseq/hg19/genome_ercc.fa \
      -V $FILEBASE.raw.vcf \
      -window 35 \
      -cluster 3 \
      -filterName FS -filter "FS > 30.0" \
      -filterName QD -filter "QD < 2.0" \
      -o $FILEBASE.trimmed.vcf 
