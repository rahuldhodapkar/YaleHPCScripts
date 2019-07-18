#!/usr/bin/env bash
#
# load the GATK resource bundle reference files and build
# appropriate STAR directories from scratch to prepare
# for downstream processing.
#
# . ~/genome_module_load.sh
# module load ANNOVAR/2018Apr16-foss-2016b-Perl-5.24.1
#

set -e
set -x

REF_BASE_DIR=$1
BROAD_REF_BUNDLE="b37"
BROAD_REF_PREFIX="human_g1k_v37_decoy"

HOST="ftp.broadinstitute.org"
USER="gsapubftp-anonymous"
PASSWD=""

BASE=`pwd`

# Build directory structur
mkdir -p $REF_BASE_DIR/genome
mkdir -p $REF_BASE_DIR/star

cd $REF_BASE_DIR/genome

ftp -i -n $HOST <<HERE
quote USER $USER
quote PASS $PASSWD
cd bundle/$BROAD_REF_BUNDLE
mget ${BROAD_REF_PREFIX}.fasta.gz
quit
HERE

# Load ERCC sequences
curl -OL https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip

# Concatenate to reference
zcat ${BROAD_REF_PREFIX}.fasta.gz \
    | cat - ERCC92.fa \
    > ${BROAD_REF_PREFIX}_ercc.fasta

# Concatenate anntation tracks
http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
cat ERCC92.gtf >>Homo_sapiens.GRCh37.75.gtf
mv Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75_ERCC92.gtf

# Delete original reference files
rm ${BROAD_REF_PREFIX}.fasta.gz

##########################################################################
## PROCESS REFERENCE FILES
##########################################################################

# Load required modules for index generation
. /home/rd389/genome_module_load.sh

# Prepare references for alignment and variant calling
samtools dict -o ${BROAD_REF_PREFIX}_ercc.dict ${BROAD_REF_PREFIX}_ercc.fasta
samtools faidx ${BROAD_REF_PREFIX}_ercc.fasta

# Generate STAR aligner genome directory
cd $BASE/$REF_BASE_DIR

star --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir star \
     --genomeFastaFiles genome/${BROAD_REF_PREFIX}_ercc.fasta \
     --sjdbGTFfile genome/Homo_sapiens.GRCh37.75_ERCC92.gtf \
     --sjdbOverhang 100

# Generate RSEM reference directory
mkdir -p rsem

rsem-prepare-reference \
    --gtf genome/Homo_sapiens.GRCh37.75_ERCC92.gtf \
    --star-sjdboverhang 100 \
    genome/${BROAD_REF_PREFIX}_ercc.fasta \
    rsem/${BROAD_REF_PREFIX}

##########################################################################
## GENERATE FILES FOR ANNOVAR
##########################################################################
# NOTE: requires gtfToGenePred tool.
# curl -OL http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

mkdir -p annovar

gtfToGenePred \
    -genePredExt genome/Homo_sapiens.GRCh37.75_ERCC92.gtf \
    annovar/human_g1k_v37_decoy_ercc_refGene.txt

retrieve_seq_from_fasta.pl \
    --format refGene \
    --seqfile genome/human_g1k_v37_decoy_ercc.fasta \
    annovar/human_g1k_v37_decoy_ercc_refGene.txt \
    --out annovar/human_g1k_v37_decoy_ercc_refGeneMrna.fa

