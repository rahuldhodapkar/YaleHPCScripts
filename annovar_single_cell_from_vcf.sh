#!/usr/bin/env bash
#
# Annotate a single cell from VCF output
#

set -e
set -x

TARGET_FOLDER=$1
cd $TARGET_FOLDER

annovar_files=/home/rd389/scratch60/refseq/annovar/

convert2annovar.pl -format vcf4 \
    Aligned.sortedByCoord.out.raw.vcf \
    -outfile raw.avinput

annotate_variation.pl \
    -out anno \
    -build human_g1k_v37_decoy_ercc \
    raw.avinput \
    $annovar_files

