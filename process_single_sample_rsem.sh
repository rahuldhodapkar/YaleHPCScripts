#!/usr/bin/env bash
#
# Perform RSEM quantitation on single sample

set -e
set -x

cd $1

rsem_index=/home/rd389/scratch60/refseq/rsem/human_g1k_v37_decoy

rsem-calculate-expression \
    --no-bam-output \
    --num-threads 8 \
    --bam \
    --paired-end \
    Aligned.toTranscriptome.out.bam \
    $rsem_index \
    rsem_quant
