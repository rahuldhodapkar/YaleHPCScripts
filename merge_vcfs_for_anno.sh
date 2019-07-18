#!/usr/bin/env bash
#
# Merge raw VCFs for downstream annotation.
#
# NOTE: currently hardcoded to build targt files in
#       /home/rd389/scratch60/full_anno/*
#

set -e
set -x

cd /home/rd389/project/FatTRegVariantCall/RawVCF

BATCH_NUMS=`seq 1 16`
TYPES="F B"

for b in $BATCH_NUMS
do
    for t in $TYPES
    do
        ls | grep "^$b$t" | sed -e "s/^/I=/g" | xargs java -jar -Xmx16g ~/bin/picard.jar MergeVcfs O=/home/rd389/scratch60/full_anno/$b$t-merged.vcf
    done
done

