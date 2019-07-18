#!/usr/bin/env bash
#
# Perform Quality Hard Filtering with new GATK as per Broad Documentation
#

set -e
set -x

cd $1

. ~/load_gatk.sh

gatk SelectVariants \
-V Aligned.sortedByCoord.out.raw.vcf \
-select-type SNP \
-O snps.vcf

gatk SelectVariants \
-V Aligned.sortedByCoord.out.raw.vcf \
-select-type INDEL \
-O indels.vcf

gatk VariantFiltration \
-V snps.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O snps_filtered.vcf

gatk VariantFiltration \
-V indels.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O indels_filtered.vcf


# Find Correct Ref SNPS
for f in `ls ~/scratch60/Zeynep/*-1.vcf`
do
    echo $f
    vcfintersect \
        -r ~/scratch60/refseq/genome/human_g1k_v37_decoy_ercc.fasta \
        -v -i $f \
        snps_filtered.vcf \
    | grep "PASS" \
    | wc -l
done

# Find Correct Ref INDELS
for f in `ls ~/scratch60/Zeynep/*-1.vcf`
do
    echo $f
    vcfintersect \
        -r ~/scratch60/refseq/genome/human_g1k_v37_decoy_ercc.fasta \
        -v -i $f \
        indels_filtered.vcf \
    | grep "PASS" \
    | wc -l
done
