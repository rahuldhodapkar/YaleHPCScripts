#!/usr/bin/env bash
#
# Script for processing a single batch of TCRs and writing the output to
# an appropriate place.
#

set -e
set -x

BASE_DIR=`pwd`
DATA_DIR=$1
SAMPLE_NAME=`basename $DATA_DIR`

cd $DATA_DIR
CELL_NAMES=`ls | sed -r 's/([^_]*).*$/\1/' | uniq`

for cell in $CELL_NAMES
do
    echo "===== BEGIN $cell ====="
    mkdir -p $BASE_DIR/$SAMPLE_NAME/$cell
    cd $BASE_DIR/$SAMPLE_NAME/$cell

    PE1_FQ=`ls -d $DATA_DIR/* | grep $cell | grep "R1"`
    PE2_FQ=`ls -d $DATA_DIR/* | grep $cell | grep "R2"`

    mixcr analyze shotgun \
        --species hsa \
        --starting-material rna \
        --only-productive \
        --contig-assembly \
        <(zcat $PE1_FQ) \
        <(zcat $PE2_FQ) \
        analysis
    echo "===== END $cell ====="
done


