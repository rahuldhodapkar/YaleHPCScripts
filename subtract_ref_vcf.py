#!/usr/bin/env python3
#
# Subtract variants in a reference file from variants in a target file.
# Distinguish variants only on the basis of CHROM,POS,REF,ALT tuples
#
# Usage:
#
#   ./subtract_ref_vcf.py <ref.vcf> <cell.vcf>
#

import os
import sys

ref_vcf_filename = sys.argv[1]
cell_vcf_filename = sys.argv[2]

fs="\t"
ofs="\t"

ref_tuples = set()

with open(ref_vcf_filename) as f:
    for line in f:
        if line.startswith("#"):
            continue

        toks = line.split(fs)

        # CHROM, POS, REF, ALT
        ref_tuples.add( (toks[0], toks[1], toks[2], toks[3]) )

with open (cell_vcf_filename) as f:
    for line in f:
        if line.startswith("#"):
            sys.stdout.write(line)
            continue

        toks = line.split(fs)
        if (toks[0], toks[1], toks[2], toks[3]) in ref_tuples:
            continue

        sys.stdout.write(line)

