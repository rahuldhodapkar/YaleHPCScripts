#!/usr/bin/env python
#
# find all unique chrom, pos, ref, alt combos

from __future__ import print_function
import sys

header = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample"

calls = {}
n = 0

for line in sys.stdin:
    if line.startswith("#"):
        continue

    n = n + 1
    if n % 10000 == 0:
        print("reading line {}".format(n), file=sys.stderr)

    cols = line.split("\t")

    calls[(cols[0], cols[1], cols[3], cols[4])] = line

print("writing completed set", file=sys.stderr)

print(header)
for c in calls.keys():
    sys.stdout.write(calls[c])