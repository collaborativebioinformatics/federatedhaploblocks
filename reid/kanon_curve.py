#!/usr/bin/env python3
"""k-anonymity suppression cost curve from per-block cluster assignment files.

Usage: kanon_curve.py KMAX FILE [FILE ...]

Each FILE is a two-column representative/member TSV, one row per haplotype, as
published under data.haploblocks.org/haploblock_hashes/1000G/<chr>/clusters/.
Cost is the fraction of haplotype-block observations lost when every cluster
smaller than k is suppressed.
"""
import sys
from collections import Counter

kmax = int(sys.argv[1])
hist, total, nblocks = Counter(), 0, 0

for path in sys.argv[2:]:
    sizes = Counter()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        sizes[f[0]] += 1
    for n in sizes.values():
        hist[n] += 1
        total += n
    nblocks += 1

nclust = sum(hist.values())
print("%d blocks, %d clusters, %d observations" % (nblocks, nclust, total))
print("singletons: %d (%.1f%% of clusters)\n" % (hist[1], 100.0 * hist[1] / nclust))
print("%-5s %11s %14s %9s" % ("k", "clusters", "observations", "cost"))
for k in range(2, kmax + 1):
    c = sum(v for sz, v in hist.items() if sz < k)
    o = sum(sz * v for sz, v in hist.items() if sz < k)
    print("%-5d %11d %14d %8.3f%%" % (k, c, o, 100.0 * o / total))
