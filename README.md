# Cluster assignment by BLAST

Testing whether a new sequence can be assigned to the right haploblock cluster
by BLASTing it against that cluster's representative, and what governs whether
the assignment is correct.

Context: cluster IDs come out of MMseqs2 and depend on the input set, so they
are not portable between sites. If sites could instead look up against a shared
set of cluster representatives, the IDs would mean the same thing everywhere.
This tests whether that lookup works.

## Method

For each locus:

1. `makeblastdb` on the cluster representatives
2. `blastn` each member sequence against them, `-evalue 1e-20 -dust no
   -qcov_hsp_perc 80`
3. Assign each member to its highest-scoring representative
4. Compare against ground truth from the `_members.tsv` files
5. Report accuracy and the score margin between the top and second hit, split
   by cluster class

Data is not in this repo. FASTA and TSV files come from
`https://data.haploblocks.org/bidirectional_blast_samples/`.

## Results

### Four loci, representatives as the anchor

| Locus | BIG | SMALL | Notes |
|---|---|---|---|
| chr1:10023878-10175888 | 0/40 | 4/5 | every BIG member assigned to the SMALL representative |
| chr2:100118009-100180909 | 40/40 | 1/5 | 4 of 5 SMALL went to SINGLE |
| chr3:100195058-100330759 | 40/40 | 15/15 | |
| chr1:100330740-100565735 | 97.75% | 100% | full population, see below |

chr1:10023878 is the pathological case. Representatives there are the least
separated of the sample: 25 mismatches between BIG and SMALL over 152 kb, about
0.16 per kb, against 1.14 per kb at chr2 and 0.80 at chr3.

### Size control, chr1:100330740-100565735

BIG, MEDIUM and SMALL share this locus, so cluster size is the only variable.
All 4,116 members scored, no sampling.

| Class | Members | Correct | Median margin |
|---|---|---|---|
| SMALL | 8 | 100.00% | 850 |
| BIG | 3,823 | 97.75% | 600 |
| MEDIUM | 285 | 84.21% | 400 |

Confusion matrix, rows are truth:

```
           BIG  MEDIUM  SMALL
BIG       3737      60     26
MEDIUM       0     240     45
SMALL        0       0      8
```

**Cluster size does not predict accuracy.** The largest cluster outperforms the
medium one. The margin between a member's own representative and the next best
one does predict it, and the ordering matches exactly.

Errors are directional. MEDIUM never lands on BIG in 285 attempts, while BIG
lands on MEDIUM 60 times. All 45 MEDIUM errors go to SMALL, its second choice.

Representative-to-representative distance does not predict where members land.
At this locus BIG-MEDIUM is 168 mismatches, MEDIUM-SMALL 188, BIG-SMALL 204, so
MEDIUM's nearest representative is BIG, yet none of its errors went there.

Margins throughout are around 0.1 percent of total score.

## Known limitations

**`make_consensus.py` is only valid where indels are rare.** It tallies the
majority base per position across sequences of the modal length, which assumes
equal length implies positional alignment. That is false at
chr1:100330740-100565735, where two sequences both 235,007 bp long differ at
8.4 percent of positions under direct comparison while BLAST puts them at 99.9
percent identity. Compensating indels shift the middle. Any consensus built from
that locus should be discarded. Building a real consensus needs multiple
sequence alignment, not a column tally.

**Strict bidirectional best hit does not fit this problem.** With a handful of
representatives and many queries, at most one query per representative can be
mutually best, so the reciprocal check fails for correctly assigned members by
construction. The `mutual` column in the run logs is not informative. A
reciprocal check would need to run against cluster members rather than the
representative.

**chr2 and chr3 used 40 queries per BIG cluster.** At chr1:100330740 the
40-query subset agreed closely with the full population, 97.5 against 97.75 for
BIG and 82.5 against 84.21 for MEDIUM, which suggests the sampled figures are
sound, but they have not been rescored in full.

## Scripts

| File | Purpose |
|---|---|
| `bbh.py` | main assignment test, one locus |
| `bbh_chr2.py`, `bbh_chr3.py` | same, pointed at those loci |
| `bbh_run2.py` | four-class version for the size control |
| `score_all.py` | scores a full-population BLAST output, confusion matrix and margins |
| `make_consensus.py` | majority-base consensus, see limitations |
| `separation.py` | pairwise identity between representatives within each haploblock |
| `cov.py` | HSP counts and query coverage per pair |
| `diag.py` | where queries land, relative margin between top and second hit |
| `tophsp.py` | compares summed-bitscore against best-single-HSP scoring |

`logs/` holds the raw output of each run.
