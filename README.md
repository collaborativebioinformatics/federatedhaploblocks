# Cluster assignment by BLAST

Testing whether a sequence can be assigned to the right haploblock cluster by
BLASTing it against that cluster's representative, and what governs whether the
assignment is correct.

Context: cluster IDs come out of MMseqs2 and depend on the input set, so they are
not portable between sites. If sites could instead look up against a shared set
of cluster representatives, the IDs would mean the same thing everywhere. This
tests whether that lookup works.

## Method

For each locus:

1. `makeblastdb` on the cluster representatives.
2. `blastn` each member sequence against them, with
   `-evalue 1e-20 -dust no -qcov_hsp_perc 80`.
3. **Sum bitscores across all HSPs** for each query-representative pair, then
   assign the member to whichever representative has the highest total. This is a
   choice, not a default. Re-scoring with best-single-HSP instead of the sum gave
   identical assignments at chr1:10023878, so it does not appear to drive the
   results, but it has not been checked at every locus.
4. Compare against ground truth from the `_members.tsv` files.
5. Report accuracy, the margin between the top and second total, and the number
   of exact ties, split by cluster class.

Ties are currently resolved by iteration order, not by any rule. Where they occur
they are reported alongside the accuracy figures.

Data is not in this repo. FASTA and TSV files come from
`https://data.haploblocks.org/bidirectional_blast_samples/`.

## Results

### Four loci, representatives as the anchor

| Locus | BIG | SMALL |
|---|---|---|
| chr1:10023878-10175888 | 0/40 | 4/5 |
| chr2:100118009-100180909 | 40/40 | 1/5 |
| chr3:100195058-100330759 | 40/40 | 15/15 |
| chr1:100330740-100565735 | 97.75% | 100% |

chr1:10023878 behaves very differently from chr2 and chr3. Comparing the BIG and
SMALL representatives at each locus gives 0.16 mismatches per kb at chr1, 0.80 at
chr3 and 1.14 at chr2, so that particular pair sits about seven times closer
together at chr1.

Note this is specific to the BIG-SMALL pair. Other pairs are closer still: at
chr3 the SMALL and SINGLE representatives differ by a single mismatch across
135,678 bp.

At chr2, four of the five SMALL members were assigned to the SINGLE
representative.

### chr1:100330740-100565735, all cluster sizes at one locus

Every member scored, 4,116 in total, no sampling.

| Class | Members | Correct | Median margin | Exact ties |
|---|---|---|---|---|
| SMALL | 8 | 100.00% | 850 | 0 |
| BIG | 3,823 | 97.75% | 600 | 29 |
| MEDIUM | 285 | 84.21% | 400 | 12 |

Confusion matrix, rows are truth:

```
           BIG  MEDIUM  SMALL
BIG       3737      60     26
MEDIUM       0     240     45
SMALL        0       0      8
```

41 queries had their top two totals exactly equal. In every one the correct class
was among the tied options, so these are coin flips rather than systematic
errors, but the percentages should be read with that in mind. There were no
near-misses: everything else was separated by at least 50 bits.

**Accuracy is not monotonic in cluster size here.** BIG has 13 times the members
of MEDIUM and is substantially more accurate. This is a within-locus comparison
across cluster sizes rather than a controlled size experiment, since the clusters
also differ in composition and internal diversity, but it does rule out the
simple explanation that a representative cannot stand for a large cluster.

**The margin tracks accuracy in the same ordering**: SMALL 850, BIG 600, MEDIUM
400. Three class-level points at one locus, so an association rather than a
validated predictor, but it is the only quantity found that orders the same way.

**Errors are directional.** MEDIUM never lands on BIG in 285 attempts, while BIG
lands on MEDIUM 60 times. All 45 MEDIUM errors go to SMALL, which has the
second-highest median score for MEDIUM members.

**Representative geometry does not predict member assignment.** At this locus
BIG-MEDIUM is 168 mismatches, MEDIUM-SMALL 188 and BIG-SMALL 204 over 235 kb, so
MEDIUM's nearest representative is BIG, yet none of its errors went there.

Margins are roughly 0.09 to 0.20 percent of total score.

## Known limitations

**`make_consensus.py` is only valid where indels are rare.** It tallies the
majority base per position across sequences of the modal length, which assumes
equal length implies positional alignment. That is false at
chr1:100330740-100565735. Two members both exactly 235,007 bp,
`NA18940_..._hap0` and `HG04100_..._hap1`, differ at 8.42 percent of positions
compared index by index, but BLAST aligns them at 99.989 percent identity with 14
mismatches and 12 gaps. The gaps are the point: compensating indels leave the
lengths equal while shifting everything between them. Any consensus built from
that locus should be discarded, and building a real one needs multiple sequence
alignment rather than a column tally.

**Strict bidirectional best hit does not fit this problem.** With a handful of
representatives and many queries, at most one query per representative can be
mutually best, so the reciprocal check fails for correctly assigned members by
construction. The `mutual` column in the run logs is not informative. A
reciprocal check would need to run against cluster members rather than the
representative.

**Ties are resolved arbitrarily.** `score_all.py` sorts by total and takes the
first, so an exact tie is decided by insertion order. Tie counts are reported but
the tied cases are still counted as correct or incorrect by that arbitrary
resolution.

**chr2 and chr3 used 40 queries per BIG cluster.** At chr1:100330740 the
40-member subset reproduced the full population closely, 97.5 against 97.75 for
BIG and 82.5 against 84.21 for MEDIUM. That shows subsampling worked at that
locus, but does not establish that the chr2 and chr3 subsets are representative
of their own loci.

**One locus per conclusion.** The size comparison rests on a single locus and the
separation comparison on three. Nothing here has been tested across the roughly
39,000 haploblocks genome-wide.

## Scripts

| File | Purpose |
|---|---|
| `bbh.py` | main assignment test, one locus |
| `bbh_chr2.py`, `bbh_chr3.py` | same, pointed at those loci |
| `bbh_run2.py` | variant with MEDIUM support, for the size-control locus |
| `score_all.py` | scores a full-population BLAST output, confusion matrix and margins |
| `make_consensus.py` | majority-base consensus, see limitations |
| `separation.py` | pairwise identity between representatives within each haploblock |
| `cov.py` | HSP counts and query coverage per pair |
| `diag.py` | where queries land, relative margin between top and second hit |
| `tophsp.py` | compares summed-bitscore against best-single-HSP scoring |

`logs/` holds the raw output of each run.
