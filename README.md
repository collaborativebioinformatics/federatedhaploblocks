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
3. Sum bitscores across all HSPs for each query-representative pair, then assign
   the member to whichever representative has the highest total. This is a
   scoring choice, not a default. See limitations.
4. Compare against ground truth from the `_members.tsv` files.
5. Report unambiguous assignments and exact ties separately, plus the margin
   between the top total and the next distinct total.

A query where two or more representatives achieve the same highest total is
ambiguous under this scoring rule. `score_all.py` reports these separately rather
than resolving them, so no figure below depends on arbitrary tie-breaking.

Data is not in this repo. FASTA and TSV files come from
`https://data.haploblocks.org/bidirectional_blast_samples/`.

## Results

### Four loci

| Locus | BIG | MEDIUM | SMALL |
|---|---|---|---|
| chr1:10023878-10175888 | 0/40 | — | 4/5 |
| chr2:100118009-100180909 | 40/40 | — | 1/5 |
| chr3:100195058-100330759 | 40/40 | — | 15/15 |
| chr1:100330740-100565735 | 97.3–98.0% | 82.5–86.7% | 100% |

The first three loci used 40 queries per BIG cluster. The fourth scored every
member; ranges are the accuracy bounds explained below.

chr1:10023878 behaves very differently from chr2 and chr3. Comparing the BIG and
SMALL representatives at each locus gives 0.16 mismatches per kb at chr1, 0.80 at
chr3 and 1.14 at chr2, so that particular pair sits about seven times closer
together at chr1.

This is specific to the BIG-SMALL pair. Other pairs are closer still: at chr3 the
SMALL and SINGLE representatives differ by a single mismatch across 135,678 bp.

At chr2, four of the five SMALL members were assigned to the SINGLE
representative.

### chr1:100330740-100565735, all cluster sizes at one locus

Every member scored, 4,116 in total, no sampling.

| Class | Members | Unique correct | Unique wrong | Ties | Median margin |
|---|---|---|---|---|---|
| SMALL | 8 | 8 | 0 | 0 | 850 |
| BIG | 3,823 | 3,719 | 75 | 29 | 600 |
| MEDIUM | 285 | 235 | 38 | 12 | 400 |

Accuracy depends on how ties are treated, so both bounds are given:

| Class | Ties counted wrong | Ties counted right | Unambiguous only |
|---|---|---|---|
| SMALL | 100.00% | 100.00% | 100.00% |
| BIG | 97.28% | 98.04% | 98.02% |
| MEDIUM | 82.46% | 86.67% | 86.08% |

In all 41 ties the correct class was among the tied representatives, so these are
ambiguous assignments under the current scoring rule rather than cases where the
correct representative scores lower. There were no near-misses: every
unambiguous assignment was separated by at least 50 bits.

Confusion matrix, unambiguous assignments only, rows are truth:

```
           BIG  MEDIUM  SMALL
BIG       3719      56     19
MEDIUM       0     235     38
SMALL        0       0      8
```

**Accuracy is not monotonic in cluster size here.** BIG has 13 times the members
of MEDIUM and is more accurate, and the two ranges do not overlap under any tie
treatment. This is a within-locus comparison across cluster sizes rather than a
controlled size experiment, since the clusters also differ in composition and
internal diversity, but it does rule out the simple explanation that a
representative cannot stand for a large cluster.

**The margin tracks accuracy in the same ordering**: SMALL 850, BIG 600, MEDIUM
400. Three class-level points at one locus, so an association rather than a
validated predictor, but it is the only quantity found that orders the same way.

**Errors are directional.** MEDIUM never lands on BIG across 273 unambiguous
assignments, while BIG lands on MEDIUM 56 times. All 38 unambiguous MEDIUM errors
go to SMALL, which has the second-highest median score for MEDIUM members.

**Pairwise representative distance does not explain error direction at this
locus.** BIG-MEDIUM differs by 168 mismatches, MEDIUM-SMALL by 188 and BIG-SMALL
by 204 over 235 kb, so BIG is MEDIUM's nearest representative by this measure,
yet none of the MEDIUM errors go to BIG.

Margins are roughly 0.09 to 0.20 percent of total score.

## Known limitations

**`make_consensus.py` assumes positional correspondence between sequences.** It
tallies the majority base at each raw sequence index among sequences of the modal
length. Equal length does not guarantee positional alignment when compensating
indels are present. Two members of the BIG cluster at chr1:100330740-100565735,
`NA18940_..._hap0` and `HG04100_..._hap1`, are both exactly 235,007 bp and differ
at 8.42 percent of positions compared index by index, yet BLAST aligns them at
99.989 percent identity with 14 mismatches and 12 gaps. The gaps are the point.
Any consensus built from that locus should be discarded, and building a real one
requires an alignment that establishes homologous positions rather than a raw
column tally.

**Summing HSP bitscores is a scoring choice.** Multiple overlapping or
alternative local alignments may contribute to the same query-representative
total, which can overweight repetitive regions. At chr1:10023878, best-single-HSP
scoring produced identical assignments to summed scoring, but that comparison has
not been repeated at the other loci.

**Strict bidirectional best hit does not fit this problem.** With a handful of
representatives and many queries, at most one query per representative can be
mutually best, so the reciprocal check fails for correctly assigned members by
construction. The `mutual` column in the run logs is not informative. A
reciprocal check would need to run against cluster members rather than the
representative.

**chr2 and chr3 used 40 queries per BIG cluster.** At chr1:100330740 the
40-member subset reproduced the full population closely, 97.5 against 97.75 under
the earlier order-dependent scorer. That shows subsampling worked at that locus,
but does not establish that the chr2 and chr3 subsets are representative of their
own loci.

**One locus per conclusion.** The size comparison rests on a single locus and the
separation comparison on three. Nothing here has been tested across the roughly
39,000 haploblocks genome-wide.

## Scripts

| File | Purpose |
|---|---|
| `bbh.py` | main assignment test, one locus |
| `bbh_chr2.py`, `bbh_chr3.py` | same, pointed at those loci |
| `bbh_run2.py` | variant with MEDIUM support, for the size-control locus |
| `score_all.py` | scores a full-population BLAST output; reports unambiguous assignments, ties and accuracy bounds |
| `checks.py` | counts exact ties and near-misses in a scored run |
| `make_consensus.py` | majority-base consensus, see limitations |
| `separation.py` | pairwise identity between representatives within each haploblock |
| `cov.py` | HSP counts and query coverage per pair |
| `diag.py` | where queries land, relative margin between top and second hit |
| `tophsp.py` | compares summed-bitscore against best-single-HSP scoring |

`logs/` holds the raw output of each run.
