# Cluster assignment by BLAST

Testing whether a sequence can be assigned to the correct haploblock cluster by BLASTing it against a shared set of cluster representatives, and what factors are associated with correct assignment.

Cluster IDs produced by MMseqs2 depend on the input set, so the numeric or internal cluster labels are not portable between sites. If different sites could instead assign sequences against a shared panel of representative haplotypes, the resulting labels could be made consistent across datasets. This experiment tests how well that representative-based lookup works.

## Method

For each locus:

1. Run `makeblastdb` on the cluster representatives.

2. Run `blastn` for each member sequence against the representative database using:

   ```bash
   -evalue 1e-20 -dust no -qcov_hsp_perc 80 -max_target_seqs 20
   ```

3. Sum BLAST bitscores across all reported HSPs for each query-representative pair.

4. For each query, identify the representative or representatives with the highest total bitscore.

5. If exactly one representative has the highest total, treat the assignment as unambiguous.

6. If two or more representatives share the same highest total, report the query as an exact tie rather than resolving it arbitrarily.

7. Compare unambiguous assignments against ground truth from the corresponding `_members.tsv` files.

8. Report:

   * unique correct assignments,
   * unique incorrect assignments,
   * exact ties,
   * accuracy bounds under alternative tie treatment,
   * and the margin between the highest and second-highest totals for unambiguous assignments.

Summing HSP bitscores is a scoring choice rather than a BLAST default. See the limitations section below.

Data is not stored in this repository. FASTA and TSV inputs are available from:

`https://data.haploblocks.org/bidirectional_blast_samples/`

## Results

### Four loci

| Locus                    |          BIG |       MEDIUM | SMALL |
| ------------------------ | -----------: | -----------: | ----: |
| chr1:10023878-10175888   |         0/40 |            — |   4/5 |
| chr2:100118009-100180909 |        40/40 |            — |   1/5 |
| chr3:100195058-100330759 |        40/40 |            — | 15/15 |
| chr1:100330740-100565735 | 97.28–98.04% | 82.46–86.67% |  100% |

The first three loci are the original sampled runs and used 40 queries per BIG cluster. They have not yet been rescored with explicit tie handling.

The fourth locus was rescored across the full population with exact ties separated from unambiguous assignments. Its reported ranges give the lower and upper accuracy bounds obtained by counting all tied queries as incorrect or correct, respectively.

### Representative separation across the first three loci

chr1:10023878-10175888 behaves very differently from chr2 and chr3.

For the BIG-SMALL representative pair:

* chr1:10023878-10175888: approximately **0.16 mismatches per kb**
* chr3:100195058-100330759: approximately **0.80 mismatches per kb**
* chr2:100118009-100180909: approximately **1.14 mismatches per kb**

The BIG and SMALL representatives are therefore about seven times closer at the pathological chr1 locus than at chr2.

This comparison is specific to the BIG-SMALL pair. Other representative pairs can be closer still: at chr3, the SMALL and SINGLE representatives differ by only one mismatch across 135,678 aligned bases.

At chr2, four of the five sampled SMALL members were assigned to the SINGLE representative.

## chr1:100330740-100565735

### Full-population assignment

This locus contains BIG, MEDIUM, and SMALL clusters at the same genomic region.

Every member was scored:

**4,116 total sequences, no sampling.**

| Class  | Members | Unique correct | Unique wrong | Exact ties | Median margin |
| ------ | ------: | -------------: | -----------: | ---------: | ------------: |
| SMALL  |       8 |              8 |            0 |          0 |           850 |
| BIG    |   3,823 |          3,719 |           75 |         29 |           600 |
| MEDIUM |     285 |            235 |           38 |         12 |           400 |

The 41 tied queries are kept separate because the representative scores are exactly equal.

In every one of those 41 ties, the ground-truth class is among the tied highest-scoring representatives.

Therefore, the ties are ambiguous under the scoring rule rather than cases where the correct representative scores lower.

There were no nonzero near-ties below 50 bits: every unambiguous assignment was separated from the runner-up by at least 50 bits.

### Accuracy bounds

Because the tied queries cannot be assigned uniquely under the current scoring rule, the true observed accuracy lies between two bounds:

* lower bound: treat every tie as incorrect,
* upper bound: treat every tie as correct.

Accuracy among unambiguous queries is also reported separately.

| Class  | Ties counted wrong | Ties counted right | Unambiguous only |
| ------ | -----------------: | -----------------: | ---------------: |
| SMALL  |            100.00% |            100.00% |          100.00% |
| BIG    |             97.28% |             98.04% |           98.02% |
| MEDIUM |             82.46% |             86.67% |           86.08% |

### Unambiguous confusion matrix

Rows are ground-truth classes and columns are predicted classes:

```text
           BIG  MEDIUM  SMALL
BIG       3719      56     19
MEDIUM       0     235     38
SMALL        0       0      8
```

This matrix contains only the 4,075 queries with a unique highest-scoring representative.

### Cluster size

**Accuracy is not monotonic in cluster size at this locus.**

BIG contains 3,823 members, approximately 13 times as many as MEDIUM, yet its assignment accuracy is substantially higher.

The BIG and MEDIUM accuracy ranges do not overlap under either extreme treatment of ties.

This is a within-locus comparison across cluster sizes rather than a controlled experiment in which cluster size is literally the only changing variable. The clusters also differ in sequence composition, internal diversity, and representative geometry.

The result therefore does not establish a general causal relationship between cluster size and assignment quality.

It does, however, show that **large cluster size by itself is not sufficient to produce poor representative-based assignment at this locus**.

### Assignment margin

The median assignment margin follows the same class ordering as accuracy:

* SMALL: **850**
* BIG: **600**
* MEDIUM: **400**

SMALL has the largest median margin and highest observed accuracy, BIG is intermediate, and MEDIUM has the smallest median margin and lowest accuracy.

This is a descriptive association based on three class-level observations at one locus, not a fitted or externally validated predictor.

Margins are small relative to the total BLAST scores, roughly **0.09% to 0.20%** of totals near 433,000 bits.

### Direction of errors

Errors are strongly directional.

Among unambiguous assignments:

* MEDIUM never assigns to BIG across 273 queries.
* BIG assigns to MEDIUM 56 times.
* All 38 unambiguous MEDIUM errors assign to SMALL.
* SMALL has no unambiguous errors.

For MEDIUM members, SMALL also has the second-highest median representative score after MEDIUM itself.

This indicates that the MEDIUM-SMALL boundary is substantially more problematic than the MEDIUM-BIG boundary.

### Representative geometry

**Pairwise representative distance does not explain the direction of errors at this locus.**

Representative-to-representative comparisons give:

| Pair            | Identity | Aligned bases | Mismatches | Gaps |
| --------------- | -------: | ------------: | ---------: | ---: |
| BIG vs MEDIUM   |  99.902% |       235,065 |        168 |   63 |
| BIG vs SMALL    |  99.892% |       235,059 |        204 |   51 |
| MEDIUM vs SMALL |  99.886% |       235,069 |        188 |   80 |

By mismatch count, BIG is MEDIUM's nearest representative.

Nevertheless, none of the unambiguous MEDIUM errors assign to BIG; all 38 assign to SMALL.

Pairwise representative distance alone therefore does not explain error direction at this locus.

## Known limitations

### `make_consensus.py` assumes positional correspondence between sequences

`make_consensus.py` constructs a simple majority-base consensus by selecting sequences of the modal length and tallying the nucleotide at each raw sequence index.

That assumes that the same sequence index corresponds to the same homologous position across members.

Equal sequence length does not guarantee this when compensating insertions and deletions are present.

At chr1:100330740-100565735, two BIG-cluster members,

* `NA18940_chr1_region_100330740-100565735_hap0`
* `HG04100_chr1_region_100330740-100565735_hap1`

are both exactly **235,007 bp** long.

A direct index-by-index comparison produces:

* **19,797 positional mismatches**
* **8.424% apparent disagreement**

Yet BLAST aligns the same pair across 235,013 bases at:

* **99.989% identity**
* **14 mismatches**
* **12 gaps**

The gaps are the important part: compensating indels can preserve overall sequence length while shifting the correspondence of bases in the interior of the sequence.

Therefore, equal-length sequences cannot safely be combined by raw positional voting at this locus.

Consensus sequences produced by `make_consensus.py` for this region should be discarded.

A proper consensus requires an alignment that establishes homologous positions before column-wise voting.

### Summing HSP bitscores is a scoring choice

The assignment procedure sums bitscores across all BLAST HSPs belonging to each query-representative pair.

This is not equivalent to selecting the single best HSP.

Multiple overlapping or alternative local alignments may contribute to the same total and can potentially overweight repetitive regions.

At chr1:10023878-10175888, re-scoring with the best single HSP produced the same assignments as summed-HSP scoring.

That comparison has not yet been repeated across the other loci.

### Exact ties are genuinely ambiguous under the current score

At chr1:100330740-100565735, 41 of 4,116 queries produced exactly equal highest total bitscores for two or more representatives.

The correct class was among the tied winners in all 41 cases.

These queries therefore cannot be uniquely classified using the current summed-bitscore rule alone.

The updated scoring procedure reports them separately rather than resolving them by insertion or iteration order.

### Strict reciprocal best hit does not fit this assignment problem

A strict bidirectional-best-hit test is not informative in this many-to-one setting.

There are many member queries but only a small number of representatives.

For a representative to be a strict reciprocal best hit, only one member can be its best reverse match.

Most correctly assigned members therefore cannot satisfy reciprocal best-hit status by construction.

The `mutual` column in earlier run logs should not be interpreted as an assignment-quality metric.

A reciprocal test would need a different design, such as comparison against cluster members rather than only cluster representatives.

### Some loci were evaluated by sampling

chr2 and chr3 used 40 queries from each BIG cluster rather than scoring their complete BIG populations.

At chr1:100330740-100565735, the original 40-member BIG subset produced 97.5% accuracy under the earlier order-dependent scorer, compared with 97.75% for the full population under that same scoring approach.

This shows that the BIG sample happened to approximate the full-population result at that locus.

It does not establish that the chr2 and chr3 samples are similarly representative of their complete populations.

### Limited number of loci

The current experiments cover four loci.

The cluster-size comparison is based on one locus.

The representative-separation comparison is based primarily on three loci.

These results should therefore be treated as observations from the tested regions rather than general rules across the genome.

Approximately 39,000 haploblocks exist in the broader dataset, and the current analysis has not yet been extended across that full set.

## Scripts

| File                         | Purpose                                                                                                                          |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `bbh.py`                     | main representative-based assignment test for one locus                                                                          |
| `bbh_chr2.py`, `bbh_chr3.py` | locus-specific variants for chr2 and chr3                                                                                        |
| `bbh_run2.py`                | assignment variant with MEDIUM-class support for chr1:100330740-100565735                                                        |
| `score_all.py`               | scores a full-population BLAST output and reports unique assignments, exact ties, accuracy bounds, confusion matrix, and margins |
| `checks.py`                  | counts exact ties, nonzero near-ties, and clearly separated queries                                                              |
| `make_consensus.py`          | simple modal-length majority-base consensus; see limitations                                                                     |
| `separation.py`              | compares representative sequences within each locus                                                                              |
| `cov.py`                     | examines HSP counts and query coverage per query-representative pair                                                             |
| `diag.py`                    | examines assignment destinations and margins between competing representatives                                                   |
| `tophsp.py`                  | compares summed-HSP bitscore scoring against best-single-HSP scoring                                                             |

## Logs

`logs/` contains retained raw output from the experimental runs.

The current logs include the initial chr1 experiments plus the chr2 and chr3 sampled runs. Full-population tie-aware scoring for chr1:100330740-100565735 is produced by `score_all.py`.
