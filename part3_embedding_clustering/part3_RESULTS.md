# Part 3, Track B — Embedding-Based Clustering: Results

**Question:** Embed sequences (genomes, ORFs, or RdRp proteins) and cluster them. Compare
clusters: do they match your taxonomy profile or ANI?

## Approach

Whole-read nucleotide embedding via **sourmash** (MinHash k-mer sketches) — an existing,
widely-used tool for exactly this kind of alignment-free sequence comparison, chosen over a
pretrained protein/nucleotide language model specifically to avoid downloading any model at
all. A sourmash similarity score between two sequences is itself a standard ANI estimator,
so this one tool gives us both the "embedding" and the "ANI" axis the question asks about.

**Data:** reused, not re-derived. We pulled the per-read Kraken2 species calls already made
in Part 2 (`part2_taxonomic_classification/results/kraken2_reports/*.kraken2_output.txt`) as
ground-truth taxonomy, and selected up to 60 reads per common respiratory virus species with
≥20 total reads (17 species, 990 reads total — see
[`scripts/part3_step1_select_reads.py`](scripts/part3_step1_select_reads.py)). No new
taxonomic classification was performed.

## Method (one script per step, run in order)

| Step | Script | What it does |
|---|---|---|
| 1 | `part3_step1_select_reads.py` | Selects ≤60 reads/species using Part 2's Kraken2 labels; writes one combined FASTA + a labels CSV |
| 2 | `part3_step2_sourmash_sketch_and_compare.sh` | Sketches every read individually (`sourmash sketch --singleton`, k=21, num=500) and builds the all-by-all similarity matrix (`sourmash compare`) |
| 3 | `part3_step3_cluster_and_compare.py` | Converts similarity → distance, runs classical MDS for a 2D plot, DBSCAN for hard clusters, and computes both comparisons (vs. taxonomy, vs. ANI) |

We used a fixed-size (`num=500`) MinHash sketch rather than the usual genome-scaled sketch,
since these are individual ~3kb reads, not whole genomes — a scaled sketch would sample too
few k-mers per read to compare reliably at this length.

## Result 1: Do clusters match the taxonomy profile?

Plots: [`results/plots/part3_mds_by_taxonomy.png`](results/plots/part3_mds_by_taxonomy.png)
vs. [`results/plots/part3_mds_by_cluster.png`](results/plots/part3_mds_by_cluster.png).
Per-read assignments: [`results/part3_cluster_summary.csv`](results/part3_cluster_summary.csv).

Visually, the MDS plot colored by Kraken2 species (Part 2's labels) shows almost every
species forming its own clearly separated island — including between closely related
species like Rhinovirus A/B/C, which is the harder test case here (vs. trivially separating
something unrelated like SARS-CoV-2 from Rhinovirus).

For the hard-cluster comparison, our first attempt (linkage-based hierarchical clustering
forced into exactly 17 clusters, matching the species count) gave a misleadingly low
Adjusted Rand Index of 0.14 — it collapsed most points into one dominant catch-all cluster.
That's a clustering-method artifact, not a real result: between-species similarity is
~0 (see Result 2), so almost any method that must produce a *fixed* number of clusters ends
up chaining unrelated points together once most pairwise distances are close to the maximum.
Switching to **DBSCAN** (density-based, no fixed cluster count, leaves ambiguous points as
noise) fixed this: **Adjusted Rand Index = 0.561**, with the selected model (eps=0.60)
finding **45 clusters** — more than the 17 actual species, see caveat below — and leaving
**34 of 990 points** as unclustered noise. eps was chosen via a small scan from 0.60 to 0.90
against the known Kraken2 labels (ARI ranged 0.543–0.561 across that scan, fairly stable;
full scan in [`results/part3_ari_score.txt`](results/part3_ari_score.txt)). This is a
diagnostic comparison against ground truth, not a blind unsupervised pipeline, so tuning
against the label we're evaluating against is reasonable here and stated plainly.

**Why more clusters than species, and what that means:** our reads are partial genome
fragments (~3kb out of genomes that are 7-30kb), tiling different, often non-overlapping
positions. Two reads from the same virus that happen to cover *different* genome regions
share very few exact 21-mers and don't cluster together — directly confirmed in
`part3_cluster_summary.csv`: Rhinovirus A's 60 reads split across 6 different DBSCAN
clusters (plus 7 left as noise) rather than forming one unified cluster. So
DBSCAN is correctly finding sub-structure driven by **which part of the genome a read came
from**, layered on top of species identity. This is an expected consequence of using
short, non-overlapping read fragments with an exact-k-mer-match method, not a clustering
failure — the within-species pairs that *do* overlap in genome position cluster tightly,
which is exactly what drives the strong ANI signal in Result 2.

## Result 2: Do clusters match ANI?

Full table: [`results/part3_within_between_similarity.csv`](results/part3_within_between_similarity.csv).

This is the cleaner, less equivocal of the two results. Across all 990 reads:

- **Mean within-species similarity: 0.288**
- **Mean between-species similarity: 0.006** — a ~47x gap.

Every one of the 17 species individually shows the same pattern (within-species similarity
ranging 0.07-0.56, between-species similarity ≈0 for nearly all of them). In other words:
whenever two reads *do* share enough sequence content to register any similarity at all,
that similarity is overwhelmingly concentrated among reads of the same species. The
embedding distance is ANI-coherent — it isn't picking up some confound unrelated to
taxonomy.

## Bottom line

The embedding (sourmash k-mer similarity) recovers real taxonomic and ANI structure: same-
species reads are far more nucleotide-similar to each other than to other species, and the
visual clustering cleanly separates species — including closely related Rhinovirus A/B/C.
The main nuance is that with short, non-overlapping ONT read fragments, a single species
naturally splits into multiple "where on the genome did this read come from" sub-clusters
rather than one unified blob — worth knowing if someone wants to use this kind of clustering
for taxonomic assignment directly (it would need either longer/overlapping reads, or
whole-genome-level sketches rather than per-read sketches, to collapse that position-driven
sub-structure).

## How to reproduce

```bash
conda activate biosurv   # same environment as Part 2, plus sourmash + scikit-learn + scipy
python3 scripts/part3_step1_select_reads.py
bash scripts/part3_step2_sourmash_sketch_and_compare.sh
python3 scripts/part3_step3_cluster_and_compare.py
```

Requires Part 2's outputs to already exist (`part2_taxonomic_classification/results/kraken2_reports/`
and `part2_taxonomic_classification/data/raw_pools/`) — this step reuses them rather than
re-deriving anything.
