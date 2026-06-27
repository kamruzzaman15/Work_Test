#!/usr/bin/env bash
# Part 2, Step 6: map every pool's reads against the combined reference of
# detected common respiratory viruses, and compute genome coverage stats.
#
# We map ALL reads (not just the ones Kraken2 classified) - minimap2 alignment
# identity is a more direct, independent check on Kraken2's k-mer calls, and
# this is also simply how you measure coverage: you need read-to-genome
# alignments, not just a taxonomic label.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data/raw_pools"
REF="${SCRIPT_DIR}/../data/reference_genomes/part2_combined_reference.fasta"
OUT_DIR="${SCRIPT_DIR}/../results/coverage"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biosurv

for fq in "${DATA_DIR}"/*.respiratory.fasta.gz; do
  pool=$(basename "$fq" .respiratory.fasta.gz)
  echo "=== mapping pool: ${pool} ==="
  bam="${OUT_DIR}/${pool}.bam"

  minimap2 -ax map-ont -t "$THREADS" "$REF" "$fq" 2>"${OUT_DIR}/${pool}.minimap2.log" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"

  samtools coverage "$bam" > "${OUT_DIR}/${pool}.coverage_summary.tsv"
  samtools depth -a "$bam" > "${OUT_DIR}/${pool}.depth.tsv"
done

echo "done. per-pool coverage summaries and per-base depth in ${OUT_DIR}"
