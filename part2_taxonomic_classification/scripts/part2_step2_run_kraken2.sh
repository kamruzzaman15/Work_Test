#!/usr/bin/env bash
# Part 2, Step 2: Taxonomic classification of each swab pool with Kraken2,
# followed by Bracken re-estimation of species-level abundance.
#
# Database: prebuilt Kraken2 "Viral" RefSeq index (genome-idx.s3.amazonaws.com),
# chosen so we don't have to build or download a large general-purpose database.
#
# NOTE on --confidence: Kraken2's confidence score is (k-mers matching the
# called taxon) / (total k-mers in the read). That ratio is well-calibrated for
# short (~150bp) reads but misleading for our ~3.2kb ONT reads: a real, strong
# viral hit can still have most of its k-mers fail to match (sequencing error,
# regions outside the conserved/reference-covered part of the genome), so a
# 10% threshold suppressed ~99% of true hits in testing (0.23% vs 49.6%
# classified on the same pool with/without it). We use the Kraken2 default
# (no confidence filter) and instead gauge confidence post-hoc from read
# counts and minimap2 mapping identity - see part2_RESULTS.md.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data/raw_pools"
DB_DIR="${SCRIPT_DIR}/../data/kraken2_db"
OUT_DIR="${SCRIPT_DIR}/../results/kraken2_reports"
THREADS="${THREADS:-8}"
# Bracken ships precomputed k-mer distributions only up to 300bp, but our ONT
# reads average ~3.2kb (seqkit stats). We use the largest available (300) as
# the closest approximation - see the long-read caveat in part2_RESULTS.md.
READ_LEN="${READ_LEN:-300}"

mkdir -p "$OUT_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biosurv

for fq in "${DATA_DIR}"/*.respiratory.fasta.gz; do
  pool=$(basename "$fq" .respiratory.fasta.gz)
  echo "=== classifying pool: ${pool} ==="

  kraken2 \
    --db "$DB_DIR" \
    --threads "$THREADS" \
    --gzip-compressed \
    --report "${OUT_DIR}/${pool}.kraken2_report.txt" \
    --output "${OUT_DIR}/${pool}.kraken2_output.txt" \
    "$fq"

  bracken \
    -d "$DB_DIR" \
    -i "${OUT_DIR}/${pool}.kraken2_report.txt" \
    -o "${OUT_DIR}/${pool}.bracken_species.txt" \
    -w "${OUT_DIR}/${pool}.bracken_report.txt" \
    -r "$READ_LEN" -l S \
    || echo "WARN: bracken failed for ${pool} (likely too few classified reads), keeping kraken2-only results"
done

echo "done. reports in ${OUT_DIR}"
