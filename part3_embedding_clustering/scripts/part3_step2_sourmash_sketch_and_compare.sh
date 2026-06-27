#!/usr/bin/env bash
# Part 3, Step 2: sketch every selected read individually with sourmash, then
# build the all-by-all pairwise similarity matrix.
#
# We use a fixed-size (num=500) MinHash sketch rather than the usual
# scaled/genome-sized sketch: these are short (~3kb) individual reads, not
# whole genomes, and a scaled sketch would sample too few k-mers per read to
# compare reliably at this length. k=21 is sourmash's standard species-level
# k-mer size.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FASTA="${SCRIPT_DIR}/../data/selected_reads/part3_selected_reads.fasta"
SIG_DIR="${SCRIPT_DIR}/../data/sourmash_sigs"
SIGS_OUT="${SIG_DIR}/part3_reads.sig.zip"
COMPARE_OUT="${SCRIPT_DIR}/../results/part3_similarity_matrix"

mkdir -p "$SIG_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biosurv

echo "sketching reads individually..."
sourmash sketch dna -p k=21,num=500,noabund --singleton \
  -o "$SIGS_OUT" --force \
  "$FASTA"

# Plain Jaccard similarity (not containment - containment estimation needs
# scaled sketches; our num=500 sketches are sized for short individual reads).
echo "building all-by-all similarity matrix..."
sourmash compare --ksize 21 \
  -o "${COMPARE_OUT}.cmp" \
  --csv "${COMPARE_OUT}.csv" \
  "$SIGS_OUT"

echo "done. similarity matrix -> ${COMPARE_OUT}.csv"
