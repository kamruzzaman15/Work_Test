#!/usr/bin/env bash
# Part 2, Step 5: download the reference genomes picked in step 4 from NCBI
# and concatenate into one combined FASTA for read mapping (step 6).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PICKS_CSV="${SCRIPT_DIR}/../results/part2_reference_genome_picks.csv"
OUT_DIR="${SCRIPT_DIR}/../data/reference_genomes"
COMBINED="${OUT_DIR}/part2_combined_reference.fasta"

mkdir -p "$OUT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate biosurv

: > "$COMBINED"

# skip header line; one row per accession/segment
tail -n +2 "$PICKS_CSV" | while IFS=',' read -r species species_taxid total_reads leaf_taxid leaf_name accession segment_index n_segments note; do
  if [ -z "$accession" ]; then
    echo "skip (no accession): ${species}"
    continue
  fi
  out_fa="${OUT_DIR}/${accession}.fasta"
  if [ ! -f "$out_fa" ]; then
    echo "fetching ${accession} (${species}, segment ${segment_index}/${n_segments})"
    efetch -db nuccore -id "$accession" -format fasta > "$out_fa"
    sleep 0.4  # be polite to NCBI eutils
  fi
  # contig name "accession|species_name|segN" so coverage tables can group segments per virus.
  # FASTA/SAM sequence names are truncated at the first whitespace by minimap2/samtools,
  # so spaces in species names must become underscores or the name gets silently cut short.
  sp_nospace=$(echo "$species" | tr ' ' '_')
  seg_suffix=""
  if [ -n "$n_segments" ] && [ "$n_segments" -gt 1 ] 2>/dev/null; then
    seg_suffix="|seg${segment_index}"
  fi
  awk -v sp="$sp_nospace" -v acc="$accession" -v seg="$seg_suffix" \
    'NR==1{print ">"acc"|"sp seg; next} {print}' "$out_fa" >> "$COMBINED"
done

echo "combined reference: $COMBINED"
grep -c "^>" "$COMBINED"
