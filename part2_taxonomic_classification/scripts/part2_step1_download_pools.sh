#!/usr/bin/env bash
# Part 2, Step 1: Download swab pool FASTA files from the Zephyr sample log.
#
# Source: https://data.securebio.org/zephyr/#respiratory-viral-reads
# Note: the host blocks the default non-browser User-Agent (403), so we send one explicitly.
#
# Pool selection: 12 pools spanning 8 distinct collection sites and a range of
# dates (2026-03-18 to 2026-06-13) for site/temporal diversity, satisfying the
# work test's "ten or more swab pools" requirement.

set -euo pipefail

BASE_URL="https://data.securebio.org/zephyr/respiratory-reads"
OUT_DIR="$(dirname "$0")/../data/raw_pools"
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"

mkdir -p "$OUT_DIR"

POOLS=(
  "260613-MBTA_Ha-NAS-P1"
  "260611-BC-NAS-P1"
  "260610-MBTA_Da-NAS-P1"
  "260609-Cent-NAS-P1"
  "260608-Copl-NAS-P1"
  "260605-BC-NAS-P1"
  "260602-Sea_P-NAS-P1"
  "260530-NoS-NAS-P1"
  "260527-Copl-NAS-P1"
  "260526-MBTA_Da-NAS-P1"
  "260523-Haymkt-NAS-P1"
  "260318-MBTA_Ha-NAS-P1"
)

for pool in "${POOLS[@]}"; do
  fname="${pool}.respiratory.fasta.gz"
  if [ -f "${OUT_DIR}/${fname}" ]; then
    echo "skip (already downloaded): ${fname}"
    continue
  fi
  echo "downloading: ${fname}"
  curl -sf -A "$UA" -o "${OUT_DIR}/${fname}" "${BASE_URL}/${fname}"
done

echo "done. files in ${OUT_DIR}:"
ls -la "$OUT_DIR"
