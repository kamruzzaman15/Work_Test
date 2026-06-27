#!/usr/bin/env python3
"""Part 2, Step 3: aggregate per-pool Kraken2/Bracken reports into one summary.

Answers: "what viruses are present, and how confident are you?"

Confidence here is NOT Kraken2's per-read confidence score (we disabled that
filter - see part2_step2_run_kraken2.sh for why it's miscalibrated for long
ONT reads). Instead we assign a simple, explainable per-call confidence tier
from two things any reader can sanity-check themselves:
  - n_reads: how many reads were assigned to the species (more = less likely
    to be a single misclassified read)
  - pct_of_pool: what fraction of the pool's classified reads this represents
    (a species that's 1 of 50,000 reads is less certain than 1 of 5)

Tiers (simple thresholds, stated plainly so they're easy to argue with):
  high   : n_reads >= 50  AND pct_of_pool >= 0.01 (>=1% of classified reads)
  medium : n_reads >= 10  AND pct_of_pool >= 0.001
  low    : everything else that was called at all
"""
import csv
import glob
import os

REPORTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results", "kraken2_reports")
OUT_CSV = os.path.join(os.path.dirname(__file__), "..", "results", "part2_virus_summary.csv")

# Common respiratory virus genera/species (by name substring) used later for
# the genome-coverage step - flagged here too so the summary table makes
# clear which hits are "common respiratory viruses" vs. other vertebrate
# viruses incidentally present in the pool.
COMMON_RESP_VIRUS_KEYWORDS = [
    "Influenza", "Orthopneumovirus", "Respirovirus", "Rhinovirus", "Enterovirus",
    "Metapneumovirus", "coronavirus", "Adenovirus", "Orthorubulavirus",
    "Mastadenovirus", "parainfluenza",
]


def is_common_respiratory(name: str) -> bool:
    return any(kw.lower() in name.lower() for kw in COMMON_RESP_VIRUS_KEYWORDS)


def confidence_tier(n_reads: int, pct_of_pool: float) -> str:
    if n_reads >= 50 and pct_of_pool >= 0.01:
        return "high"
    if n_reads >= 10 and pct_of_pool >= 0.001:
        return "medium"
    return "low"


def parse_kraken_report(path):
    """Parse Kraken2's standard report format, return species-level (S) rows."""
    rows = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                continue
            pct, n_clade, n_direct, rank, taxid, name = fields
            if rank == "S":
                rows.append({
                    "name": name.strip(),
                    "taxid": taxid,
                    "n_reads_clade": int(n_clade),
                })
    return rows


def main():
    records = []
    report_files = sorted(glob.glob(os.path.join(REPORTS_DIR, "*.kraken2_report.txt")))
    for report_path in report_files:
        pool = os.path.basename(report_path).replace(".kraken2_report.txt", "")
        species_rows = parse_kraken_report(report_path)
        total_reads_in_pool = sum_total_reads(report_path)
        total_classified = sum(r["n_reads_clade"] for r in species_rows)

        for row in species_rows:
            pct_of_pool = row["n_reads_clade"] / total_classified if total_classified else 0.0
            records.append({
                "pool": pool,
                "species": row["name"],
                "taxid": row["taxid"],
                "n_reads": row["n_reads_clade"],
                "pct_of_classified_reads": round(pct_of_pool * 100, 4),
                "pct_of_total_pool_reads": round(row["n_reads_clade"] / total_reads_in_pool * 100, 4) if total_reads_in_pool else 0.0,
                "confidence": confidence_tier(row["n_reads_clade"], pct_of_pool),
                "is_common_respiratory_virus": is_common_respiratory(row["name"]),
                "total_reads_in_pool": total_reads_in_pool,
                "total_classified_in_pool": total_classified,
            })

    records.sort(key=lambda r: (r["pool"], -r["n_reads"]))

    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)

    print(f"wrote {len(records)} species-level calls across {len(report_files)} pools -> {OUT_CSV}")


def sum_total_reads(report_path):
    """Root line ('unclassified' + root 'classified') gives total reads in the pool."""
    total = 0
    with open(report_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                continue
            rank = fields[3]
            if rank in ("U", "R"):
                total += int(fields[1])
    return total


if __name__ == "__main__":
    main()
