#!/usr/bin/env python3
"""Part 2, Step 4: pick one reference genome accession per detected common
respiratory virus, for the genome-coverage analysis (step 5).

Kraken2's report aggregates reads at the species rank (e.g. "Rhinovirus B"),
but each individual read is actually assigned to a specific leaf taxid in the
NCBI taxonomy (a particular strain/genotype with its own RefSeq accession).
To map reads for coverage, we need one concrete reference sequence per virus,
so here we:
  1. Walk every classified read's leaf taxid up the taxonomy tree to find
     which of our target species it belongs to.
  2. For each species, take the single most-frequent leaf taxid across all
     12 pools combined (i.e. the strain our reads actually look most like).
  3. Look up that leaf taxid's RefSeq accession in the Kraken2 DB's
     seqid2taxid.map (the exact sequence Kraken2 itself used).

This deliberately picks "the genotype most consistent with our reads" rather
than an arbitrary representative strain - important for diverse genera like
Rhinovirus/Enterovirus where strain choice materially affects mapping
identity and coverage breadth.
"""
import csv
import glob
import os
import re
from collections import Counter, defaultdict

BASE = os.path.join(os.path.dirname(__file__), "..")
DB_DIR = os.path.join(BASE, "data", "kraken2_db")
REPORTS_DIR = os.path.join(BASE, "results", "kraken2_reports")
SUMMARY_CSV = os.path.join(BASE, "results", "part2_virus_summary.csv")
OUT_CSV = os.path.join(BASE, "results", "part2_reference_genome_picks.csv")

MIN_TOTAL_READS = 20  # ignore species too sparse for a meaningful coverage plot


def load_taxonomy():
    parent_of = {}
    name_of = {}
    with open(os.path.join(DB_DIR, "ktaxonomy.tsv")) as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split("|")]
            taxid, parent, name = parts[0], parts[1], parts[4]
            parent_of[taxid] = parent
            name_of[taxid] = name
    return parent_of, name_of


def load_seqid2taxid():
    taxid_to_accessions = defaultdict(list)
    with open(os.path.join(DB_DIR, "seqid2taxid.map")) as fh:
        for line in fh:
            seqid, taxid = line.rstrip("\n").split("\t")
            m = re.search(r"\|([A-Z]{1,2}_?\d+\.\d+)$", seqid)
            accession = m.group(1) if m else seqid.split("|")[-1]
            taxid_to_accessions[taxid].append(accession)
    return taxid_to_accessions


def ancestor_in(taxid, parent_of, target_set, max_depth=15):
    t = taxid
    for _ in range(max_depth):
        if t in target_set:
            return t
        if t not in parent_of or parent_of[t] == t:
            return None
        t = parent_of[t]
    return None


def main():
    parent_of, name_of = load_taxonomy()
    taxid_to_accessions = load_seqid2taxid()

    target_species = {}  # taxid -> (name, total_reads)
    with open(SUMMARY_CSV) as fh:
        for row in csv.DictReader(fh):
            if row["is_common_respiratory_virus"] == "True":
                target_species[row["taxid"]] = row["species"]

    species_totals = Counter()
    with open(SUMMARY_CSV) as fh:
        for row in csv.DictReader(fh):
            if row["taxid"] in target_species:
                species_totals[row["taxid"]] += int(row["n_reads"])
    target_species = {t: n for t, n in target_species.items() if species_totals[t] >= MIN_TOTAL_READS}

    leaf_counts_per_species = defaultdict(Counter)
    for out_path in glob.glob(os.path.join(REPORTS_DIR, "*.kraken2_output.txt")):
        with open(out_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if fields[0] != "C":
                    continue
                leaf_taxid = fields[2]
                species_taxid = ancestor_in(leaf_taxid, parent_of, target_species)
                if species_taxid:
                    leaf_counts_per_species[species_taxid][leaf_taxid] += 1

    rows = []
    for species_taxid, species_name in target_species.items():
        leaf_counts = leaf_counts_per_species.get(species_taxid)
        if not leaf_counts:
            rows.append({"species": species_name, "species_taxid": species_taxid,
                         "total_reads": species_totals[species_taxid],
                         "chosen_leaf_taxid": "", "chosen_accession": "",
                         "note": "no leaf taxid resolved"})
            continue
        best_leaf_taxid, n_reads_for_leaf = leaf_counts.most_common(1)[0]
        # NCBI taxonomy assigns one taxid per strain, not per segment, so a
        # segmented genome (e.g. influenza) has multiple accessions under the
        # same leaf taxid. Emit one row per accession/segment so step 5 fetches
        # the whole genome, not just whichever segment happened to be listed first.
        accessions = taxid_to_accessions.get(best_leaf_taxid, [])
        if not accessions:
            rows.append({"species": species_name, "species_taxid": species_taxid,
                         "total_reads": species_totals[species_taxid],
                         "chosen_leaf_taxid": best_leaf_taxid,
                         "chosen_leaf_name": name_of.get(best_leaf_taxid, ""),
                         "chosen_accession": "", "segment_index": "", "n_segments": 0,
                         "note": "no accession found in seqid2taxid.map"})
            continue
        for i, accession in enumerate(sorted(accessions), start=1):
            rows.append({
                "species": species_name,
                "species_taxid": species_taxid,
                "total_reads": species_totals[species_taxid],
                "chosen_leaf_taxid": best_leaf_taxid,
                "chosen_leaf_name": name_of.get(best_leaf_taxid, ""),
                "chosen_accession": accession,
                "segment_index": i,
                "n_segments": len(accessions),
                "note": "",
            })

    rows.sort(key=lambda r: (-r["total_reads"], r.get("segment_index") or 0))
    fieldnames = ["species", "species_taxid", "total_reads", "chosen_leaf_taxid",
                  "chosen_leaf_name", "chosen_accession", "segment_index", "n_segments", "note"]
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} reference picks (across {len(target_species)} species) -> {OUT_CSV}")
    for r in rows:
        seg = f" [segment {r['segment_index']}/{r['n_segments']}]" if r.get("n_segments", 0) > 1 else ""
        print(f"  {r['total_reads']:6d} reads  {r['species']:55s} -> {r['chosen_accession']}{seg} ({r['chosen_leaf_name']})")


if __name__ == "__main__":
    main()
