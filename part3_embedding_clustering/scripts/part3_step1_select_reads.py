#!/usr/bin/env python3
"""Part 3, Step 1: select a balanced set of reads per common respiratory virus
for embedding/clustering, reusing Part 2's existing Kraken2 taxonomy calls as
ground truth (no new classification work).

Output:
  - data/selected_reads/part3_selected_reads.fasta   (one record per chosen read)
  - data/selected_reads/part3_read_labels.csv         (read_id, pool, species, taxid)
"""
import glob
import gzip
import os
import random
import re
from collections import defaultdict

random.seed(0)

PART2 = os.path.join(os.path.dirname(__file__), "..", "..", "part2_taxonomic_classification")
KRAKEN_REPORTS_DIR = os.path.join(PART2, "results", "kraken2_reports")
RAW_POOLS_DIR = os.path.join(PART2, "data", "raw_pools")
SUMMARY_CSV = os.path.join(PART2, "results", "part2_virus_summary.csv")

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "selected_reads")
OUT_FASTA = os.path.join(OUT_DIR, "part3_selected_reads.fasta")
OUT_LABELS = os.path.join(OUT_DIR, "part3_read_labels.csv")

MIN_TOTAL_READS = 20  
MAX_READS_PER_SPECIES = 60  

os.makedirs(OUT_DIR, exist_ok=True)


def load_common_respiratory_species():
    species_taxids = set()
    import csv
    with open(SUMMARY_CSV) as fh:
        totals = defaultdict(int)
        for row in csv.DictReader(fh):
            if row["is_common_respiratory_virus"] == "True":
                totals[row["taxid"]] += int(row["n_reads"])
        for taxid, n in totals.items():
            if n >= MIN_TOTAL_READS:
                species_taxids.add(taxid)
    return species_taxids


def collect_reads_per_species(target_taxids):
    """read_id -> (pool, species_taxid) for every classified read in any pool
    whose leaf taxid rolls up to one of our target species. We reuse the same
    taxonomy-climbing trick as Part 2 step 4."""
    parent_of, name_of = load_taxonomy()
    reads_by_species = defaultdict(list)  

    for out_path in sorted(glob.glob(os.path.join(KRAKEN_REPORTS_DIR, "*.kraken2_output.txt"))):
        pool = os.path.basename(out_path).replace(".kraken2_output.txt", "")
        with open(out_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if fields[0] != "C":
                    continue
                read_id, leaf_taxid = fields[1], fields[2]
                species_taxid = ancestor_in(leaf_taxid, parent_of, target_taxids)
                if species_taxid:
                    reads_by_species[species_taxid].append((pool, read_id))
    return reads_by_species, name_of


def load_taxonomy():
    db_dir = os.path.join(PART2, "data", "kraken2_db")
    parent_of, name_of = {}, {}
    with open(os.path.join(db_dir, "ktaxonomy.tsv")) as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split("|")]
            parent_of[parts[0]] = parts[1]
            name_of[parts[0]] = parts[4]
    return parent_of, name_of


def ancestor_in(taxid, parent_of, target_set, max_depth=15):
    t = taxid
    for _ in range(max_depth):
        if t in target_set:
            return t
        if t not in parent_of or parent_of[t] == t:
            return None
        t = parent_of[t]
    return None


def load_pool_sequences(pool, wanted_ids):
    """Pull specific reads (by ID) out of a pool's gzipped FASTA."""
    path = os.path.join(RAW_POOLS_DIR, f"{pool}.respiratory.fasta.gz")
    found = {}
    wanted = set(wanted_ids)
    if not wanted:
        return found
    with gzip.open(path, "rt") as fh:
        read_id = None
        seq_chunks = []
        for line in fh:
            if line.startswith(">"):
                if read_id in wanted:
                    found[read_id] = "".join(seq_chunks)
                read_id = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if read_id in wanted:
            found[read_id] = "".join(seq_chunks)
    return found


def main():
    target_taxids = load_common_respiratory_species()
    print(f"{len(target_taxids)} common respiratory virus species pass the >= {MIN_TOTAL_READS} read threshold")

    reads_by_species, name_of = collect_reads_per_species(target_taxids)


    chosen_by_pool = defaultdict(list)  # pool -> [(read_id, species_taxid), ...]
    label_rows = []
    for taxid, candidates in reads_by_species.items():
        random.shuffle(candidates)
        chosen = candidates[:MAX_READS_PER_SPECIES]
        for pool, read_id in chosen:
            chosen_by_pool[pool].append((read_id, taxid))

    with open(OUT_FASTA, "w") as fasta_out, open(OUT_LABELS, "w") as labels_out:
        labels_out.write("read_id,pool,species_taxid,species_name\n")
        for pool, read_list in chosen_by_pool.items():
            wanted_ids = [r for r, _ in read_list]
            seqs = load_pool_sequences(pool, wanted_ids)
            for read_id, taxid in read_list:
                seq = seqs.get(read_id)
                if not seq:
                    continue
                species_name = name_of.get(taxid, "unknown")
                safe_id = re.sub(r"\s+", "_", f"{pool}__{read_id}")
                fasta_out.write(f">{safe_id}\n{seq}\n")
                labels_out.write(f"{safe_id},{pool},{taxid},{species_name}\n")
                label_rows.append((safe_id, species_name))

    by_species_count = defaultdict(int)
    for _, sp in label_rows:
        by_species_count[sp] += 1
    print(f"wrote {len(label_rows)} reads -> {OUT_FASTA}")
    for sp, n in sorted(by_species_count.items(), key=lambda x: -x[1]):
        print(f"  {n:4d}  {sp}")


if __name__ == "__main__":
    main()
