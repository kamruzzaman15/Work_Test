#!/usr/bin/env python3
"""Part 2, Step 7: aggregate per-pool samtools-coverage tables into one
pool x virus breadth/depth matrix, and plot per-position depth tracks for
the viruses with the most supporting reads.

Answers: "what kind of coverage across the genome do we see for common
respiratory viruses?"
"""
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE = os.path.join(os.path.dirname(__file__), "..")
COV_DIR = os.path.join(BASE, "results", "coverage")
PLOTS_DIR = os.path.join(BASE, "results", "plots")
OUT_BREADTH_CSV = os.path.join(BASE, "results", "part2_coverage_breadth_depth_matrix.csv")
OUT_SEGMENT_CSV = os.path.join(BASE, "results", "part2_coverage_per_segment.csv")

os.makedirs(PLOTS_DIR, exist_ok=True)


def parse_contig_name(rname):
    """'NC_038311.1|Rhinovirus_A' or '...|Gammainfluenzavirus_influenzae|seg2' -> (accession, virus, segment)"""
    parts = rname.split("|")
    accession = parts[0]
    virus = parts[1].replace("_", " ")
    segment = parts[2] if len(parts) > 2 else None
    return accession, virus, segment


def load_all_coverage():
    rows = []
    for path in sorted(glob.glob(os.path.join(COV_DIR, "*.coverage_summary.tsv"))):
        pool = os.path.basename(path).replace(".coverage_summary.tsv", "")
        df = pd.read_csv(path, sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        df["pool"] = pool
        df[["accession", "virus", "segment"]] = df["rname"].apply(lambda r: pd.Series(parse_contig_name(r)))
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def main():
    cov = load_all_coverage()
    cov.to_csv(OUT_SEGMENT_CSV, index=False)
    print(f"wrote per-pool, per-segment coverage table -> {OUT_SEGMENT_CSV}")

    # Genome-level (not per-segment) breadth/depth: for multi-segment viruses,
    # weight by segment length so one virus = one row per pool.
    def genome_level(group):
        seg_len = group["endpos"] - group["startpos"] + 1
        total_len = seg_len.sum()
        covered = (group["coverage"] / 100.0 * seg_len).sum()
        depth_weighted = (group["meandepth"] * seg_len).sum() / total_len
        return pd.Series({
            "breadth_pct": round(covered / total_len * 100, 2),
            "mean_depth": round(depth_weighted, 2),
            "n_segments": len(group),
            "total_reads_mapped": int(group["numreads"].sum()),
            "genome_length": int(total_len),
        })

    genome_cov = cov.groupby(["pool", "virus"]).apply(genome_level, include_groups=False).reset_index()
    matrix = genome_cov.pivot_table(index="virus", columns="pool", values="breadth_pct", fill_value=0.0)
    matrix.to_csv(OUT_BREADTH_CSV)
    print(f"wrote pool x virus breadth%% matrix -> {OUT_BREADTH_CSV}")

    depth_matrix = genome_cov.pivot_table(index="virus", columns="pool", values="mean_depth", fill_value=0.0)
    depth_out = OUT_BREADTH_CSV.replace("breadth_depth_matrix", "depth_matrix")
    depth_matrix.to_csv(depth_out)
    print(f"wrote pool x virus mean-depth matrix -> {depth_out}")

    # Plot depth-vs-position for the top 6 viruses by total mapped reads,
    # using whichever single pool gave that virus the most reads (most
    # informative track for that virus).
    top_viruses = (
        genome_cov.groupby("virus")["total_reads_mapped"].sum().sort_values(ascending=False).head(6).index.tolist()
    )
    print("plotting depth tracks for:", top_viruses)

    for virus in top_viruses:
        sub = genome_cov[genome_cov["virus"] == virus]
        best_pool = sub.loc[sub["total_reads_mapped"].idxmax(), "pool"]
        plot_depth_track(virus, best_pool)


def plot_depth_track(virus, pool):
    depth_path = os.path.join(COV_DIR, f"{pool}.depth.tsv")
    depth = pd.read_csv(depth_path, sep="\t", header=None, names=["rname", "pos", "depth"])
    virus_key = virus.replace(" ", "_")
    depth = depth[depth["rname"].str.contains(re.escape(virus_key), regex=True)]
    if depth.empty:
        return

    segments = sorted(depth["rname"].unique())
    fig, axes = plt.subplots(len(segments), 1, figsize=(8, 2 * len(segments)), squeeze=False)
    for ax, seg_name in zip(axes[:, 0], segments):
        seg_depth = depth[depth["rname"] == seg_name]
        ax.fill_between(seg_depth["pos"], seg_depth["depth"], step="mid")
        ax.set_ylabel("depth")
        label = seg_name.split("|")[-1] if "|seg" in seg_name else "genome"
        ax.set_title(f"{virus} ({label}) - pool {pool}", fontsize=9)
    axes[-1, 0].set_xlabel("genome position (bp)")
    fig.tight_layout()
    out_path = os.path.join(PLOTS_DIR, f"part2_coverage_{virus_key}.png")
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  saved {out_path}")


if __name__ == "__main__":
    main()
