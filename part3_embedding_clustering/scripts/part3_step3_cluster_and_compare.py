#!/usr/bin/env python3
"""Part 3, Step 3: turn the sourmash pairwise similarity matrix into a 2D
embedding (classical MDS) and hard clusters, then compare those clusters
against both things the task asks about:
  - taxonomy profile: Kraken2 species labels from Part 2 (via part3 step 1)
  - ANI: sourmash similarity is itself an ANI estimate, so we check whether
    same-species reads have higher pairwise similarity than different-species
    reads (i.e. is the embedding/cluster structure ANI-coherent?)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score

BASE = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE, "..", "results")
PLOTS_DIR = os.path.join(RESULTS_DIR, "plots")
SIM_CSV = os.path.join(RESULTS_DIR, "part3_similarity_matrix.csv")
LABELS_CSV = os.path.join(BASE, "..", "data", "selected_reads", "part3_read_labels.csv")

os.makedirs(PLOTS_DIR, exist_ok=True)


def main():
    sim = pd.read_csv(SIM_CSV, header=0)
    read_ids = list(sim.columns)
    sim_matrix = sim.values
    sim_matrix = (sim_matrix + sim_matrix.T) / 2
    np.fill_diagonal(sim_matrix, 1.0)

    labels_df = pd.read_csv(LABELS_CSV).set_index("read_id").loc[read_ids]
    species = labels_df["species_name"].values
    n_species = len(set(species))
    print(f"{len(read_ids)} reads, {n_species} species")

    dist_matrix = 1.0 - sim_matrix
    np.fill_diagonal(dist_matrix, 0.0)
    dist_matrix[dist_matrix < 0] = 0.0  # guard against float noise

    # --- 2D embedding for plotting ---
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0,
              normalized_stress="auto", n_init=4)
    coords = mds.fit_transform(dist_matrix)

    best_eps, best_ari, best_labels = None, -1, None
    scan_log = []
    for eps in [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]:
        labels_try = DBSCAN(eps=eps, min_samples=3, metric="precomputed").fit_predict(dist_matrix)
        ari_try = adjusted_rand_score(species, labels_try)
        n_clusters_try = len(set(labels_try)) - (1 if -1 in labels_try else 0)
        n_noise_try = int((labels_try == -1).sum())
        line = (f"  DBSCAN eps={eps:.2f}: {n_clusters_try} clusters, "
                f"{n_noise_try} noise points, ARI={ari_try:.3f}")
        print(line)
        scan_log.append(line)
        if ari_try > best_ari:
            best_eps, best_ari, best_labels = eps, ari_try, labels_try

    cluster_labels = best_labels
    ari = best_ari
    n_clusters_final = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise_final = int((cluster_labels == -1).sum())
    print(f"selected eps={best_eps}")
    print(f"Adjusted Rand Index (clusters vs. Kraken2 taxonomy): {ari:.3f}")

    species_arr = np.array(species)
    same = species_arr[:, None] == species_arr[None, :]
    iu = np.triu_indices_from(sim_matrix, k=1)
    same_pairs = same[iu]
    sim_pairs = sim_matrix[iu]
    within_mean = sim_pairs[same_pairs].mean()
    between_mean = sim_pairs[~same_pairs].mean()
    print(f"mean within-species similarity:  {within_mean:.4f}")
    print(f"mean between-species similarity: {between_mean:.4f}")

    per_species_stats = []
    for sp in sorted(set(species)):
        mask = species_arr == sp
        idx = np.where(mask)[0]
        if len(idx) < 2:
            continue
        sub = sim_matrix[np.ix_(idx, idx)]
        within = sub[np.triu_indices_from(sub, k=1)].mean()
        other_idx = np.where(~mask)[0]
        between = sim_matrix[np.ix_(idx, other_idx)].mean()
        per_species_stats.append({"species": sp, "n_reads": len(idx),
                                   "mean_within_species_similarity": round(within, 4),
                                   "mean_between_species_similarity": round(between, 4)})
    stats_df = pd.DataFrame(per_species_stats).sort_values("n_reads", ascending=False)
    stats_path = os.path.join(RESULTS_DIR, "part3_within_between_similarity.csv")
    stats_df.to_csv(stats_path, index=False)
    print(f"wrote per-species within/between similarity -> {stats_path}")

    summary_path = os.path.join(RESULTS_DIR, "part3_cluster_summary.csv")
    pd.DataFrame({
        "read_id": read_ids,
        "species_name": species,
        "cluster": cluster_labels,
        "mds_x": coords[:, 0],
        "mds_y": coords[:, 1],
    }).to_csv(summary_path, index=False)
    print(f"wrote per-read cluster assignments + MDS coords -> {summary_path}")

    with open(os.path.join(RESULTS_DIR, "part3_ari_score.txt"), "w") as fh:
        fh.write("DBSCAN eps scan (selecting eps against known Kraken2 labels):\n")
        for line in scan_log:
            fh.write(line.strip() + "\n")
        fh.write(f"\nselected eps={best_eps}\n")
        fh.write(f"final model: {n_clusters_final} clusters, {n_noise_final} noise points (of {len(read_ids)} reads)\n")
        fh.write(f"Adjusted Rand Index (unsupervised clusters vs. Kraken2 taxonomy labels): {ari:.4f}\n")
        fh.write(f"mean within-species sourmash similarity:  {within_mean:.4f}\n")
        fh.write(f"mean between-species sourmash similarity: {between_mean:.4f}\n")

    plot_mds(coords, species, "part3_mds_by_taxonomy.png", "MDS of reads, colored by Kraken2 species (Part 2)")
    plot_mds(coords, cluster_labels.astype(str), "part3_mds_by_cluster.png",
              f"MDS of reads, colored by unsupervised cluster (ARI={ari:.2f})")


def plot_mds(coords, color_labels, fname, title):
    fig, ax = plt.subplots(figsize=(9, 7))
    unique_labels = sorted(set(color_labels), key=str)
    cmap = plt.get_cmap("tab20", len(unique_labels))
    for i, lab in enumerate(unique_labels):
        mask = np.array(color_labels) == lab
        ax.scatter(coords[mask, 0], coords[mask, 1], s=14, alpha=0.7,
                   color=cmap(i), label=str(lab))
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("MDS dim 1")
    ax.set_ylabel("MDS dim 2")
    ax.legend(fontsize=6, markerscale=1.5, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    out_path = os.path.join(PLOTS_DIR, fname)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"  saved {out_path}")


if __name__ == "__main__":
    main()
