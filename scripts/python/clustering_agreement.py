#!/usr/bin/env python
# ---------------------------------------------------------------------------
# Clustering agreement between the CPU (Seurat) and GPU (rapids-singlecell)
# paths of KazRNA-Pipe, at every Leiden resolution.
#
# Outputs:
#   ari_nmi_asw_by_resolution.tsv     Main numerical table
#   ari_nmi_asw_by_resolution.pdf     Figure 3C of the manuscript
#   provenance.json
#
# The ARI threshold of 0.99 between two independent GPU runs (mentioned in
# manuscript Section 4.4) is also checked here when --self_check is passed.
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import json
import time
import socket
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cpu_labels", required=True,
                   help="Long-format TSV with columns: cell_id, resolution, cluster")
    p.add_argument("--gpu_labels", required=True,
                   help="Same schema as --cpu_labels")
    p.add_argument("--embedding",  required=False,
                   help="(Optional) joint embedding TSV for ASW; cell_id + numeric columns")
    p.add_argument("--outdir",     default="results/sc/agreement")
    p.add_argument("--self_check", action="store_true",
                   help="If both label files are GPU-vs-GPU, assert ARI > 0.99 "
                        "for every resolution (manuscript Section 4.4)")
    return p.parse_args()


def read_long(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"cell_id", "resolution", "cluster"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    return df


def main() -> None:
    opt = parse_args()
    outdir = Path(opt.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cpu = read_long(opt.cpu_labels)
    gpu = read_long(opt.gpu_labels)

    common_cells = sorted(set(cpu["cell_id"]) & set(gpu["cell_id"]))
    if len(common_cells) < 100:
        raise SystemExit(f"Only {len(common_cells)} cells in common; aborting.")

    cpu = cpu[cpu["cell_id"].isin(common_cells)].set_index(["resolution", "cell_id"])
    gpu = gpu[gpu["cell_id"].isin(common_cells)].set_index(["resolution", "cell_id"])

    resolutions = sorted(set(cpu.index.get_level_values("resolution")) &
                          set(gpu.index.get_level_values("resolution")))

    embedding = None
    if opt.embedding:
        emb_df = pd.read_csv(opt.embedding, sep="\t").set_index("cell_id")
        embedding = emb_df.loc[common_cells].to_numpy()

    rows = []
    for r in resolutions:
        c = cpu.xs(r, level="resolution").reindex(common_cells)["cluster"].astype(str).values
        g = gpu.xs(r, level="resolution").reindex(common_cells)["cluster"].astype(str).values

        ari = adjusted_rand_score(c, g)
        nmi = normalized_mutual_info_score(c, g)
        asw_cpu = silhouette_score(embedding, c) if embedding is not None else np.nan
        asw_gpu = silhouette_score(embedding, g) if embedding is not None else np.nan

        rows.append({
            "resolution": r,
            "ARI": ari,
            "NMI": nmi,
            "ASW_cpu": asw_cpu,
            "ASW_gpu": asw_gpu,
            "n_clusters_cpu": len(set(c)),
            "n_clusters_gpu": len(set(g)),
        })

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "ari_nmi_asw_by_resolution.tsv", sep="\t", index=False)

    # ----- Self-check assertion --------
    if opt.self_check:
        bad = df[df["ARI"] < 0.99]
        if not bad.empty:
            raise SystemExit(
                f"Self-check failed: ARI dropped below 0.99 at resolutions {bad['resolution'].tolist()}"
            )
        print("[self-check] All resolutions have ARI > 0.99. PASS.")

    # ----- Plot (manuscript Figure 3C) --------
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["resolution"], df["ARI"], "-o", label="ARI")
    ax.plot(df["resolution"], df["NMI"], "-s", label="NMI")
    if embedding is not None:
        ax.plot(df["resolution"], df["ASW_gpu"], "-^", label="ASW (GPU)")
    ax.axhline(0.95, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("CPU vs GPU agreement")
    ax.set_ylim(0.5, 1.02)
    ax.legend(loc="lower left")
    ax.set_title("Clustering agreement between Seurat (CPU) and rapids-singlecell (GPU)")
    fig.tight_layout()
    fig.savefig(outdir / "ari_nmi_asw_by_resolution.pdf")
    plt.close(fig)

    # ----- Provenance --------
    provenance = {
        "script": "clustering_agreement.py",
        "datetime": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "hostname": socket.gethostname(),
        "arguments": vars(opt),
        "n_common_cells": len(common_cells),
        "resolutions": resolutions,
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2))

    print(f"Wrote agreement metrics to {outdir}")


if __name__ == "__main__":
    main()
