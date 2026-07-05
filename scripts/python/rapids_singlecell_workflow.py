#!/usr/bin/env python
# ---------------------------------------------------------------------------
# rapids-singlecell GPU workflow for KazRNA-Pipe.
#
# Mirrors scripts/R/seurat_workflow.R stage-for-stage so that clustering
# agreement (ARI/NMI/ASW) reported in manuscript Figure 3C is a like-for-like
# comparison.
#
# Stages: load -> QC -> normalize -> HVG -> PCA -> Harmony -> neighbors ->
#         Leiden (multi-resolution) -> UMAP -> markers.
#
# Per-stage wall-clock is recorded for manuscript Figure 3A.
#
# Container: ghcr.io/scverse/rapids-singlecell:0.10.10-cuda12.2
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Dict

import numpy as np
import pandas as pd

# rapids-singlecell pulls in cuML, cuGraph and Scanpy together. If we
# cannot import it (e.g. CI without GPU), fall back to Scanpy so the script
# is still usable for testing pipeline plumbing.
try:
    import rapids_singlecell as rsc
    GPU = True
except (ImportError, RuntimeError):
    import scanpy as sc                                # type: ignore
    rsc = None                                         # type: ignore
    GPU = False
    print("[warn] rapids-singlecell unavailable; falling back to Scanpy.",
          file=sys.stderr)

import scanpy as sc                                    # type: ignore
import anndata as ad                                   # type: ignore


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="KazRNA-Pipe rapids-singlecell GPU workflow")
    p.add_argument("--input_dir", required=True,
                   help="Directory containing one subdir per sample with 10x matrices")
    p.add_argument("--outdir", default="results/sc/gpu")
    p.add_argument("--mt_threshold", type=float, default=20.0)
    p.add_argument("--min_features", type=int, default=200)
    p.add_argument("--hvg_n", type=int, default=3000)
    p.add_argument("--n_pcs", type=int, default=50)
    p.add_argument("--resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


class StageTimer:
    """Records wall-clock for each pipeline stage, written to timing.tsv."""

    def __init__(self) -> None:
        self.records: Dict[str, float] = {}

    def __call__(self, label: str, fn: Callable, *args, **kwargs):
        t0 = time.perf_counter()
        result = fn(*args, **kwargs)
        dt = time.perf_counter() - t0
        self.records[label] = dt
        print(f"[stage] {label:<15} {dt:7.1f} s", flush=True)
        return result


def load_samples(input_dir: Path) -> ad.AnnData:
    """Load every sample under input_dir/<sample_id>/ as 10x matrices."""
    sample_dirs = sorted([p for p in input_dir.iterdir() if p.is_dir()])
    if not sample_dirs:
        raise SystemExit(f"No sample subdirectories under {input_dir}")

    adatas = []
    for d in sample_dirs:
        ad_i = sc.read_10x_mtx(d, var_names="gene_symbols", cache=False)
        ad_i.obs["sample_id"] = d.name
        ad_i.var_names_make_unique()
        adatas.append(ad_i)

    adata = ad.concat(adatas, join="outer", label="sample_id",
                      keys=[d.name for d in sample_dirs], index_unique="-")
    return adata


def qc(adata: ad.AnnData, mt_threshold: float, min_features: int) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    keep = (adata.obs["n_genes_by_counts"] >= min_features) & \
           (adata.obs["pct_counts_mt"] < mt_threshold)
    return adata[keep].copy()


def normalize(adata: ad.AnnData) -> ad.AnnData:
    """Median library-size scaling + log1p, matching the Seurat target."""
    target = float(np.median(adata.obs["total_counts"]))
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)
    return adata


def hvg(adata: ad.AnnData, n_top: int) -> ad.AnnData:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3_paper")
    return adata


def run_pca(adata: ad.AnnData, n_pcs: int, seed: int) -> ad.AnnData:
    if GPU:
        rsc.tl.pca(adata, n_comps=n_pcs, random_state=seed)
    else:
        sc.tl.pca(adata, n_comps=n_pcs, random_state=seed)
    return adata


def run_harmony(adata: ad.AnnData, batch_key: str = "sample_id") -> ad.AnnData:
    if GPU:
        rsc.pp.harmony_integrate(adata, key=batch_key)
    else:
        sc.external.pp.harmony_integrate(adata, key=batch_key)
    return adata


def neighbors(adata: ad.AnnData, n_pcs: int) -> ad.AnnData:
    if GPU:
        rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep="X_pca_harmony")
    else:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep="X_pca_harmony")
    return adata


def leiden_multi(adata: ad.AnnData, resolutions: list[float], seed: int) -> ad.AnnData:
    for r in resolutions:
        key = f"leiden_res_{r:.2f}"
        if GPU:
            rsc.tl.leiden(adata, resolution=r, key_added=key, random_state=seed)
        else:
            sc.tl.leiden(adata, resolution=r, key_added=key, random_state=seed,
                         flavor="igraph", directed=False, n_iterations=2)
    return adata


def run_umap(adata: ad.AnnData, seed: int) -> ad.AnnData:
    if GPU:
        rsc.tl.umap(adata, random_state=seed)
    else:
        sc.tl.umap(adata, random_state=seed)
    return adata


def markers(adata: ad.AnnData) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, "leiden_res_1.00", method="wilcoxon")
    rgg = adata.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)

    rows = []
    for g in groups:
        for i, gene in enumerate(rgg["names"][g]):
            rows.append({
                "cluster": g,
                "gene": gene,
                "logfoldchange": float(rgg["logfoldchanges"][g][i]),
                "p_val":   float(rgg["pvals"][g][i]),
                "p_val_adj": float(rgg["pvals_adj"][g][i]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    opt = parse_args()
    outdir = Path(opt.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    resolutions = [float(r) for r in opt.resolutions.split(",")]

    timer = StageTimer()

    adata = timer("load", load_samples, Path(opt.input_dir))
    adata = timer("qc", qc, adata, opt.mt_threshold, opt.min_features)
    adata = timer("normalize", normalize, adata)
    adata = timer("hvg", hvg, adata, opt.hvg_n)
    adata = timer("pca", run_pca, adata, opt.n_pcs, opt.seed)
    adata = timer("harmony", run_harmony, adata)
    adata = timer("neighbors", neighbors, adata, opt.n_pcs)
    adata = timer("leiden", leiden_multi, adata, resolutions, opt.seed)
    adata = timer("umap", run_umap, adata, opt.seed)
    marker_df = timer("markers", markers, adata)

    # -------- Save outputs --------
    adata.write_h5ad(outdir / "adata.h5ad", compression="gzip")
    marker_df.to_csv(outdir / "markers.tsv", sep="\t", index=False)

    # Long-format labels for downstream agreement analysis
    labels_rows = []
    for r in resolutions:
        col = f"leiden_res_{r:.2f}"
        if col not in adata.obs.columns:
            continue
        labels_rows.append(pd.DataFrame({
            "cell_id": adata.obs_names,
            "resolution": r,
            "cluster": adata.obs[col].astype(str).values,
        }))
    pd.concat(labels_rows).to_csv(outdir / "labels.tsv", sep="\t", index=False)

    # UMAP coords
    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"],
                           index=adata.obs_names).reset_index().rename(
                           columns={"index": "cell_id"})
    umap_df.to_csv(outdir / "umap_coords.tsv", sep="\t", index=False)

    # -------- Timing for Figure 3A --------
    timing_df = pd.DataFrame({
        "stage": list(timer.records.keys()),
        "wall_seconds": list(timer.records.values()),
        "path": "gpu" if GPU else "cpu_fallback",
    })
    timing_df.to_csv(outdir / "timing.tsv", sep="\t", index=False)

    # -------- Provenance --------
    def _try(cmd):
        try:
            return subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL,
                                           text=True).strip()
        except subprocess.CalledProcessError:
            return None

    provenance = {
        "script": "rapids_singlecell_workflow.py",
        "git_describe": _try("git describe --always --dirty"),
        "datetime": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "hostname": socket.gethostname(),
        "gpu_path_used": GPU,
        "arguments": vars(opt),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "library_versions": {
            "python": sys.version.split()[0],
            "rapids_singlecell": getattr(rsc, "__version__", None) if GPU else None,
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "nvidia_smi": _try("nvidia-smi --query-gpu=name,driver_version,memory.total "
                           "--format=csv,noheader") if GPU else None,
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2, default=str))

    print(f"rapids-singlecell workflow complete. Output in: {outdir}")


if __name__ == "__main__":
    main()
