#!/usr/bin/env python
# ---------------------------------------------------------------------------
# rapids-singlecell GPU workflow for KazRNA-Pipe.
#
# Mirrors scripts/R/seurat_workflow.R stage-for-stage so that clustering
# agreement (ARI/NMI/ASW) reported in Figure 3C is a like-for-like
# comparison.
#
# Stages: load -> QC -> normalize -> HVG -> PCA -> Harmony -> neighbors ->
#         Leiden (multi-resolution) -> UMAP -> markers -> annotation.
#
# Per-stage wall-clock is recorded for Figure 3A.
#
# Container: ghcr.io/scverse/rapids-singlecell:0.10.10-cuda12.2
# ---------------------------------------------------------------------------
from __future__ import annotations

import argparse
import json
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

try:
    import rapids_singlecell as rsc
    GPU = True
except (ImportError, RuntimeError):
    rsc = None                                         # type: ignore
    GPU = False
    print("[warn] rapids-singlecell unavailable; falling back to Scanpy.",
          file=sys.stderr)

import scanpy as sc                                    # type: ignore
import anndata as ad                                   # type: ignore


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="KazRNA-Pipe rapids-singlecell GPU workflow")
    # Two mutually exclusive input modes. --input_h5 is what the pipeline uses;
    # --input_dir is kept so the script remains usable on 10x MTX directories.
    p.add_argument("--input_h5", default=None,
                   help="Comma-separated list of 10x .h5 files")
    p.add_argument("--sample_ids", default=None,
                   help="Comma-separated sample IDs, same order as --input_h5")
    p.add_argument("--input_dir", default=None,
                   help="Directory containing one subdir per sample with 10x matrices")
    p.add_argument("--markers", default=None,
                   help="YAML of cell-type marker genes for cluster annotation")
    p.add_argument("--outdir", default="results/sc/gpu")
    p.add_argument("--normalization", default="lognorm",
                   choices=["lognorm", "pearson"],
                   help="lognorm matches the Seurat CPU path; use it for any "
                        "CPU-vs-GPU comparison")
    p.add_argument("--mt_threshold", type=float, default=20.0)
    p.add_argument("--min_features", type=int, default=200)
    p.add_argument("--hvg_n", type=int, default=3000)
    p.add_argument("--n_pcs", type=int, default=50)
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)

    opt = p.parse_args()
    if not opt.input_h5 and not opt.input_dir:
        p.error("one of --input_h5 or --input_dir is required")
    if opt.input_h5 and not opt.sample_ids:
        p.error("--input_h5 requires --sample_ids")
    return opt


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


def load_samples(input_dir: Optional[str] = None,
                 input_h5: Optional[str] = None,
                 sample_ids: Optional[str] = None) -> ad.AnnData:
    """Load every sample, either from a list of .h5 files or from 10x MTX dirs."""
    if input_h5:
        files = [Path(p.strip()) for p in input_h5.split(",") if p.strip()]
        ids = [s.strip() for s in sample_ids.split(",") if s.strip()]
        if len(files) != len(ids):
            raise SystemExit(
                f"--input_h5 has {len(files)} entries but --sample_ids has {len(ids)}")
        missing = [str(f) for f in files if not f.exists()]
        if missing:
            raise SystemExit(f"Input file(s) not found: {', '.join(missing)}")
        adatas = [sc.read_10x_h5(f) for f in files]
    else:
        dirs = sorted([p for p in Path(input_dir).iterdir() if p.is_dir()])
        if not dirs:
            raise SystemExit(f"No sample subdirectories under {input_dir}")
        ids = [d.name for d in dirs]
        adatas = [sc.read_10x_mtx(d, var_names="gene_symbols", cache=False) for d in dirs]

    for a, s in zip(adatas, ids):
        a.var_names_make_unique()
        a.obs["sample_id"] = s

    return ad.concat(adatas, join="outer", label="sample_id",
                     keys=ids, index_unique="-")


def qc(adata: ad.AnnData, mt_threshold: float, min_features: int) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    keep = (adata.obs["n_genes_by_counts"] >= min_features) & \
           (adata.obs["pct_counts_mt"] < mt_threshold)
    return adata[keep].copy()


def normalize(adata: ad.AnnData, method: str = "lognorm") -> ad.AnnData:
    """
    lognorm: library-size scaling to 1e4 + log1p. Identical to the Seurat CPU
             path under --normalization lognorm, so the two are comparable.
    pearson: analytic Pearson residuals, the closest Scanpy analogue of
             SCTransform. Do NOT use this for the CPU-vs-GPU comparison unless
             the CPU path is also run with --normalization sctransform.
    """
    adata.layers["counts"] = adata.X.copy()
    if method == "pearson":
        sc.experimental.pp.normalize_pearson_residuals(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


def hvg(adata: ad.AnnData, n_top: int) -> ad.AnnData:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3",
                                layer="counts")
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


def neighbors(adata: ad.AnnData, n_pcs: int, n_neighbors: int) -> ad.AnnData:
    if GPU:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                         use_rep="X_pca_harmony")
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        use_rep="X_pca_harmony")
    return adata


def leiden_multi(adata: ad.AnnData, resolutions: List[float], seed: int) -> ad.AnnData:
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


def markers(adata: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")
    rgg = adata.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)

    rows = []
    for g in groups:
        for i, gene in enumerate(rgg["names"][g]):
            rows.append({
                "cluster": g,
                "gene": gene,
                "logfoldchange": float(rgg["logfoldchanges"][g][i]),
                "p_val": float(rgg["pvals"][g][i]),
                "p_val_adj": float(rgg["pvals_adj"][g][i]),
            })
    return pd.DataFrame(rows)


def annotate(adata: ad.AnnData, marker_yaml: Optional[str],
             cluster_key: str) -> pd.DataFrame:
    """
    Score every cluster against each cell-type marker set and assign the
    highest-scoring type. Returns a cluster -> cell_type table. Falls back to
    'unassigned' when no marker file is supplied, so the output file always
    exists and the process contract is stable.
    """
    clusters = sorted(adata.obs[cluster_key].astype(str).unique(),
                      key=lambda x: (len(x), x))

    if not marker_yaml or not Path(marker_yaml).exists():
        print(f"[warn] no marker file at {marker_yaml}; clusters left unassigned",
              file=sys.stderr)
        return pd.DataFrame({"cluster": clusters, "cell_type": "unassigned",
                             "score": np.nan})

    import yaml                                          # type: ignore
    marker_sets = yaml.safe_load(Path(marker_yaml).read_text())

    score_cols = {}
    for cell_type, genes in marker_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        col = f"score_{cell_type}"
        sc.tl.score_genes(adata, gene_list=present, score_name=col)
        score_cols[cell_type] = col

    if not score_cols:
        print("[warn] no marker genes found in the data; clusters left unassigned",
              file=sys.stderr)
        return pd.DataFrame({"cluster": clusters, "cell_type": "unassigned",
                             "score": np.nan})

    mean_scores = adata.obs.groupby(cluster_key, observed=True)[
        list(score_cols.values())].mean()
    inv = {v: k for k, v in score_cols.items()}

    return pd.DataFrame({
        "cluster": mean_scores.index.astype(str),
        "cell_type": [inv[c] for c in mean_scores.idxmax(axis=1)],
        "score": mean_scores.max(axis=1).values,
    }).reset_index(drop=True)


def main() -> None:
    opt = parse_args()
    outdir = Path(opt.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    resolutions = [float(r) for r in opt.resolutions.split(",")]
    primary_res = 1.0 if 1.0 in resolutions else resolutions[len(resolutions) // 2]
    cluster_key = f"leiden_res_{primary_res:.2f}"

    sc.settings.n_jobs = opt.threads
    np.random.seed(opt.seed)

    timer = StageTimer()

    adata = timer("load", load_samples, opt.input_dir, opt.input_h5, opt.sample_ids)
    adata = timer("qc", qc, adata, opt.mt_threshold, opt.min_features)
    adata = timer("normalize", normalize, adata, opt.normalization)
    adata = timer("hvg", hvg, adata, opt.hvg_n)
    adata = timer("pca", run_pca, adata, opt.n_pcs, opt.seed)
    adata = timer("harmony", run_harmony, adata)
    adata = timer("neighbors", neighbors, adata, opt.n_pcs, opt.n_neighbors)
    adata = timer("leiden", leiden_multi, adata, resolutions, opt.seed)
    adata = timer("umap", run_umap, adata, opt.seed)
    marker_df = timer("markers", markers, adata, cluster_key)
    celltype_df = timer("annotate", annotate, adata, opt.markers, cluster_key)

    # -------- Save outputs --------
    adata.write_h5ad(outdir / "adata.h5ad", compression="gzip")
    marker_df.to_csv(outdir / "markers.tsv", sep="\t", index=False)
    celltype_df.to_csv(outdir / "celltype_assignments.tsv", sep="\t", index=False)

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
        "primary_cluster_key": cluster_key,
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
