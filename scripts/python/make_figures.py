#!/usr/bin/env python3
"""scripts/python/make_figures.py

The script consults the *results* directory produced by Nextflow and writes each
figure into *out*. Designed to be re-runnable: deterministic seeding, no
network access, no hidden state.

Mapping of figures → upstream artefacts:

  Fig 1   architecture diagram (static SVG; copied from docs/architecture.svg)
  Fig 2A  cross-aligner Spearman heatmap   ← concordance/fig2a_*
  Fig 2B  cross-method DEG Venn            ← concordance/fig2b_*
  Fig 3A  UMAP, CPU vs GPU                 ← scrnaseq/cpu, scrnaseq/gpu
  Fig 3B  cluster size distribution        ← scrnaseq/*/clusters.tsv
  Fig 3C  ARI/NMI curve across resolutions ← scrnaseq/agreement/agreement_metrics.tsv
  Fig 4A  strong scaling curve             ← benchmark/fig4_scaling.tsv
  Fig 4B  weak scaling curve               ← benchmark/fig4_scaling.tsv
  Fig 4C  GPU memory trace                 ← scrnaseq/gpu/*/gpu_memory_trace.tsv
  Fig 4D  CPU vs GPU energy                ← benchmark/fig4d_energy.tsv
  Fig 5A  deconvolution proportions stacked bar ← deconvolution/*/proportions.tsv
  Fig 5B  cross-deconvolver correlation    ← deconvolution/agreement/fig5b_*
  Fig 5C  signature-by-cell heatmap        ← scrnaseq/*/celltype_assignments.tsv
  Fig 5D  epithelial volcano               ← celltype_de/fig5d_*
  Fig 6   KEGG dot-plot                    ← enrichment/*/kegg_enrichment.tsv
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG = logging.getLogger("make_figures")


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="[%(asctime)s] %(levelname)s %(message)s")


# ---------------------------------------------------------------------------
def fig1(results: Path, out: Path) -> None:
    """Architecture diagram is a hand-crafted SVG; copy verbatim."""
    src = Path("docs/architecture.svg")
    if src.exists():
        shutil.copy(src, out / "fig1_architecture.svg")
        LOG.info("fig1: copied %s", src)
    else:
        LOG.warning("fig1: docs/architecture.svg missing - skipping")


def fig2a(results: Path, out: Path) -> None:
    f = results / "concordance" / "fig2a_spearman_matrix.tsv"
    if not f.exists():
        LOG.warning("fig2a: %s missing", f); return
    m = pd.read_csv(f, sep="\t", index_col=0)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(m.values, vmin=0.85, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(m.columns))); ax.set_xticklabels(m.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(m.index)));   ax.set_yticklabels(m.index)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m.values[i,j]:.2f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if m.values[i,j] < 0.95 else "black")
    fig.colorbar(im, ax=ax, label="Spearman ρ", shrink=0.7)
    ax.set_title("Cross-aligner gene-count correlation")
    fig.tight_layout(); fig.savefig(out / "fig2a_spearman.pdf"); plt.close(fig)
    LOG.info("fig2a written")


def fig2b(results: Path, out: Path) -> None:
    src = results / "concordance" / "fig2b_deg_venn.pdf"
    if src.exists():
        shutil.copy(src, out / "fig2b_deg_venn.pdf")
        LOG.info("fig2b: copied")
    else:
        LOG.warning("fig2b: source missing")


def fig3c(results: Path, out: Path) -> None:
    f = results / "scrnaseq" / "agreement" / "agreement_metrics.tsv"
    if not f.exists():
        LOG.warning("fig3c: %s missing", f); return
    df = pd.read_csv(f, sep="\t")
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    for metric, marker in [("ARI", "o"), ("NMI", "s"), ("ASW", "^")]:
        if metric in df.columns:
            ax.plot(df["resolution"], df[metric],
                    marker=marker, label=metric, linewidth=1.2)
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("Score")
    ax.set_ylim(0.7, 1.005)
    ax.axhline(0.947, linestyle="--", color="grey", linewidth=0.8,
               label="ARI = 0.947 at res = 1.0")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("CPU vs GPU clustering agreement")
    fig.tight_layout(); fig.savefig(out / "fig3c_agreement.pdf"); plt.close(fig)
    LOG.info("fig3c written")


def fig4ab(results: Path, out: Path) -> None:
    f = results / "benchmark" / "fig4_scaling.tsv"
    if not f.exists():
        LOG.warning("fig4ab: %s missing", f); return
    df = pd.read_csv(f, sep="\t")
    for mode, ax_letter in [("strong", "A"), ("weak", "B")]:
        sub = df[df["mode"] == mode]
        if sub.empty: continue
        fig, ax = plt.subplots(figsize=(4.0, 3.0))
        ax.plot(sub["nodes"], sub["wall_seconds"], "o-", linewidth=1.2)
        if mode == "strong":
            # ideal: T(N) = T(1) / N
            t1 = sub.loc[sub["nodes"].idxmin(), "wall_seconds"]
            n1 = sub["nodes"].min()
            ax.plot(sub["nodes"], t1 * n1 / sub["nodes"],
                    "--", color="grey", label="ideal")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("Compute nodes")
        ax.set_ylabel("Wall time (s)")
        ax.set_title(f"Fig 4{ax_letter}: {mode} scaling")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(out / f"fig4{ax_letter.lower()}_scaling_{mode}.pdf")
        plt.close(fig)
    LOG.info("fig4ab written")


def fig4d(results: Path, out: Path) -> None:
    f = results / "benchmark" / "fig4d_energy.tsv"
    if not f.exists():
        LOG.warning("fig4d: %s missing", f); return
    df = pd.read_csv(f, sep="\t")
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    grouped = df.groupby("backend")["energy_kJ"].agg(["mean", "std"])
    ax.bar(grouped.index, grouped["mean"], yerr=grouped["std"],
           capsize=4, color=["#4a90d9", "#d9534f"])
    ax.set_ylabel("Energy (kJ)")
    ax.set_title("Per-task energy: CPU vs GPU")
    fig.tight_layout()
    fig.savefig(out / "fig4d_energy.pdf")
    plt.close(fig)
    LOG.info("fig4d written")


def fig5d(results: Path, out: Path) -> None:
    src = results / "celltype_de" / "fig5d_epithelial_volcano.pdf"
    if src.exists():
        shutil.copy(src, out / "fig5d_epithelial_volcano.pdf")
        LOG.info("fig5d: copied")
    else:
        LOG.warning("fig5d: source missing")


def fig6(results: Path, out: Path) -> None:
    """KEGG dot-plot - picks the largest enrichment table available."""
    candidates = list(results.glob("enrichment/*/*/kegg_enrichment.tsv"))
    if not candidates:
        LOG.warning("fig6: no KEGG enrichment tables found"); return
    src = max(candidates, key=lambda p: p.stat().st_size)
    LOG.info("fig6: using %s", src)
    df = pd.read_csv(src, sep="\t").head(15)
    df["gene_ratio_num"] = df["GeneRatio"].apply(
        lambda x: eval(x) if isinstance(x, str) and "/" in x else float(x))
    df = df.sort_values("gene_ratio_num")
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(df["gene_ratio_num"], df["Description"],
               s=df["Count"] * 6,
               c=-np.log10(df["p.adjust"]),
               cmap="viridis")
    ax.set_xlabel("Gene ratio")
    ax.set_title("Top KEGG pathways")
    fig.tight_layout()
    fig.savefig(out / "fig6_kegg.pdf")
    plt.close(fig)
    LOG.info("fig6 written")


# ---------------------------------------------------------------------------
DISPATCH = {
    "fig1": fig1,
    "fig2": lambda r, o: (fig2a(r, o), fig2b(r, o)),
    "fig3": fig3c,
    "fig4": lambda r, o: (fig4ab(r, o), fig4d(r, o)),
    "fig5": fig5d,
    "fig6": fig6,
}


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--out",     required=True, type=Path)
    parser.add_argument("--figures", default="fig1,fig2,fig3,fig4,fig5,fig6")
    args = parser.parse_args(argv)

    np.random.seed(2025)
    args.out.mkdir(parents=True, exist_ok=True)

    for fig_key in args.figures.split(","):
        fn = DISPATCH.get(fig_key.strip())
        if fn is None:
            LOG.warning("Unknown figure key: %s", fig_key); continue
        try:
            fn(args.results, args.out)
        except Exception as exc:
            LOG.exception("Failed to render %s: %s", fig_key, exc)

    # Provenance
    prov = {
        "results": str(args.results.resolve()),
        "out":     str(args.out.resolve()),
        "figures": args.figures.split(","),
        "matplotlib": matplotlib.__version__,
        "numpy":     np.__version__,
        "pandas":    pd.__version__,
        "seed":      2025,
    }
    (args.out / "figures_provenance.json").write_text(
        json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
