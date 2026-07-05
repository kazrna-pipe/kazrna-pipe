#!/usr/bin/env python3
"""scripts/python/diff_results.py

Compares two output trees produced by two runs of the same pipeline and
reports any non-trivial differences. Exit code 0 iff the two runs agree.

Tolerances:
  - .tsv / .csv: numeric columns compared with rtol=1e-5, atol=1e-7
                 integer columns must match exactly
  - .yml / .yaml: parsed and compared as dicts (key order ignored)
  - .json: parsed and compared as dicts (key order ignored, timestamps allowed
           to differ)
  - .h5 / .h5ad / .rds / .bam / .pdf: existence-only check
  - everything else: SHA-256 must match

The CI uses this script to enforce that two consecutive `nextflow run` calls
produce byte-equivalent (or numerically equivalent) outputs. Stochastic stages
seed their RNGs (seed = 2025 throughout); any failure here is a real
reproducibility regression and should be investigated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import yaml
import numpy as np
import pandas as pd

SKIP_SUFFIXES = {".log", ".html"}
EXIST_ONLY    = {".h5", ".h5ad", ".rds", ".bam", ".bai", ".pdf", ".svg", ".png"}
TIMESTAMP_KEYS = {"timestamp", "captured_at", "built_at"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def diff_table(a: Path, b: Path, sep: str) -> list[str]:
    try:
        da = pd.read_csv(a, sep=sep, low_memory=False)
        db = pd.read_csv(b, sep=sep, low_memory=False)
    except Exception as exc:
        return [f"failed to parse {a.name} ({exc})"]
    if list(da.columns) != list(db.columns):
        return [f"columns differ in {a.name}"]
    if da.shape != db.shape:
        return [f"shape differs in {a.name}: {da.shape} vs {db.shape}"]
    diffs: list[str] = []
    for col in da.columns:
        if da[col].dtype.kind in "fc":
            if not np.allclose(da[col], db[col],
                               rtol=1e-5, atol=1e-7, equal_nan=True):
                diffs.append(f"{a.name}:{col} (float) differs beyond tol")
        elif da[col].dtype.kind in "iu":
            if not (da[col] == db[col]).all():
                diffs.append(f"{a.name}:{col} (int) differs")
        else:
            if not (da[col].astype(str) == db[col].astype(str)).all():
                diffs.append(f"{a.name}:{col} (str) differs")
    return diffs


def diff_structured(a: Path, b: Path, loader) -> list[str]:
    try:
        da = loader(a.read_text())
        db = loader(b.read_text())
    except Exception as exc:
        return [f"failed to parse {a.name} ({exc})"]
    if not isinstance(da, dict) or not isinstance(db, dict):
        return [] if da == db else [f"{a.name} contents differ"]
    for k in list(da.keys()) + list(db.keys()):
        if k in TIMESTAMP_KEYS:
            da.pop(k, None); db.pop(k, None)
    if da != db:
        return [f"{a.name} structured contents differ"]
    return []


def walk(root: Path) -> dict[Path, Path]:
    return {p.relative_to(root): p for p in root.rglob("*") if p.is_file()}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("a", type=Path)
    parser.add_argument("b", type=Path)
    args = parser.parse_args(argv)

    map_a = walk(args.a)
    map_b = walk(args.b)
    keys_a = set(map_a); keys_b = set(map_b)

    diffs: list[str] = []
    for missing in (keys_a - keys_b):
        diffs.append(f"only in A: {missing}")
    for extra in (keys_b - keys_a):
        diffs.append(f"only in B: {extra}")

    for rel in sorted(keys_a & keys_b):
        a = map_a[rel]; b = map_b[rel]
        suf = a.suffix.lower()
        if suf in SKIP_SUFFIXES:
            continue
        if suf in EXIST_ONLY:
            continue
        if suf in {".tsv", ".txt"}:
            diffs.extend(diff_table(a, b, sep="\t"))
            continue
        if suf == ".csv":
            diffs.extend(diff_table(a, b, sep=","))
            continue
        if suf in {".yml", ".yaml"}:
            diffs.extend(diff_structured(a, b, yaml.safe_load))
            continue
        if suf == ".json":
            diffs.extend(diff_structured(a, b, json.loads))
            continue
        # Fall back to byte comparison
        if sha256(a) != sha256(b):
            diffs.append(f"sha256 differs: {rel}")

    if diffs:
        print("REPRODUCIBILITY CHECK FAILED", file=sys.stderr)
        for d in diffs:
            print("  -", d, file=sys.stderr)
        return 1
    print("REPRODUCIBILITY CHECK OK ({} files compared)".format(len(keys_a & keys_b)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
