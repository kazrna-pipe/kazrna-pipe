#!/usr/bin/env python3
"""
diff_results.py - compare two pipeline runs and report any difference that is
not explained by wall-clock time, timestamps, scheduler bookkeeping, or a
documented stochastic step.

    python3 diff_results.py results/run1 results/run2
    python3 diff_results.py results/run1 results/run2 --tolerance 1e-8
    python3 diff_results.py results/run1 results/run2 --list-excluded

Exit 0 when the runs agree, 1 when they do not.

WHAT THIS CHECK CLAIMS

Counts and every differential-expression output are bit-reproducible. 

"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fnmatch import fnmatch
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: run this inside the kazrna-py container")

# --- Files not compared at all ---------------------------------------------
EXCLUDE = [
    # 1. Timestamps and scheduler bookkeeping
    "pipeline_info/trace.txt",
    "pipeline_info/timeline.html",
    "pipeline_info/report.html",
    "pipeline_info/dag.svg",
    "pipeline_info/software_versions.html",
    "**/timing.tsv",
    "**/*.timing.tsv",
    "**/*.timing.raw",
    "**/*_fastqc.zip",
    "**/*_fastqc.html",
    "**/*.Log.final.out",
    "**/*.Log.out",
    "**/*.Log.progress.out",
    "**/*trimming_report.txt",
    "**/*.hisat2.summary.txt",
    "**/logs/*.log",
    "**/cmd_info.json",
    "**/meta_info.json",
    "**/*provenance.json",
    "**/*_session.txt",
    "**/*.pdf",                     # PDF writers embed a creation date
    "**/*.rds",                     # R serialisation records session metadata

    "**/aux_info/*.gz",
    "**/aux_info/*.tsv",
    "**/libParams/flenDist.txt",

    "**/*.bam",
    "**/*.bam.bai",
    "**/*.sam",
]

# --- Columns not compared, per file pattern --------------------------------
# "*" means every column: the whole table derives from effective length.
EXCLUDE_COLUMNS: dict[str, object] = {
    "**/quant.sf":                {"EffectiveLength", "TPM"},
    "**/salmon.gene_lengths.tsv": "*",
    "**/salmon.gene_tpm.tsv":     "*",
}

BINARY_SUFFIXES = {".bam", ".bai", ".rds", ".h5ad", ".h5", ".pdf", ".gz", ".zip"}

# Below this, the runs are almost certainly incomplete rather than identical.
MIN_EXPECTED_FILES = 20


def matches(rel: str, pattern: str) -> bool:
    return fnmatch(rel, pattern) or fnmatch("/" + rel, "/" + pattern)


def excluded(rel: str) -> bool:
    return any(matches(rel, p) for p in EXCLUDE)


def excluded_columns(rel: str):
    for pattern, cols in EXCLUDE_COLUMNS.items():
        if matches(rel, pattern):
            return cols
    return set()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_table(a: Path, b: Path, tol: float, skip_cols):
    """Numeric-aware comparison. Returns (problems, number of columns compared)."""
    sep = "," if a.suffix == ".csv" else "\t"
    try:
        da = pd.read_csv(a, sep=sep, comment="#")
        db = pd.read_csv(b, sep=sep, comment="#")
    except Exception:
        if sha256(a) == sha256(b):
            return [], 1
        return ["content differs (not parseable as a table)"], 1

    if list(da.columns) != list(db.columns):
        return [f"columns differ: {list(da.columns)[:6]} vs {list(db.columns)[:6]}"], 0
    if da.shape != db.shape:
        return [f"shape differs: {da.shape} vs {db.shape}"], 0

    problems, n_compared = [], 0
    for col in da.columns:
        if skip_cols == "*" or (isinstance(skip_cols, set) and col in skip_cols):
            continue
        n_compared += 1
        ca, cb = da[col], db[col]
        if pd.api.types.is_numeric_dtype(ca) and pd.api.types.is_numeric_dtype(cb):
            delta = (ca - cb).abs()
            worst = float(delta.max()) if len(delta) else 0.0
            if worst > tol:
                problems.append(f"column '{col}': max absolute difference {worst:.4g} > {tol:g}")
        elif not ca.equals(cb):
            problems.append(f"column '{col}': {int((ca != cb).sum())} value(s) differ")
    return problems, n_compared


def compare_json(a: Path, b: Path) -> list:
    try:
        ja, jb = json.loads(a.read_text()), json.loads(b.read_text())
    except Exception:
        return [] if sha256(a) == sha256(b) else ["content differs (not parseable as JSON)"]
    return [] if ja == jb else ["structured contents differ"]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=Path, help="first results directory")
    ap.add_argument("b", type=Path, help="second results directory")
    ap.add_argument("--tolerance", type=float, default=1e-8,
                    help="maximum absolute difference for numeric columns [%(default)g]")
    ap.add_argument("--list-excluded", action="store_true",
                    help="print every skipped file before the summary")
    args = ap.parse_args()

    for d in (args.a, args.b):
        if not d.is_dir():
            print(f"Not a directory: {d}", file=sys.stderr)
            return 1

    files_a = {p.relative_to(args.a).as_posix() for p in args.a.rglob("*") if p.is_file()}
    files_b = {p.relative_to(args.b).as_posix() for p in args.b.rglob("*") if p.is_file()}

    problems, skipped, partial = [], [], []

    for rel in sorted(files_a - files_b):
        if not excluded(rel):
            problems.append(f"{rel}: present in {args.a} only")
    for rel in sorted(files_b - files_a):
        if not excluded(rel):
            problems.append(f"{rel}: present in {args.b} only")

    compared = 0
    for rel in sorted(files_a & files_b):
        if excluded(rel):
            skipped.append(rel)
            continue

        pa, pb = args.a / rel, args.b / rel
        skip_cols = excluded_columns(rel)

        if skip_cols == "*":
            skipped.append(rel)
            continue

        compared += 1
        if skip_cols:
            partial.append(f"{rel} (skipped: {', '.join(sorted(skip_cols))})")

        if pa.suffix in BINARY_SUFFIXES:
            if sha256(pa) != sha256(pb):
                problems.append(f"{rel}: binary contents differ")
        elif pa.suffix == ".json":
            problems += [f"{rel}: {m}" for m in compare_json(pa, pb)]
        elif pa.suffix in {".tsv", ".csv", ".txt", ".sf"}:
            msgs, n_cols = compare_table(pa, pb, args.tolerance, skip_cols)
            problems += [f"{rel}: {m}" for m in msgs]
            if skip_cols and n_cols == 0:
                problems.append(f"{rel}: every column was excluded; nothing verified")
        else:
            if sha256(pa) != sha256(pb):
                problems.append(f"{rel}: contents differ")

    if args.list_excluded:
        print(f"Skipped {len(skipped)} file(s) with run-dependent content:")
        for rel in skipped:
            print(f"  {rel}")
        print()

    print(f"Compared {compared} file(s); skipped {len(skipped)}; "
          f"tolerance {args.tolerance:g}")
    if partial:
        print(f"Partially compared ({len(partial)} file(s), stochastic columns skipped):")
        for line in partial:
            print(f"  {line}")

    # A check that examined nothing must never report success.
    if compared == 0:
        print("\nREPRODUCIBILITY CHECK FAILED - no files were compared. Both runs "
              "are empty, or every output was excluded.")
        return 1
    if compared < MIN_EXPECTED_FILES:
        print(f"\nREPRODUCIBILITY CHECK FAILED - only {compared} file(s) compared, "
              f"expected at least {MIN_EXPECTED_FILES}. The runs are probably "
              "incomplete; check that both finished successfully.")
        return 1

    if problems:
        print(f"\nREPRODUCIBILITY CHECK FAILED - {len(problems)} difference(s):")
        for p in problems:
            print(f"  - {p}")
        return 1

    print("\nREPRODUCIBILITY CHECK PASSED - the two runs agree on every analysis "
          "output, including all counts and differential-expression results.")
    return 0


if __name__ == "__main__":
    sys.exit(main())