#!/usr/bin/env python3
"""
check_interfaces.py - static audit of Nextflow module <-> analysis script contracts.

Catches the two failure classes that stop kazrna-pipe from completing a run,
without executing R, Python or Nextflow:

  1. CLI mismatch    - a .nf module passes a flag the script does not define
  2. OUTPUT mismatch - a .nf module declares an output the script never writes

Usage:
    python3 tools/check_interfaces.py            # from the repository root
    python3 tools/check_interfaces.py --quiet    # failures only
Exit 0 = every module consistent; exit 1 = at least one mismatch.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

G, R, Y, D, X = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"

INVOCATION = re.compile(
    r'(?:Rscript|python3?|bash)\s+\$\{projectDir\}/(scripts/[^\s\\"\']+)([\s\S]*?)'
    r'(?=(?:Rscript|python3?|bash)\s+\$\{projectDir\}|\n\s*cat\s*<<|\n\s*"""|\Z)'
)


def script_options(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix == ".R":
        return set(re.findall(r'make_option\(\s*"--([A-Za-z0-9_.]+)"', text))
    if path.suffix == ".py":
        return set(re.findall(r'add_argument\(\s*"--([A-Za-z0-9_-]+)"', text))
    if path.suffix == ".sh":
        return set(re.findall(r'--([A-Za-z0-9_-]+)\)', text))
    return set()


def script_writes(path: Path) -> tuple[set[str], set[str], bool]:
    """(literal filenames, --out_prefix suffixes, filenames_come_from_arguments)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    literals = set(re.findall(
        r'["\']([A-Za-z0-9_.\-]+\.(?:tsv|csv|txt|pdf|png|json|rds|h5ad|yml|yaml|h5|log))["\']', text))
    suffixes = set(re.findall(r'["\'](_[A-Za-z0-9_.\-]+\.(?:tsv|csv|txt|pdf|json|rds))["\']', text))
    # If no literal filename appears, the script is told where to write via its
    # arguments (--out_yml, --out) and cannot be audited statically.
    arg_driven = (not literals) or path.suffix == ".sh"
    return literals, suffixes, arg_driven


def parse_module(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")

    body = text.split("script:", 1)[-1]
    calls = [(m.group(1), set(re.findall(r'--([A-Za-z0-9_-]+)', m.group(2))))
             for m in INVOCATION.finditer(body)]

    out_block = ""
    om = re.search(r'\boutput:\s*(.*?)(?:\n\s*(?:script|stub|shell|when):)', text, re.S)
    if om:
        out_block = om.group(1)
    outputs = set(re.findall(r'path\s*\(?\s*["\']([^"\'${}]+)["\']', out_block))

    pm = re.search(r'--out_prefix\s+([A-Za-z0-9_.\-]+)', body)
    return calls, outputs, (pm.group(1) if pm else None)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--root", default=".")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    modules = sorted((root / "modules").glob("*.nf"))
    if not modules:
        print(f"{R}No modules/*.nf under {root}{X}")
        return 1

    failures = warnings = 0
    for mod in modules:
        calls, outputs, prefix = parse_module(mod)
        name = mod.name

        if not calls:
            if not args.quiet:
                print(f"{D}SKIP {name:26} no in-repo script invoked{X}")
            continue

        problems, notes = [], []

        for script_rel, flags in calls:
            script = root / script_rel
            if not script.exists():
                problems.append(f"  MISSING invokes {script_rel}, which does not exist")
                continue
            accepted = script_options(script)
            if not accepted:
                continue
            unknown = sorted(f for f in flags if f not in accepted)
            if unknown:
                problems.append(
                    f"  CLI     {script_rel}\n"
                    f"          passes  {', '.join('--' + u for u in unknown)}\n"
                    f"          defines {', '.join('--' + a for a in sorted(accepted))}")

        main_script = root / calls[0][0]
        if main_script.exists():
            literals, suffixes, arg_driven = script_writes(main_script)
            missing = [o for o in sorted(outputs)
                       if o != "versions.yml" and "*" not in o and o not in literals
                       and not (prefix and any(o == prefix + s for s in suffixes))]
            if missing and arg_driven:
                notes.append(f"  ?       declares {', '.join(missing)}; script names outputs "
                             f"from its arguments, so verify by hand")
            elif missing:
                problems.append(
                    f"  OUTPUT  declares {', '.join(missing)}\n"
                    f"          script writes {', '.join(sorted(literals)) or '(nothing detected)'}")

        if problems:
            failures += 1
            print(f"{R}FAIL{X} {name:26} {D}-> {calls[0][0]}{X}")
            print("\n".join(problems))
        elif notes:
            warnings += 1
            if not args.quiet:
                print(f"{Y}WARN{X} {name:26} {D}-> {calls[0][0]}{X}")
                print("\n".join(notes))
        elif not args.quiet:
            print(f"{G}OK  {X} {name:26} {D}-> {calls[0][0]}{X}")

    print()
    if failures:
        print(f"{R}{failures} module(s) inconsistent{X}" +
              (f"{Y}, {warnings} to verify by hand{X}" if warnings else ""))
        return 1
    print(f"{G}All modules consistent{X}" +
          (f"{Y} ({warnings} to verify by hand){X}" if warnings else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
