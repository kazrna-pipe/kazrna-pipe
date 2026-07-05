#!/usr/bin/env python3
"""scripts/python/collect_versions.py

Merge per-process versions.yml files into one canonical report (YAML, TSV,
HTML). Used by modules/software_versions.nf at the end of every pipeline run.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from collections import defaultdict
from pathlib import Path

import yaml


def parse(paths: list[Path]) -> dict[str, dict[str, str]]:
    merged: dict[str, dict[str, str]] = defaultdict(dict)
    for p in paths:
        try:
            data = yaml.safe_load(p.read_text())
        except yaml.YAMLError as exc:
            print(f"WARNING: failed to parse {p}: {exc}", file=sys.stderr)
            continue
        if not isinstance(data, dict):
            continue
        for process, tools in data.items():
            if isinstance(tools, dict):
                for tool, version in tools.items():
                    merged[process][tool] = str(version).strip()
    return merged


def write_yml(out: Path, data: dict) -> None:
    out.write_text(yaml.safe_dump(dict(sorted(data.items())),
                                  sort_keys=True,
                                  default_flow_style=False))


def write_tsv(out: Path, data: dict) -> None:
    lines = ["process\ttool\tversion"]
    for process in sorted(data):
        for tool, ver in sorted(data[process].items()):
            lines.append(f"{process}\t{tool}\t{ver}")
    out.write_text("\n".join(lines) + "\n")


def write_html(out: Path, data: dict) -> None:
    rows = []
    for process in sorted(data):
        for tool, ver in sorted(data[process].items()):
            rows.append(f"<tr><td>{process}</td><td>{tool}</td><td><code>{ver}</code></td></tr>")
    html = f"""<!DOCTYPE html>
<html><head><title>kazrna-pipe software versions</title>
<style>
body {{ font-family: -apple-system, sans-serif; max-width: 900px; margin: 2em auto; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ text-align: left; padding: .4em .7em; border-bottom: 1px solid #ddd; }}
th {{ background: #f4f4f4; }}
code {{ background: #f4f4f4; padding: 1px 4px; border-radius: 3px; }}
</style></head>
<body>
<h1>kazrna-pipe - software versions</h1>
<p>Generated {dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}</p>
<table>
<thead><tr><th>Process</th><th>Tool</th><th>Version</th></tr></thead>
<tbody>
{chr(10).join(rows)}
</tbody></table>
</body></html>
"""
    out.write_text(html)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs",   nargs="+", required=True, type=Path)
    parser.add_argument("--out_yml",  required=True, type=Path)
    parser.add_argument("--out_tsv",  required=True, type=Path)
    parser.add_argument("--out_html", required=True, type=Path)
    args = parser.parse_args(argv)

    data = parse(args.inputs)
    write_yml(args.out_yml, data)
    write_tsv(args.out_tsv, data)
    write_html(args.out_html, data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
