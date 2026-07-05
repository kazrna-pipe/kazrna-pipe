#!/usr/bin/env python3
"""scripts/python/write_provenance.py

Writes a per-task provenance JSON file. Used by tasks that wrap closed-source
tools (e.g. CIBERSORTx) where standard library introspection is unavailable.
"""

from __future__ import annotations

import argparse
import datetime as dt
import getpass
import json
import os
import platform
import socket
import subprocess
import sys
from pathlib import Path


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool",  required=True)
    parser.add_argument("--image", default=None)
    parser.add_argument("--out",   required=True, type=Path)
    parser.add_argument("--extra", action="append", default=[],
                        help="key=value pairs to add to the provenance record")
    args = parser.parse_args()

    record = {
        "tool":      args.tool,
        "image":     args.image,
        "host":      socket.gethostname(),
        "user":      getpass.getuser(),
        "platform":  platform.platform(),
        "python":    sys.version.split()[0],
        "cwd":       os.getcwd(),
        "git_commit": git_commit(),
        "timestamp": dt.datetime.now(dt.timezone.utc)
                       .isoformat(timespec="seconds"),
    }
    for kv in args.extra:
        if "=" in kv:
            k, v = kv.split("=", 1)
            record[k] = v

    args.out.write_text(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
