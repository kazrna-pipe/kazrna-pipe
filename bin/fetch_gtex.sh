#!/usr/bin/env bash
# bin/fetch_gtex.sh - pull GTEx v8 esophagus mucosa gene-reads file and verify.
set -euo pipefail

OUT_DIR="data/raw/GTEx-v8"
CHECK="data/snapshots/GTEx-v8/frozen_manifest.tsv"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)   OUT_DIR="$2"; shift 2 ;;
        --check) CHECK="$2";   shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

while IFS=$'\t' read -r file url bytes sha256; do
    [[ "$file" == "file" || "$file" == \#* ]] && continue
    if [[ ! -f "$file" ]]; then
        echo "[$(date -Is)] fetch: $file"
        curl -sSL --retry 5 -o "$file" "$url"
    fi
    got=$(sha256sum "$file" | awk '{print $1}')
    if [[ "$got" != "$sha256" ]]; then
        echo "ERROR: SHA-256 mismatch for $file ($got vs $sha256)" >&2
        exit 3
    fi
    echo "verified: $file"
done < "../../../${CHECK}"
