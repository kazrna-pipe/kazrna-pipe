#!/usr/bin/env bash
# bin/reproduce_figures.sh - creating figures from a
# completed pipeline output. Useful when iterating on plotting only.
#
# Usage:
#   bash bin/reproduce_figures.sh --results results/ --out figures/
set -euo pipefail

RESULTS="results"
OUT="figures"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results) RESULTS="$2"; shift 2 ;;
        --out)     OUT="$2";     shift 2 ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "${OUT}"

python scripts/python/make_figures.py \
    --results "${RESULTS}" \
    --out     "${OUT}" \
    --figures fig1,fig2,fig3,fig4,fig5,fig6

echo "[$(date -Is)] Figures written to ${OUT}"
