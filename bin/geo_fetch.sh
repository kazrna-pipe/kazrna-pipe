#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# geo_fetch.sh
#
# Download the filtered count matrices for a GEO Series. Used to fetch
# GSE160269 (208,659-cell ESCC single-cell atlas) for KazRNA-Pipe.
# Usage:
#   bash bin/geo_fetch.sh --accession GSE160269 --target data/raw/GSE160269 --verify-sha256
# ---------------------------------------------------------------------------
set -euo pipefail

ACCESSION=""
TARGET=""
VERIFY_SHA256="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession)     ACCESSION="$2";       shift 2 ;;
        --target)        TARGET="$2";          shift 2 ;;
        --verify-sha256) VERIFY_SHA256="true"; shift   ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ----$/p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$ACCESSION" ]] && { echo "Missing --accession" >&2; exit 1; }
[[ -z "$TARGET"   ]] && { echo "Missing --target"   >&2; exit 1; }

mkdir -p "$TARGET"/{matrices,snapshots}
cd "$TARGET"

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/${ACCESSION:0:$((${#ACCESSION}-3))}nnn/${ACCESSION}"

# ---- 1. Mirror the metadata snapshot --------------------------------------
echo "[$(date -Iseconds)] Mirroring metadata for ${ACCESSION}..."
curl --silent --fail "${GEO_BASE}/matrix/${ACCESSION}_series_matrix.txt.gz" \
    -o "snapshots/${ACCESSION}_series_matrix_$(date +%F).txt.gz"
curl --silent --fail "${GEO_BASE}/soft/${ACCESSION}_family.soft.gz" \
    -o "snapshots/${ACCESSION}_family_$(date +%F).soft.gz" || true

echo "GEO metadata snapshot for ${ACCESSION}" > snapshot.txt
echo "Captured at: $(date -Iseconds)"        >> snapshot.txt
echo "Source:      ${GEO_BASE}"              >> snapshot.txt

# ---- 2. List the supplementary files --------------------------------------
echo "[$(date -Iseconds)] Listing supplementary files..."
curl --silent --fail "${GEO_BASE}/suppl/" > supp_listing.html
grep -oE "${ACCESSION}_[^\"]+\.(tar\.gz|h5|tar)" supp_listing.html \
    | sort -u > supp_files.txt
NUM_FILES=$(wc -l < supp_files.txt)
echo "[$(date -Iseconds)] ${NUM_FILES} supplementary files found."

# ---- 3. Download each supplementary file ----------------------------------
while IFS= read -r fname; do
    if [[ -f "matrices/${fname}" ]]; then
        echo "    ${fname} already present, skipping."
        continue
    fi
    echo "    -> ${fname}"
    curl --silent --fail --remote-time \
        "${GEO_BASE}/suppl/${fname}" \
        -o "matrices/${fname}"
done < supp_files.txt

# ---- 4. Unpack the per-patient matrices -----------------------------------
echo "[$(date -Iseconds)] Unpacking..."
(cd matrices && for f in *.tar.gz; do
    [[ -f "$f" ]] && tar -xzf "$f"
done)

# ---- 5. SHA-256 manifest --------------------------------------------------
(cd matrices && find . -type f -name "*.gz" -o -name "*.h5" -o -name "*.mtx" \
    -o -name "barcodes.tsv" -o -name "features.tsv" -o -name "genes.tsv" \
    | sort | xargs sha256sum) > sha256.manifest 2>/dev/null

if [[ "$VERIFY_SHA256" == "true" ]]; then
    FROZEN="${0%/*}/../data/${ACCESSION}.sha256"
    if [[ -f "$FROZEN" ]] && ! diff -q sha256.manifest "$FROZEN" > /dev/null; then
        echo "WARNING: SHA-256 manifest differs from frozen reference."
        diff sha256.manifest "$FROZEN" > sha256.diff || true
        exit 2
    fi
fi

echo "[$(date -Iseconds)] Done. Output: $(pwd)"
