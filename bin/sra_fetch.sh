#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# sra_fetch.sh
#
# Resolves an NCBI BioProject (or SRA Study) accession to its constituent runs
# via the EUtils API, downloads them with SRA Toolkit, converts to paired-end
# FASTQ, and verifies SHA-256 against the frozen manifest.
#
# Used to fetch PRJNA608223 (22 Kazakh ESCC bulk RNA-seq) for KazRNA-Pipe.
#
# Usage:
#   bash bin/sra_fetch.sh --accession PRJNA608223 --target data/raw/PRJNA608223 --verify-sha256
# ---------------------------------------------------------------------------
set -euo pipefail

ACCESSION=""
TARGET=""
VERIFY_SHA256="false"
THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession)     ACCESSION="$2";       shift 2 ;;
        --target)        TARGET="$2";          shift 2 ;;
        --verify-sha256) VERIFY_SHA256="true"; shift   ;;
        --threads)       THREADS="$2";         shift 2 ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ----$/p' "$0"
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$ACCESSION" ]] && { echo "Missing --accession" >&2; exit 1; }
[[ -z "$TARGET"   ]] && { echo "Missing --target"   >&2; exit 1; }

mkdir -p "$TARGET"
cd "$TARGET"

# ---- Tool checks ----------------------------------------------------------
for t in prefetch fasterq-dump pigz curl jq sha256sum; do
    command -v "$t" >/dev/null 2>&1 || { echo "Missing tool: $t" >&2; exit 1; }
done

# ---- 1. Resolve the BioProject to its run accessions ---------------------
# Uses the EBI ENA REST API (more reliable than NCBI EUtils for bulk).
echo "[$(date -Iseconds)] Resolving $ACCESSION to run accessions..."
RUNS_TSV="${ACCESSION}.runs.tsv"
curl --silent --fail \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession,sample_alias,library_strategy,library_layout,read_count,base_count,fastq_md5,fastq_ftp&format=tsv" \
    > "$RUNS_TSV"

NUM_RUNS=$(($(wc -l < "$RUNS_TSV") - 1))
echo "[$(date -Iseconds)] Found ${NUM_RUNS} runs for ${ACCESSION}."

# Persist the snapshot date so the reviewer can verify upstream changes.
echo "ENA filereport snapshot for ${ACCESSION}" > snapshot.txt
echo "Captured at: $(date -Iseconds)"           >> snapshot.txt
echo "Endpoint:    https://www.ebi.ac.uk/ena/portal/api/filereport" >> snapshot.txt
echo "Run count:   ${NUM_RUNS}"                 >> snapshot.txt

# ---- 2. Download each run via prefetch + fasterq-dump --------------------
mkdir -p fastq
tail -n +2 "$RUNS_TSV" | while IFS=$'\t' read -r run_acc sample_alias strategy layout reads bases fastq_md5 fastq_ftp; do
    echo "[$(date -Iseconds)] -> $run_acc  (sample: $sample_alias, strategy: $strategy)"
    if [[ -f "fastq/${run_acc}_1.fastq.gz" && -f "fastq/${run_acc}_2.fastq.gz" ]]; then
        echo "    Already present, skipping."
        continue
    fi

    prefetch --max-size 100g "$run_acc" -O .

    # fasterq-dump produces split FASTQs in the current directory; pigz to gzip.
    fasterq-dump \
        --split-files \
        --threads "$THREADS" \
        --skip-technical \
        --outdir fastq \
        "${run_acc}/${run_acc}.sra"

    pigz -p "$THREADS" "fastq/${run_acc}_1.fastq" "fastq/${run_acc}_2.fastq"
    rm -rf "${run_acc}"
done

# ---- 3. Build a SHA-256 manifest ------------------------------------------
echo "[$(date -Iseconds)] Computing SHA-256 manifest..."
(cd fastq && sha256sum *.fastq.gz | sort) > sha256.manifest

if [[ "$VERIFY_SHA256" == "true" ]]; then
    FROZEN_MANIFEST="${0%/*}/../data/${ACCESSION}.sha256"
    if [[ -f "$FROZEN_MANIFEST" ]]; then
        if diff -q sha256.manifest "$FROZEN_MANIFEST" > /dev/null; then
            echo "[$(date -Iseconds)] SHA-256 manifest matches the frozen reference."
        else
            echo "WARNING: SHA-256 manifest differs from frozen reference."
            echo "         Differences saved to sha256.diff for inspection."
            diff sha256.manifest "$FROZEN_MANIFEST" > sha256.diff || true
            exit 2
        fi
    else
        echo "NOTE: No frozen manifest at $FROZEN_MANIFEST; cannot verify."
    fi
fi

echo "[$(date -Iseconds)] Done. Output: $(pwd)"
