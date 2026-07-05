#!/usr/bin/env bash
# bin/prepare_references.sh
#
# Downloads GENCODE v44 GRCh38.p14 reference assets, verifies SHA-256 against
# data/snapshots/references/frozen_manifest.tsv, and builds the STAR, HISAT2,
# and Salmon indices used by the pipeline.
#
# Usage:
#   bash bin/prepare_references.sh --out_dir data/refs --threads 16
#
# The resulting `data/refs/` tree is referenced by `nextflow.config` defaults.

set -euo pipefail

# ----- defaults ----------------------------------------------------------------
OUT_DIR="data/refs"
THREADS=8
RELEASE=44
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}/GRCh38.p14.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}/gencode.v${RELEASE}.primary_assembly.annotation.gtf.gz"
TX_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}/gencode.v${RELEASE}.transcripts.fa.gz"
MANIFEST="data/snapshots/references/frozen_manifest.tsv"

# ----- argparse ----------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out_dir)  OUT_DIR="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --release)  RELEASE="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "${OUT_DIR}"/{star_index,hisat2_index,salmon_index}
cd "${OUT_DIR}"

# ----- download + verify -------------------------------------------------------
download_and_check () {
    local url="$1" expected="$2"
    local file; file=$(basename "$url")
    if [[ ! -f "$file" ]]; then
        echo "[$(date -Is)] download: $file"
        curl -sSL --retry 5 -o "$file" "$url"
    fi
    local got; got=$(sha256sum "$file" | awk '{print $1}')
    if [[ "$got" != "$expected" ]]; then
        echo "ERROR: SHA-256 mismatch for $file" >&2
        echo "  expected: $expected" >&2
        echo "  observed: $got"      >&2
        echo "  upstream may have re-versioned. Re-pin frozen_manifest.tsv if intentional." >&2
        exit 3
    fi
    echo "[$(date -Is)] verified: $file"
}

GENOME_SHA=$(awk -F'\t' '$1=="GRCh38.p14.genome.fa.gz"{print $4}'                          "../../${MANIFEST}")
GTF_SHA=$(   awk -F'\t' '$1=="gencode.v44.primary_assembly.annotation.gtf.gz"{print $4}'   "../../${MANIFEST}")
TX_SHA=$(    awk -F'\t' '$1=="gencode.v44.transcripts.fa.gz"{print $4}'                    "../../${MANIFEST}")

download_and_check "$GENOME_URL" "$GENOME_SHA"
download_and_check "$GTF_URL"    "$GTF_SHA"
download_and_check "$TX_URL"     "$TX_SHA"

gunzip -k -f GRCh38.p14.genome.fa.gz                              # → GRCh38.p14.genome.fa
gunzip -k -f gencode.v44.primary_assembly.annotation.gtf.gz       # → ...gtf
gunzip -k -f gencode.v44.transcripts.fa.gz                        # → gencode.v44.transcripts.fa

# ----- tx2gene -----------------------------------------------------------------
echo "[$(date -Is)] build tx2gene"
awk 'BEGIN{OFS="\t"}
     $3=="transcript" {
       match($0, /transcript_id "([^"]+)"/, t);
       match($0, /gene_id "([^"]+)"/,       g);
       match($0, /gene_name "([^"]+)"/,     n);
       print t[1], g[1], (n[1]?n[1]:g[1])
     }' gencode.v44.primary_assembly.annotation.gtf > tx2gene.v44.tsv

# ----- STAR --------------------------------------------------------------------
echo "[$(date -Is)] STAR genomeGenerate"
STAR --runThreadN "${THREADS}" \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.p14.genome.fa \
     --sjdbGTFfile      gencode.v44.primary_assembly.annotation.gtf \
     --sjdbOverhang     99

# ----- HISAT2 ------------------------------------------------------------------
echo "[$(date -Is)] hisat2-build"
hisat2-build -p "${THREADS}" GRCh38.p14.genome.fa hisat2_index/GRCh38_v44

# ----- Salmon (decoy-aware) ----------------------------------------------------
echo "[$(date -Is)] Salmon decoy + index"
grep '^>' GRCh38.p14.genome.fa | sed 's/^>//;s/ .*//' > decoys.txt
cat gencode.v44.transcripts.fa GRCh38.p14.genome.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -k 31 \
             -p "${THREADS}" -i salmon_index --gencode

echo "[$(date -Is)] DONE. Indices in ${OUT_DIR}"
