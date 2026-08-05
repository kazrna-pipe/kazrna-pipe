#!/usr/bin/env bash
# Builds the tiny self-consistent test fixtures used by -profile test.
# Uses the same biocontainers the pipeline pins, so index versions match.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TD="${ROOT}/test_data"
REFS="${TD}/refs"
FQ="${TD}/bulk_fastq"
SC="${TD}/sc_h5"
mkdir -p "${REFS}" "${FQ}" "${SC}"

echo ">> [1/7] Generating synthetic reference (200 genes, ~304 kb)"
python3 "${ROOT}/test_data/make_reference.py" "${REFS}"

echo ">> [2/7] Generating paired-end reads (4 samples, 40 DE genes)"
python3 "${ROOT}/test_data/make_reads.py" "${REFS}/mini.fa" "${FQ}"

echo ">> [3/7] Building STAR index"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/star:2.7.11a--h0033a41_0 \
    STAR --runMode genomeGenerate --genomeDir /refs/star_index \
         --genomeFastaFiles /refs/mini.fa --sjdbGTFfile /refs/gencode.v44.annotation.gtf \
         --genomeSAindexNbases 7 --sjdbOverhang 99 --runThreadN 2 || true
test -f "${REFS}/star_index/SAindex" || { echo "STAR index incomplete"; exit 1; }

echo ">> [4/7] Building HISAT2 index"
mkdir -p "${REFS}/hisat2_index"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/hisat2:2.2.1--h87f3376_5 \
    hisat2-build /refs/mini.fa /refs/hisat2_index/genome

echo ">> [5/7] Building Salmon index (transcriptome)"
python3 "${ROOT}/test_data/make_transcriptome.py" "${REFS}"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/salmon:1.10.2--hecfa306_0 \
    salmon index -t /refs/mini.transcripts.fa -i /refs/salmon_index -k 21

echo ">> [6/7] Generating tiny single-cell 10x .h5 matrices"
# Run in the pipeline's own python image rather than on the host: this script
# needs numpy, scipy and h5py, and requiring them on every developer's machine
# would make fixture regeneration depend on the host environment. The reference
# and read generators use only the standard library, so they run directly.
docker run --rm --user $(id -u):$(id -g) -v "${SC}:/sc" -v "${ROOT}/test_data:/src:ro" \
    ghcr.io/kazrna-pipe/kazrna-py@sha256:17bc730681ca90d649f1bf25840f4672ded381747338085e0e6645a1c63f6832 python3 /src/make_sc_h5.py /sc

echo ">> [7/7] Recording checksums"
cd "${TD}"
find . -type f ! -name SHA256SUMS ! -name README.md -print0 \
    | sort -z | xargs -0 sha256sum > SHA256SUMS
echo ">> Wrote ${TD}/SHA256SUMS ($(wc -l < SHA256SUMS) files)"

echo ">> Done. Fixtures written under ${TD}"