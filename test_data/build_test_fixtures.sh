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

echo ">> [1/6] Generating a small synthetic reference (3 genes, ~30 kb)"
python3 "${ROOT}/test_data/make_reference.py" "${REFS}"

echo ">> [2/6] Generating matched paired-end reads (4 samples)"
python3 "${ROOT}/test_data/make_reads.py" "${REFS}/mini.fa" "${FQ}"

echo ">> [3/6] Building STAR index"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/star:2.7.11a--h0033a41_0 \
    STAR --runMode genomeGenerate --genomeDir /refs/star_index \
         --genomeFastaFiles /refs/mini.fa --sjdbGTFfile /refs/gencode.v44.annotation.gtf \
         --genomeSAindexNbases 5 --sjdbOverhang 99 --runThreadN 2 || true
test -f "${REFS}/star_index/SAindex" || { echo "STAR index incomplete"; exit 1; }

echo ">> [4/6] Building HISAT2 index"
mkdir -p "${REFS}/hisat2_index"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/hisat2:2.2.1--h87f3376_5 \
    hisat2-build /refs/mini.fa /refs/hisat2_index/genome

echo ">> [5/6] Building Salmon index (transcriptome)"
python3 "${ROOT}/test_data/make_transcriptome.py" "${REFS}"
docker run --rm --user $(id -u):$(id -g) -w /refs -v "${REFS}:/refs" quay.io/biocontainers/salmon:1.10.2--hecfa306_0 \
    salmon index -t /refs/mini.transcripts.fa -i /refs/salmon_index -k 21

echo ">> [6/6] Generating tiny single-cell 10x .h5 matrices"
python3 "${ROOT}/test_data/make_sc_h5.py" "${SC}"

echo ">> Done. Fixtures written under ${TD}"
