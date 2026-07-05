#!/usr/bin/env bash
# bin/run_bulk.sh - convenience wrapper for the bulk RNA-seq sub-pipeline.
# Logs the exact invocation and Nextflow version into the output directory.
set -euo pipefail

PROFILE="${PROFILE:-slurm}"
OUTDIR="${OUTDIR:-results/bulk}"
SAMPLES="${SAMPLES:-data/samples_PRJNA608223.csv}"

mkdir -p "${OUTDIR}"
{
    echo "# Run started: $(date -Is)"
    echo "# Host:        $(hostname)"
    echo "# User:        ${USER}"
    echo "# Nextflow:    $(nextflow -v 2>&1 | head -1)"
    echo "# Git commit:  $(git -C "$(dirname "$0")/.." rev-parse HEAD)"
    echo "# Command:"
    echo "nextflow run main.nf -profile ${PROFILE} \\"
    echo "    --samplesheet ${SAMPLES} \\"
    echo "    --workflows bulk \\"
    echo "    --outdir ${OUTDIR}"
} > "${OUTDIR}/.invocation.log"

exec nextflow run main.nf -profile "${PROFILE}" \
    --samplesheet "${SAMPLES}" \
    --workflows bulk \
    --outdir "${OUTDIR}" \
    -with-report  "${OUTDIR}/report.html" \
    -with-trace   "${OUTDIR}/trace.tsv" \
    -with-timeline "${OUTDIR}/timeline.html" \
    -with-dag      "${OUTDIR}/dag.svg"
