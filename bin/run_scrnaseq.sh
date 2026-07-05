#!/usr/bin/env bash
# bin/run_scrnaseq.sh - run both CPU (Seurat) and GPU (rapids-singlecell) sc
# workflows in parallel and emit the clustering-agreement metrics for Fig 3C.
set -euo pipefail

PROFILE="${PROFILE:-slurm}"
OUTDIR="${OUTDIR:-results/scrnaseq}"
SAMPLES="${SAMPLES:-data/samples_GSE160269.csv}"
RUN_GPU="${RUN_GPU:-true}"

mkdir -p "${OUTDIR}"

nextflow run main.nf -profile "${PROFILE}" \
    --samplesheet "${SAMPLES}" \
    --workflows scrnaseq \
    --run_gpu "${RUN_GPU}" \
    --outdir "${OUTDIR}" \
    -with-report  "${OUTDIR}/report.html" \
    -with-trace   "${OUTDIR}/trace.tsv" \
    -with-timeline "${OUTDIR}/timeline.html"
