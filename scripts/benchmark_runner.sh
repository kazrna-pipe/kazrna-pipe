#!/usr/bin/env bash
# scripts/benchmark_runner.sh
#
# Runs one scaling-benchmark point and emits two TSVs:
#   timings.tsv   nodes  mode  wall_seconds  cpu_seconds  peak_rss_GB
#   energy.tsv    nodes  mode  backend       energy_kJ    samples
#
# Wraps:
#   - Nextflow main.nf restricted to the bulk_rnaseq workflow on a representative
#     input, with `process.maxForks` scaled to *nodes*.
#   - perf-stat for CPU package RAPL counters when present.
#   - nvidia-smi sampling at 1 Hz for GPU power.

set -euo pipefail

NODES=1
MODE="strong"
INPUT=""
OUT_DIR="."
LABEL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --nodes)   NODES="$2"; shift 2 ;;
        --mode)    MODE="$2";  shift 2 ;;
        --input)   INPUT="$2"; shift 2 ;;
        --out_dir) OUT_DIR="$2"; shift 2 ;;
        --label)   LABEL="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# Sample GPU power in the background, if we have GPUs
NVS_PID=""
if command -v nvidia-smi >/dev/null 2>&1; then
    (nvidia-smi --query-gpu=timestamp,index,power.draw,memory.used \
                --format=csv,noheader,nounits \
                -l 1 > "nvidia_smi.log") &
    NVS_PID=$!
fi

# RAPL CPU energy (sudo-less via perf if permitted)
RAPL_PID=""
if perf list 2>/dev/null | grep -q "power/energy-pkg/"; then
    (perf stat -a -e power/energy-pkg/ -I 1000 -- sleep 100000 \
        2> "rapl_${LABEL}.log") &
    RAPL_PID=$!
fi

START=$(date +%s)

# Cap parallelism to the nodes-equivalent number of forks. For "weak", scale
# input size with N; for "strong", keep input constant.
case "${MODE}" in
    weak)   FORKS=$((NODES * 4)) ;;
    strong) FORKS=$((NODES * 4)) ;;
    *) echo "Unknown mode: ${MODE}"; exit 2 ;;
esac

nextflow run "$(dirname "$0")/../main.nf" \
    -profile slurm \
    --benchmark_only true \
    --samplesheet "${INPUT}" \
    --max_cpus 64 \
    --max_memory 256.GB \
    --outdir "ci_bench_${LABEL}" \
    -ansi-log false \
    -process.maxForks "${FORKS}" \
    > "nextflow_${LABEL}.log" 2>&1

END=$(date +%s)
WALL=$((END - START))

# Stop the samplers
if [[ -n "${NVS_PID}" ]]; then kill "${NVS_PID}"  2>/dev/null || true; fi
if [[ -n "${RAPL_PID}" ]]; then kill "${RAPL_PID}" 2>/dev/null || true; fi
wait 2>/dev/null || true

# Process samplers
gpu_energy_kJ=0
if [[ -s "nvidia_smi.log" ]]; then
    gpu_energy_kJ=$(awk -F, '{s += $3} END {print s/1000}' nvidia_smi.log)
fi
cpu_energy_kJ=0
if [[ -s "rapl_${LABEL}.log" ]]; then
    cpu_energy_kJ=$(awk '/energy-pkg/ {s+=$2} END {print s/1000}' "rapl_${LABEL}.log")
fi

# Emit timings.tsv (header + one row)
{
    echo -e "nodes\tmode\twall_seconds\tcpu_energy_kJ\tgpu_energy_kJ\tlabel"
    echo -e "${NODES}\t${MODE}\t${WALL}\t${cpu_energy_kJ}\t${gpu_energy_kJ}\t${LABEL}"
} > timings.tsv

# Emit energy.tsv (CPU and GPU as separate rows so plotting code can group)
{
    echo -e "nodes\tmode\tbackend\tenergy_kJ\tlabel"
    echo -e "${NODES}\t${MODE}\tCPU\t${cpu_energy_kJ}\t${LABEL}"
    echo -e "${NODES}\t${MODE}\tGPU\t${gpu_energy_kJ}\t${LABEL}"
} > energy.tsv

echo "[$(date -Is)] benchmark point done: ${LABEL}, wall=${WALL}s, " \
     "CPU=${cpu_energy_kJ} kJ, GPU=${gpu_energy_kJ} kJ"
