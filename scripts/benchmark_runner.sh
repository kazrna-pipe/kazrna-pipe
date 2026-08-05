#!/usr/bin/env bash
# scripts/benchmark_runner.sh
#
# Runs one scaling-benchmark point and emits two TSVs:
#   timings.tsv   nodes  mode  replicate  wall_seconds  cpu_energy_kJ  gpu_energy_kJ  n_samples  label
#   energy.tsv    nodes  mode  replicate  backend       energy_kJ      label
#
# Each point is repeated --replicates times so that mean and standard
# deviation can be computed from the raw measurements rather than asserted.
#
# strong scaling: total work is held constant while parallelism grows with N.
# weak scaling  : work grows in proportion to N, so per-node load is constant.
#
# Wraps:
#   - Nextflow main.nf on a representative input, with process.maxForks
#     scaled to *nodes*.
#   - perf-stat for CPU package RAPL counters when present.
#   - nvidia-smi sampling at 1 Hz for GPU power.

set -euo pipefail

NODES=1
MODE="strong"
INPUT=""
OUT_DIR="."
LABEL=""
REPLICATES=3
BASE_SAMPLES=8          # strong-scaling base problem size, in samples
CORES_PER_NODE=16
BENCH_PROFILE="slurm"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --nodes)          NODES="$2";          shift 2 ;;
        --mode)           MODE="$2";           shift 2 ;;
        --input)          INPUT="$2";          shift 2 ;;
        --out_dir)        OUT_DIR="$2";        shift 2 ;;
        --label)          LABEL="$2";          shift 2 ;;
        --replicates)     REPLICATES="$2";     shift 2 ;;
        --base_samples)   BASE_SAMPLES="$2";   shift 2 ;;
        --cores_per_node) CORES_PER_NODE="$2"; shift 2 ;;
        --profile)        BENCH_PROFILE="$2";  shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -n "${INPUT}" ]] || { echo "--input is required" >&2; exit 2; }
[[ -s "${INPUT}"  ]] || { echo "Input samplesheet not found: ${INPUT}" >&2; exit 2; }

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="$(cd "$(dirname "${INPUT}")" && pwd)/$(basename "${INPUT}")"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ---- Workload definition ---------------------------------------------------
# FORKS is the parallelism ceiling; SAMPLES is the size of the problem.
case "${MODE}" in
    strong) FORKS=$((NODES * CORES_PER_NODE)); SAMPLES=${BASE_SAMPLES} ;;
    weak)   FORKS=$((NODES * CORES_PER_NODE)); SAMPLES=$((BASE_SAMPLES * NODES)) ;;
    *) echo "Unknown mode: ${MODE}" >&2; exit 2 ;;
esac

# Take the first SAMPLES rows of the samplesheet, keeping the header.
head -n $((SAMPLES + 1)) "${INPUT}" > "bench_samplesheet_${LABEL}.csv"
ACTUAL_SAMPLES=$(($(wc -l < "bench_samplesheet_${LABEL}.csv") - 1))
if [[ ${ACTUAL_SAMPLES} -lt ${SAMPLES} ]]; then
    echo "[warn] samplesheet has only ${ACTUAL_SAMPLES} rows; ${MODE} scaling at" \
         "${NODES} nodes wanted ${SAMPLES}. Weak-scaling efficiency will be wrong." >&2
fi

# ---- Headers ---------------------------------------------------------------
printf 'nodes\tmode\treplicate\twall_seconds\tcpu_energy_kJ\tgpu_energy_kJ\tn_samples\tlabel\n' > timings.tsv
printf 'nodes\tmode\treplicate\tbackend\tenergy_kJ\tlabel\n' > energy.tsv

# ---- Replicate loop --------------------------------------------------------
for REP in $(seq 1 "${REPLICATES}"); do
    REP_LABEL="${LABEL}_rep${REP}"
    echo "[$(date -Is)] starting ${REP_LABEL} (${ACTUAL_SAMPLES} samples, maxForks=${FORKS})"

    NVS_PID=""
    if command -v nvidia-smi >/dev/null 2>&1; then
        (nvidia-smi --query-gpu=timestamp,index,power.draw,memory.used \
                    --format=csv,noheader,nounits \
                    -l 1 > "nvidia_smi_rep${REP}.log") &
        NVS_PID=$!
    fi

    RAPL_PID=""
    if perf list 2>/dev/null | grep -q "power/energy-pkg/"; then
        (perf stat -a -e power/energy-pkg/ -I 1000 -- sleep 100000 \
            2> "rapl_${REP_LABEL}.log") &
        RAPL_PID=$!
    fi

    START=$(date +%s)

    nextflow run "${PIPELINE_DIR}/main.nf" \
        -profile "${BENCH_PROFILE}" \
        --input       "bench_samplesheet_${LABEL}.csv" \
        --skip_sc     true \
        --skip_deconv true \
        --max_cpus    "${CORES_PER_NODE}" \
        --max_memory  256.GB \
        --outdir      "bench_${REP_LABEL}" \
        -work-dir     "work_${REP_LABEL}" \
        -ansi-log     false \
        -process.maxForks "${FORKS}" \
        > "nextflow_${REP_LABEL}.log" 2>&1

    END=$(date +%s)
    WALL=$((END - START))

    if [[ -n "${NVS_PID}"  ]]; then kill "${NVS_PID}"  2>/dev/null || true; fi
    if [[ -n "${RAPL_PID}" ]]; then kill "${RAPL_PID}" 2>/dev/null || true; fi
    wait 2>/dev/null || true

    gpu_energy_kJ=0
    if [[ -s "nvidia_smi_rep${REP}.log" ]]; then
        gpu_energy_kJ=$(awk -F, '{s += $3} END {printf "%.3f", s/1000}' "nvidia_smi_rep${REP}.log")
    fi
    cpu_energy_kJ=0
    if [[ -s "rapl_${REP_LABEL}.log" ]]; then
        cpu_energy_kJ=$(awk '/energy-pkg/ {s+=$2} END {printf "%.3f", s/1000}' "rapl_${REP_LABEL}.log")
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${NODES}" "${MODE}" "${REP}" "${WALL}" \
        "${cpu_energy_kJ}" "${gpu_energy_kJ}" "${ACTUAL_SAMPLES}" "${LABEL}" >> timings.tsv

    printf '%s\t%s\t%s\tCPU\t%s\t%s\n' "${NODES}" "${MODE}" "${REP}" "${cpu_energy_kJ}" "${LABEL}" >> energy.tsv
    printf '%s\t%s\t%s\tGPU\t%s\t%s\n' "${NODES}" "${MODE}" "${REP}" "${gpu_energy_kJ}" "${LABEL}" >> energy.tsv

    echo "[$(date -Is)] ${REP_LABEL} done: wall=${WALL}s," \
         "CPU=${cpu_energy_kJ} kJ, GPU=${gpu_energy_kJ} kJ"

    rm -rf "work_${REP_LABEL}"
done

echo "[$(date -Is)] benchmark point complete: ${LABEL}, ${REPLICATES} replicates"
