# Benchmarking methodology

This document specifies the scaling experiments precisely enough to be
reproduced. It exists because "the benchmark was repeated three times" is not
a reproducible statement without the workload, the hardware, and the
definitions of the reported quantities.

Every number below marked ⟦…⟧ must be filled in from your actual runs before
submission.

## Workload

| | Strong scaling | Weak scaling |
|---|---|---|
| Problem size | Fixed at `--base_samples` (default 8) | `base_samples × N` |
| Parallelism | `N × cores_per_node` | `N × cores_per_node` |
| Per-node load | Falls as N rises | Constant |
| Question answered | Does more hardware finish the same job faster? | Does a proportionally larger job take the same time? |

The base problem at N = 1 is ⟦N⟧ paired-end samples, ⟦X⟧ GB of compressed
FASTQ, ⟦Y⟧ M read pairs total. Weak scaling draws additional samples from
`data/benchmark_samplesheet.csv`, which must contain at least
`base_samples × max(nodes)` rows — otherwise the runner emits a warning and the
weak-scaling efficiency is invalid.

Only the bulk path is benchmarked (`--skip_sc --skip_deconv`), so that the
measurement is of the alignment and quantification workload rather than a mix
of unrelated stages.

## Hardware

| | |
|---|---|
| CPU | Dual Intel Xeon Gold, 16 physical cores per node |
| Memory | 128 GB per node |
| Interconnect |  |
| Shared storage |  |
| GPU (where used) | NVIDIA A100 80 GB |
| Scheduler |  |
| Nextflow | , executor `slurm` |

## Held constant across all points

- Reference indices, built once and reused
- Container images, pinned by digest (see `containers/MANIFEST.md`)
- Random seed (42)
- `--max_memory`, `--max_cpus` per node
- Nextflow work directory on ⟦shared storage / node-local disk⟧

## What is measured

Wall time is measured around the entire `nextflow run` invocation, so it
**includes** scheduler queueing within the allocation, container startup, and
all I/O to shared storage. It excludes time spent waiting for the allocation
itself. Reporting the inclusive figure is deliberate: it is what a user
experiences, and excluding I/O would flatter the scaling curve.

Energy is sampled at 1 Hz: CPU package power via `perf stat -e
power/energy-pkg/` where RAPL is readable, GPU power via `nvidia-smi
--query-gpu=power.draw`. Both are integrated over the run and reported in kJ.
Where RAPL is not readable, the CPU column is 0 and must be reported as
unavailable rather than as zero energy.

## Replication and dispersion

Each `(nodes, mode)` point is run `--replicates` times (default 3). The runner
writes one row per replicate:

```
nodes  mode  replicate  wall_seconds  cpu_energy_kJ  gpu_energy_kJ  n_samples  label
```

Figures show the **mean across replicates with error bars at ±1 standard
deviation**, computed from these raw rows. The per-replicate table is published
as `results/benchmark/fig4_scaling.tsv` and deposited with the software release,
so the dispersion can be recomputed independently.

## Definitions

For strong scaling, with `T(N)` the mean wall time at N nodes:

- Speed-up: `S(N) = T(1) / T(N)`
- Parallel efficiency: `E(N) = S(N) / N`

For weak scaling, where the problem grows with N:

- Weak-scaling efficiency: `E_w(N) = T(1) / T(N)`

`E(N) = 1` is perfect scaling; values above 1 indicate a measurement artefact
(usually caching between replicates) and should be investigated, not reported.

## Running

```bash
nextflow run main.nf -profile slurm,singularity \
    --run_benchmark true \
    --benchmark_nodes '1,2,4,8' \
    --benchmark_modes 'strong,weak' \
    --benchmark_replicates 3 \
    --benchmark_input data/benchmark_samplesheet.csv
```

Or a single point directly:

```bash
bash scripts/benchmark_runner.sh \
    --nodes 4 --mode strong \
    --input data/benchmark_samplesheet.csv \
    --replicates 3 --base_samples 8 --cores_per_node 16 \
    --out_dir results/benchmark/strong/4n --label 4n_strong
```

## Known limitations

- Node counts are simulated through `process.maxForks` rather than through
  distinct physical allocations, so inter-node communication cost is
  represented only insofar as it appears in shared-storage I/O.
- Energy figures depend on RAPL and NVML availability and are not comparable
  across sites.
- The scaling ceiling is set by the number of samples: beyond ⟦N⟧ nodes, the
  per-sample tasks no longer saturate the available cores.
