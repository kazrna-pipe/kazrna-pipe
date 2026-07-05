// modules/benchmark_run.nf
// Runs a single benchmark point: N nodes x mode (strong | weak), records
// wall time, CPU/GPU utilisation, and energy (RAPL + nvidia-smi).

process BENCHMARK_RUN {
    tag        "${nodes}node_${mode}"
    label      'process_benchmark'
    publishDir "${params.outdir}/benchmark/${mode}/${nodes}n", mode: 'copy'

    input:
    tuple val(nodes), val(mode)
    path  bench_input

    output:
    path "timings.tsv",   emit: timings
    path "energy.tsv",    emit: energy
    path "nvidia_smi.log", emit: gpu_log, optional: true
    path "rapl_*.log",     emit: cpu_log, optional: true

    script:
    """
    bash ${projectDir}/scripts/benchmark_runner.sh \\
        --nodes      ${nodes} \\
        --mode       ${mode} \\
        --input      ${bench_input} \\
        --out_dir    . \\
        --label      "${nodes}n_${mode}"
    """
}
