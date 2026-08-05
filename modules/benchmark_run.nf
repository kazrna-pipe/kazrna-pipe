// modules/benchmark_run.nf
// Runs a single benchmark point: N nodes x mode (strong | weak), records
// wall time, CPU/GPU utilisation, and energy (RAPL + nvidia-smi).
//
// publishDir uses closure form because its path depends on input values.
// A plain string is resolved at process-definition time, when `mode` and
// `nodes` do not yet exist; a closure is evaluated per task, when they do.

process BENCHMARK_RUN {
    tag        "${nodes}node_${mode}"
    label      'process_benchmark'
    publishDir path: { "${params.outdir}/benchmark/${mode}/${nodes}n" }, mode: 'copy'

    input:
    tuple val(nodes), val(mode)
    path  bench_input

    output:
    path "timings.tsv",    emit: timings
    path "energy.tsv",     emit: energy
    path "nvidia_smi.log", emit: gpu_log, optional: true
    path "rapl_*.log",     emit: cpu_log, optional: true
    path "versions.yml",   emit: versions

    script:
    """
    bash ${projectDir}/scripts/benchmark_runner.sh \\
        --nodes      ${nodes} \\
        --mode       ${mode} \\
        --input      ${bench_input} \\
        --replicates ${params.benchmark_replicates} \\
        --out_dir    . \\
        --label      "${nodes}n_${mode}"

    cat <<-VER > versions.yml
    "${task.process}":
        nextflow: \$(nextflow -version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "\${NXF_VER:-unknown}")
        bash:     \$(bash --version | head -1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    VER
    """
}
