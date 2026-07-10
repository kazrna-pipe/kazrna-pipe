// workflows/benchmark.nf

include { BENCHMARK_RUN } from '../modules/benchmark_run.nf'

workflow HPC_BENCHMARK {

    main:
    Channel
        .from( params.benchmark_nodes.toString().split(',') )
        .map { it.toInteger() }
        .set { node_counts }

    Channel
        .from( params.benchmark_modes.toString().split(',') )    // strong,weak
        .set { modes }

    benchmark_jobs = node_counts.combine(modes)

    BENCHMARK_RUN(
        benchmark_jobs,
        file(params.benchmark_input)
    )

    BENCHMARK_RUN.out.timings
        .collectFile(
            name:       'fig4_scaling.tsv',
            keepHeader: true,
            skip:       1,
            storeDir:   "${params.outdir}/benchmark"
        )

    BENCHMARK_RUN.out.energy
        .collectFile(
            name:       'fig4d_energy.tsv',
            keepHeader: true,
            skip:       1,
            storeDir:   "${params.outdir}/benchmark"
        )

    emit:
    versions = BENCHMARK_RUN.out.versions
}