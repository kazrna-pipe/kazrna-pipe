// workflows/benchmark.nf
// Strong- and weak-scaling benchmarks on the KazHPC SLURM cluster, including
// RAPL CPU energy counters and nvidia-smi GPU power sampling. Produces the
// data backing Figure 4 and Table 2 in the manuscript.
//
// Run with: nextflow run main.nf -profile slurm --benchmark_only true
//                                --benchmark_nodes 1,2,4,8,16

include { BENCHMARK_RUN } from '../modules/benchmark_run.nf'

workflow BENCHMARK {

    take:
    benchmark_input        // path to a small representative dataset

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
        benchmark_input
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
}
