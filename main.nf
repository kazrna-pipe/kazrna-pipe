#!/usr/bin/env nextflow
/*
 * KazRNA-Pipe: integrated bulk + single-cell RNA-seq pipeline.
 *
 * Entry point. Delegates to workflows/ for the per-module orchestration.
 *
 * Usage:
 *   nextflow run main.nf -profile test,singularity
 *   nextflow run main.nf -profile slurm,singularity --input samples.csv --sc_input sc_samples.csv
 *
 * See REPRODUCIBILITY.md for a step-by-step guide.
 *
 * NOTE on structure: Nextflow's current script parser does not allow bare
 * statements at the top level of a script. Parameter validation and the
 * launch banner therefore live inside functions called from the workflow
 * body, rather than executing at parse time.
 */

include { BULK_RNASEQ          } from './workflows/bulk_rnaseq.nf'
include { SCRNASEQ_CPU         } from './workflows/scrnaseq_cpu.nf'
include { SCRNASEQ_GPU         } from './workflows/scrnaseq_gpu.nf'
include { DECONVOLUTION        } from './workflows/deconvolution.nf'
include { HPC_BENCHMARK        } from './workflows/benchmark.nf'
include { CLUSTERING_AGREEMENT } from './modules/clustering_agreement.nf'
include { SOFTWARE_VERSIONS    } from './modules/software_versions.nf'

// ---------------------------------------------------------------------------
// Parameter validation
// ---------------------------------------------------------------------------

def validate_params() {
    def required = []
    if (!params.skip_bulk) required << 'input'
    if (!params.skip_sc)   required << 'sc_input'

    def missing = required.findAll { !params[it] }
    if (missing) {
        error "Missing required parameter(s): ${missing.collect { "--${it}" }.join(', ')}"
    }
}

// ---------------------------------------------------------------------------
// Launch banner - part of the run's audit trail
// ---------------------------------------------------------------------------

def print_banner() {
    log.info """
=============================================================
KazRNA-Pipe   v${workflow.manifest.version}
=============================================================
input              : ${params.input ?: '(skipped)'}
sc_input           : ${params.sc_input ?: '(skipped)'}
tcga_input         : ${params.tcga_input ?: '(skipped)'}
refs_dir           : ${params.refs_dir}
outdir             : ${params.outdir}
profile            : ${workflow.profile}
container engine   : ${workflow.containerEngine ?: 'none'}
nextflow version   : ${nextflow.version}
launchDir          : ${workflow.launchDir}
workDir            : ${workflow.workDir}
=============================================================
""".stripIndent()
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------

workflow {

    validate_params()
    print_banner()

    versions_ch = Channel.empty()

    // 1. Bulk RNA-seq
    if (!params.skip_bulk) {
        ch_bulk_samples = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, row.condition, file(row.fastq_1), file(row.fastq_2)) }

        BULK_RNASEQ(ch_bulk_samples, params.refs_dir)
        versions_ch = versions_ch.mix(BULK_RNASEQ.out.versions)
    }

    // 2. Single-cell RNA-seq (CPU + GPU paths run in parallel for benchmarking)
    if (!params.skip_sc) {
        ch_sc_samples = Channel
            .fromPath(params.sc_input)
            .splitCsv(header: true)
            .map { row -> tuple([id: row.sample_id, condition: row.condition], file(row.matrix_h5)) }

        ch_marker_yaml = file(params.marker_yaml)

        SCRNASEQ_CPU(ch_sc_samples, ch_marker_yaml)
        versions_ch = versions_ch.mix(SCRNASEQ_CPU.out.versions)

        if (params.run_gpu) {
            SCRNASEQ_GPU(ch_sc_samples, ch_marker_yaml)
            versions_ch = versions_ch.mix(SCRNASEQ_GPU.out.versions)

            // Figure 3C: ARI / NMI / ASW between the CPU and GPU clusterings.
            CLUSTERING_AGREEMENT(
                SCRNASEQ_CPU.out.clusters,
                SCRNASEQ_GPU.out.clusters
            )
            versions_ch = versions_ch.mix(CLUSTERING_AGREEMENT.out.versions)
        }
    }

    // 3. Cross-modality deconvolution (requires both bulk and sc outputs)
    if (!params.skip_bulk && !params.skip_sc && !params.skip_deconv) {
        DECONVOLUTION(
            BULK_RNASEQ.out.counts,
            SCRNASEQ_CPU.out.reference
        )
        versions_ch = versions_ch.mix(DECONVOLUTION.out.versions)
    }

    // 4. HPC benchmark (strong/weak scaling, energy)
    if (params.run_benchmark) {
        HPC_BENCHMARK()
        versions_ch = versions_ch.mix(HPC_BENCHMARK.out.versions)
    }

    // 5. Record every tool version actually invoked, alongside container digests.
    //
    // Every module emits a file literally named versions.yml, so collecting
    // them stages 20+ identically named files into one directory and Nextflow
    // refuses. collectFile() concatenates them into a single document first;
    // unique() drops the duplicates that arise when a process runs per sample.
    ch_versions_collated = versions_ch
        .unique()
        .collectFile(name: 'collated_versions.yml', newLine: true, sort: true)

    SOFTWARE_VERSIONS(ch_versions_collated)
}
