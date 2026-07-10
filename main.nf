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
 */

nextflow.enable.dsl = 2

include { BULK_RNASEQ      } from './workflows/bulk_rnaseq.nf'
include { SCRNASEQ_CPU     } from './workflows/scrnaseq_cpu.nf'
include { SCRNASEQ_GPU     } from './workflows/scrnaseq_gpu.nf'
include { DECONVOLUTION    } from './workflows/deconvolution.nf'
include { HPC_BENCHMARK    } from './workflows/benchmark.nf'
include { SOFTWARE_VERSIONS } from './modules/software_versions.nf'

// ---------------------------------------------------------------------------
// Parameter validation
// ---------------------------------------------------------------------------

def required_params = []
if (!params.skip_bulk)   required_params << 'input'
if (!params.skip_sc)     required_params << 'sc_input'

required_params.each { p ->
    if (!params[p]) {
        log.error "Missing required parameter --${p}"
        System.exit(1)
    }
}

// Pretty-print the launch configuration for the audit trail.
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

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------

workflow {

    versions_ch = Channel.empty()

    // 1. Bulk RNA-seq
    if (!params.skip_bulk) {
        ch_bulk_samples = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, row.condition, file(row.fastq_1), file(row.fastq_2)) }

        BULK_RNASEQ(ch_bulk_samples, params.refs_dir)
        bulk_counts = BULK_RNASEQ.out.counts
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
    SOFTWARE_VERSIONS(versions_ch.collect())
}

// ---------------------------------------------------------------------------
// Completion handler
// ---------------------------------------------------------------------------

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    log.info """
    =============================================================
    KazRNA-Pipe completed with status: ${status}
    Duration         : ${workflow.duration}
    Results directory: ${params.outdir}
    Command          : ${workflow.commandLine}
    =============================================================
    """.stripIndent()
}
