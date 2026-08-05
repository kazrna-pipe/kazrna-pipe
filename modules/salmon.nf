/*
 * Salmon v1.10.2
 *
 * Selective alignment to a decoy-aware GENCODE v44 transcriptome.
 * Salmon's lower agreement with STAR/HISAT2 in the manuscript (rho ~ 0.95)
 * is the expected, well-documented behaviour and motivates its inclusion.
 *
 * NOTE: --numBootstraps is deliberately not set. Bootstrap replicates are the
 * only stochastic component of salmon quant and salmon has no --seed option in
 * mapping-based mode, so enabling them makes quant.sf differ between identical
 * runs. Nothing downstream consumes the bootstraps: tximport_aggregate.R reads
 * only quant.sf. Add --numBootstraps back if inferential replicates are ever
 * needed, and exclude the salmon bootstrap directories from the determinism
 * check.
 *
 * NOTE: the whole per-sample directory is emitted, not <sample>/quant.sf.
 * Collecting bare quant.sf files would stage several identically named files
 * into one directory, which Nextflow silently renames (quant.sf, quant.2.sf,
 * ...), destroying the sample-to-file mapping tximport depends on.
 */

process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/bulk/salmon", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/salmon:1.10.2--hecfa306_0'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  salmon_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}"), emit: quant
    path 'versions.yml',                                        emit: versions

    script:
    """
    salmon quant \\
        --index ${salmon_index} \\
        --libType A \\
        --threads ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        -o ${sample_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/^salmon //' )
    END_VERSIONS
    """
}
