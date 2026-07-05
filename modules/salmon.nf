/*
 * Salmon v1.10.2
 *
 * Selective alignment to a decoy-aware GENCODE v44 transcriptome.
 * Salmon's lower agreement with STAR/HISAT2 in the manuscript (ρ ≈ 0.95)
 * is the expected, well-documented behaviour and motivates its inclusion.
 */

process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/bulk/salmon", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/salmon:1.10.2--h7e5ed60_1'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  salmon_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}/quant.sf"), emit: quant
    path "${sample_id}/cmd_info.json",  emit: cmd_info
    path "${sample_id}/lib_format_counts.json", emit: lib
    path 'versions.yml',                emit: versions

    script:
    """
    salmon quant \\
        --index ${salmon_index} \\
        --libType A \\
        --threads ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        --numBootstraps 100 \\
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        -o ${sample_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/^salmon //' )
    END_VERSIONS
    """
}
