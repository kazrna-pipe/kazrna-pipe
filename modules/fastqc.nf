/*
 * FastQC v0.12.1
 *
 * Per-sample read quality reports. Outputs are aggregated by MultiQC downstream.
 */

process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/bulk/fastqc", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip"), emit: reports
    path 'versions.yml', emit: versions

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --quiet \\
        ${fastq_1} ${fastq_2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version 2>&1 | sed 's/^FastQC v//' )
    END_VERSIONS
    """
}
