/*
 * Trim Galore v0.6.10
 *
 * Adapter and quality trimming. Defaults match the nf-core/rnaseq configuration
 * to keep cross-pipeline comparisons honest.
 */

process TRIMGALORE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/bulk/trimgalore", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.fq.gz') ? null : fn }

    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_R1_val_1.fq.gz"),
                                          path("${sample_id}_R2_val_2.fq.gz"),
        emit: reads
    path "*_trimming_report.txt", emit: report
    path 'versions.yml',          emit: versions

    script:
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --output_dir . \\
        --basename ${sample_id} \\
        ${fastq_1} ${fastq_2}

    # rename to predictable filenames
    mv ${sample_id}_val_1.fq.gz ${sample_id}_R1_val_1.fq.gz 2>/dev/null || true
    mv ${sample_id}_val_2.fq.gz ${sample_id}_R2_val_2.fq.gz 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trim_galore: \$( trim_galore --version 2>&1 | sed -n 's/^.*Trim Galore version: //p' )
        cutadapt:    \$( cutadapt --version )
    END_VERSIONS
    """
}
