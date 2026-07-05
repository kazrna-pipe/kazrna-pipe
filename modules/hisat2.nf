/*
 * HISAT2 v2.2.1
 *
 * Graph-based alignment as an independent cross-check on STAR. The high
 * Spearman ρ between STAR and HISAT2 in the manuscript (0.974) is partly an
 * artefact of their shared design philosophy; both are spliced alignment-based
 * methods. The Salmon comparison is the more informative one.
 */

process HISAT2_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/bulk/hisat2", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/hisat2:2.2.1--h87f3376_5'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  hisat2_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam"), emit: bam
    path "${sample_id}.hisat2.summary.txt", emit: summary
    path "${sample_id}.featurecounts.txt",  emit: counts
    path 'versions.yml',                    emit: versions

    script:
    def index_base = "${hisat2_index}/grch38"
    """
    /usr/bin/time -f '%e\\t%M' -o ${sample_id}.timing.raw \\
    bash -c "
        hisat2 \\
            -x ${index_base} \\
            -1 ${fastq_1} \\
            -2 ${fastq_2} \\
            -p ${task.cpus} \\
            --summary-file ${sample_id}.hisat2.summary.txt \\
            --new-summary | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    "

    samtools index ${sample_id}.bam

    # gene-level counts via featureCounts for parity with STAR
    featureCounts \\
        -T ${task.cpus} \\
        -p --countReadPairs \\
        -a ${hisat2_index}/../gencode.v44.annotation.gtf \\
        -o ${sample_id}.featurecounts.txt \\
        ${sample_id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2:       \$( hisat2 --version | head -1 | awk '{print \$3}' )
        samtools:     \$( samtools --version | head -1 | awk '{print \$2}' )
        featureCounts: \$( featureCounts -v 2>&1 | head -2 | tail -1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
