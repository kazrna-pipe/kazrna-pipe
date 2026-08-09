/*
 * HISAT2 v2.2.1 - alignment only.
 */

process HISAT2_ALIGN {
    tag        "${sample_id}"
    label      'process_high'
    publishDir "${params.outdir}/bulk/hisat2", mode: params.publish_dir_mode,
               pattern: "*.summary.txt"

    container 'community.wave.seqera.io/library/hisat2_samtools@sha256:c34a62b3e0c61c75a10f869eaf75061a464d2177c32bb697ac864db4039a2ff4'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  hisat2_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    path "${sample_id}.hisat2.summary.txt",                         emit: summary
    path 'versions.yml',                                            emit: versions

    script:
    """
    set -o pipefail

    INDEX_BASE=\$(ls ${hisat2_index}/*.1.ht2 | head -1 | sed 's/\\.1\\.ht2\$//')
    echo "Using HISAT2 index basename: \${INDEX_BASE}"

    hisat2 \\
        -x "\${INDEX_BASE}" \\
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        -p ${task.cpus} \\
        --summary-file ${sample_id}.hisat2.summary.txt \\
        --new-summary \\
        --seed ${params.seed} \\
        | samtools sort -@ 4 -m 2G -o ${sample_id}.bam -

    samtools index -@ 4 ${sample_id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -1 | awk '{print \$3}' )
        samtools: \$( samtools --version | head -1 | awk '{print \$2}' )
    END_VERSIONS
    """
}
