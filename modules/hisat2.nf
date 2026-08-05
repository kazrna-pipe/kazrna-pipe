/*
 * HISAT2 v2.2.1 - alignment only.
 *
 * Graph-based alignment as an independent cross-check on STAR. The high
 * Spearman rho between STAR and HISAT2 in the manuscript (0.974) is partly an
 * artefact of their shared design philosophy; both are spliced alignment-based
 * methods. The Salmon comparison is the more informative one.
 *
 * This process runs only hisat2. Sorting is done by SAMTOOLS_SORT and counting
 * by the shared FEATURECOUNTS process, so that:
 *   - each process needs only a single-tool container, which is what
 *     BioContainers publishes;
 *   - HISAT2 and STAR both flow through the *same* counting process, so their
 *     concordance reflects the aligners rather than differing count settings.
 *
 * The index basename is derived from the *.1.ht2 file present in the index
 * directory rather than hard-coded, so any index naming convention works.
 */

process HISAT2_ALIGN {
    tag        "${sample_id}"
    label      'process_high'
    publishDir "${params.outdir}/bulk/hisat2", mode: params.publish_dir_mode,
               pattern: "*.summary.txt"

    container 'quay.io/biocontainers/hisat2:2.2.1--h87f3376_5'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  hisat2_index

    output:
    // The SAM is not published: with -p > 1 its record order depends on thread
    // completion order, so two runs differ byte-for-byte while containing the
    // same alignments. SAMTOOLS_SORT publishes the coordinate-sorted BAM,
    // which is deterministic and is the artefact worth keeping.
    tuple val(sample_id), val(condition), path("${sample_id}.sam"), emit: sam
    path "${sample_id}.hisat2.summary.txt",                         emit: summary
    path 'versions.yml',                                            emit: versions

    script:
    """
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
        -S ${sample_id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -1 | awk '{print \$3}' )
    END_VERSIONS
    """
}
