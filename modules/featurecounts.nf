// modules/featurecounts.nf
// Subread featureCounts 2.0.6 - gene-level quantification from STAR/HISAT2 BAMs.
//
// Aliased twice in workflows/bulk_rnaseq.nf (FEATURECOUNTS_STAR and
// FEATURECOUNTS_HISAT2) so both aligners are counted by identical code with
// identical settings. publishDir is per-aligner to keep the outputs separate.


process FEATURECOUNTS {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir path: { "${params.outdir}/quant/featurecounts/${meta.aligner}" }, mode: 'copy'

    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'

    input:
    tuple val(meta), path(bam), path(bai)
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.counts.tsv"), emit: counts
    path  "${meta.id}.counts.tsv.summary",          emit: summary
    path  "versions.yml",                           emit: versions

    script:
    def strand_flag = meta.strandedness == 'reverse' ? 2 : meta.strandedness == 'forward' ? 1 : 0
    def paired      = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -F GTF \\
        -t exon \\
        -g gene_id \\
        -s ${strand_flag} \\
        ${paired} \\
        -o ${meta.id}.counts.tsv \\
        ${bam}

    cat <<-VER > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oP 'v[0-9.]+')
    VER
    """
}
