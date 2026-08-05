// modules/samtools_sort.nf
// Coordinate-sort a SAM/BAM and build its index.
//
// Split out of HISAT2_ALIGN so every process runs in a single-tool container.
// Emits the same (meta, bam, bai) shape that FEATURECOUNTS expects, so both
// aligner paths converge on one counting process.

process SAMTOOLS_SORT {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir path: { "${params.outdir}/bulk/${meta.aligner}" }, mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "versions.yml",                                                 emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${alignment}
    samtools index -@ ${task.cpus} ${meta.id}.bam

    cat <<-VER > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -1 | awk '{print \$2}' )
    VER
    """
}
