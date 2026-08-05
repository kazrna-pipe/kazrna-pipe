// modules/limma_voom.nf

process LIMMA_VOOM {
    tag        "${meta.aligner}"
    label      'process_medium'
    publishDir path: { "${params.outdir}/de/limma_voom/${meta.aligner}" }, mode: 'copy'

    input:
    tuple val(meta), path(counts)
    path  sample_sheet

    output:
    tuple val(meta), path("limma_voom_results.tsv"),     emit: results
    tuple val(meta), path("limma_voom_significant.tsv"), emit: significant
    path  "limma_voom_provenance.json",                   emit: provenance
    path  "limma_voom_session.txt",                       emit: session
    path  "versions.yml",                                 emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/limma_voom_analysis.R \\
        --counts ${counts} \\
        --meta   ${sample_sheet} \\
        --fdr    ${params.fdr_threshold} \\
        --lfc    ${params.lfc_threshold} \\
        --outdir .

    cat <<-VER > versions.yml
    "${task.process}":
        limma: \$(Rscript -e 'cat(as.character(packageVersion("limma")))')
        edgeR: \$(Rscript -e 'cat(as.character(packageVersion("edgeR")))')
    VER
    """
}
