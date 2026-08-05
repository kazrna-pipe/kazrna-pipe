// modules/edger.nf

process EDGER {
    tag        "${meta.aligner}"
    label      'process_medium'
    publishDir path: { "${params.outdir}/de/edger/${meta.aligner}" }, mode: 'copy'

    input:
    tuple val(meta), path(counts)
    path  sample_sheet

    output:
    tuple val(meta), path("edger_results.tsv"),     emit: results
    tuple val(meta), path("edger_significant.tsv"), emit: significant
    path  "edger_provenance.json",                   emit: provenance
    path  "edger_session.txt",                       emit: session
    path  "versions.yml",                            emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/edger_analysis.R \\
        --counts ${counts} \\
        --meta   ${sample_sheet} \\
        --fdr    ${params.fdr_threshold} \\
        --lfc    ${params.lfc_threshold} \\
        --outdir .

    cat <<-VER > versions.yml
    "${task.process}":
        edgeR: \$(Rscript -e 'cat(as.character(packageVersion("edgeR")))')
    VER
    """
}
