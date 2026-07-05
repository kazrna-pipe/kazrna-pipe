// modules/edger.nf

process EDGER {
    tag        "${meta.aligner}"
    label      'process_medium'
    publishDir "${params.outdir}/de/edger/${meta.aligner}", mode: 'copy'

    input:
    tuple val(meta), path(counts)
    path  sample_sheet

    output:
    tuple val(meta), path("edger_results.tsv"), emit: results
    path  "edger_provenance.json",              emit: provenance
    path  "edger_session_info.txt",             emit: session
    path  "versions.yml",                       emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/edger_analysis.R \\
        --counts        ${counts} \\
        --metadata      ${sample_sheet} \\
        --condition_col condition \\
        --contrast      tumor,normal \\
        --covariates    age,sex \\
        --fdr           ${params.fdr_threshold} \\
        --lfc           ${params.lfc_threshold} \\
        --out_prefix    edger

    cat <<-VER > versions.yml
    "${task.process}":
        edgeR: \$(Rscript -e 'cat(as.character(packageVersion("edgeR")))')
    VER
    """
}
