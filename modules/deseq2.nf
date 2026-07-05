// modules/deseq2.nf
// Wraps scripts/R/deseq2_analysis.R for use inside the bulk workflow.

process DESEQ2 {
    tag        "${meta.aligner}"
    label      'process_medium'
    publishDir "${params.outdir}/de/deseq2/${meta.aligner}", mode: 'copy'

    input:
    tuple val(meta), path(counts)
    path  sample_sheet

    output:
    tuple val(meta), path("deseq2_results.tsv"),         emit: results
    tuple val(meta), path("deseq2_normalised_vst.tsv"),  emit: vst
    path  "deseq2_provenance.json",                       emit: provenance
    path  "deseq2_session_info.txt",                      emit: session
    path  "versions.yml",                                 emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/deseq2_analysis.R \\
        --counts        ${counts} \\
        --metadata      ${sample_sheet} \\
        --condition_col condition \\
        --contrast      tumor,normal \\
        --covariates    age,sex \\
        --fdr           ${params.fdr_threshold} \\
        --lfc           ${params.lfc_threshold} \\
        --out_prefix    deseq2

    cat <<-VER > versions.yml
    "${task.process}":
        DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    VER
    """
}
