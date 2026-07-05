// modules/clusterprofiler.nf
// clusterProfiler 4.12.6 - KEGG + GO over-representation on DE gene sets.

process CLUSTERPROFILER {
    tag        "${meta.aligner}_${meta.method}"
    label      'process_low'
    publishDir "${params.outdir}/enrichment/${meta.method}/${meta.aligner}", mode: 'copy'

    input:
    tuple val(meta), path(de_results)

    output:
    path "kegg_enrichment.tsv", emit: kegg, optional: true
    path "go_bp_enrichment.tsv", emit: go,   optional: true
    path "enrichment_provenance.json", emit: provenance
    path "versions.yml",        emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/enrichment_analysis.R \\
        --de_results    ${de_results} \\
        --organism      hsa \\
        --fdr           ${params.fdr_threshold} \\
        --lfc           ${params.lfc_threshold} \\
        --out_dir       .

    cat <<-VER > versions.yml
    "${task.process}":
        clusterProfiler: \$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))')
        org.Hs.eg.db:    \$(Rscript -e 'cat(as.character(packageVersion("org.Hs.eg.db")))')
    VER
    """
}
