// modules/concordance.nf
// Cross-aligner and cross-method concordance.
//
// Takes one merged gene x sample matrix per aligner and one DEG table per
// DE method, as six explicit paths, so that staged filenames are deterministic
// and the script can tell which file is which.

process CONCORDANCE {
    label      'process_medium'
    publishDir "${params.outdir}/concordance", mode: 'copy'

    input:
    path star_counts
    path hisat2_counts
    path salmon_counts
    path deseq2_degs
    path edger_degs
    path limma_degs

    output:
    path "spearman_matrix.tsv",         emit: spearman
    path "spearman_matrix.pdf",         emit: spearman_plot
    path "core_degs_3way.tsv",          emit: overlap
    path "venn_DESeq2_edgeR_limma.pdf", emit: venn
    path "method_summary.tsv",          emit: summary
    path "concordance_provenance.json", emit: provenance
    path "versions.yml",                emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/concordance_analysis.R \\
        --star_counts   ${star_counts} \\
        --hisat2_counts ${hisat2_counts} \\
        --salmon_counts ${salmon_counts} \\
        --deseq2_degs   ${deseq2_degs} \\
        --edger_degs    ${edger_degs} \\
        --limma_degs    ${limma_degs} \\
        --fdr           ${params.fdr_threshold} \\
        --lfc           ${params.lfc_threshold} \\
        --outdir        .

    cat <<-VER > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | grep -oP '[0-9.]+' | head -1)
        VennDiagram: \$(Rscript -e 'cat(as.character(packageVersion("VennDiagram")))')
        pheatmap: \$(Rscript -e 'cat(as.character(packageVersion("pheatmap")))')
    VER
    """
}
