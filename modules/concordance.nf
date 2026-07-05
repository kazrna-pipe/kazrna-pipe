// modules/concordance.nf
// Wraps scripts/R/concordance_analysis.R - produces cross-aligner Spearman
// matrix and cross-method DEG overlap (Figure 2A and 2B).

process CONCORDANCE {
    label      'process_low'
    publishDir "${params.outdir}/concordance", mode: 'copy'

    input:
    path counts_files     // e.g. star.counts.tsv, hisat2.counts.tsv, salmon.gene_counts.tsv
    path de_results       // e.g. deseq2/star/deseq2_results.tsv, ...
    path sample_sheet

    output:
    path "fig2a_spearman_matrix.tsv",        emit: spearman
    path "fig2a_spearman_heatmap.pdf",       emit: spearman_plot
    path "fig2b_deg_overlap.tsv",            emit: overlap
    path "fig2b_deg_venn.pdf",               emit: venn
    path "concordance_provenance.json",      emit: provenance
    path "versions.yml",                     emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/concordance_analysis.R \\
        --counts_dir   . \\
        --de_dir       . \\
        --samples      ${sample_sheet} \\
        --out_dir      .

    cat <<-VER > versions.yml
    "${task.process}":
        VennDiagram:    \$(Rscript -e 'cat(as.character(packageVersion("VennDiagram")))')
        ComplexHeatmap: \$(Rscript -e 'cat(as.character(packageVersion("ComplexHeatmap")))')
    VER
    """
}
