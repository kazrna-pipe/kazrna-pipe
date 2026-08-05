// modules/merge_counts.nf
// Merge per-sample featureCounts output into a single gene x sample matrix.
//
// Required by both the DE modules and CONCORDANCE: featureCounts emits one
// file per sample, but DESeq2/edgeR/limma-voom and concordance_analysis.R all
// expect a single matrix whose column names are sample IDs.

process MERGE_COUNTS {
    tag        "${aligner}"
    label      'process_single'
    publishDir "${params.outdir}/quant/merged", mode: 'copy'

    input:
    val   aligner 
    path  count_files
    val   sample_ids

    output:
    tuple val(aligner), path("${aligner}.gene_counts.tsv"), emit: counts
    path  "versions.yml",                                   emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/merge_counts.R \\
        --counts     ${count_files.join(',')} \\
        --sample_ids ${sample_ids} \\
        --out        ${aligner}.gene_counts.tsv

    cat <<-VER > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | grep -oP '[0-9.]+' | head -1)
    VER
    """
}
