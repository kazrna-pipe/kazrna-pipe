// modules/celltype_de.nf


process CELLTYPE_DE {
    label      'process_high_memory'
    publishDir "${params.outdir}/celltype_de", mode: 'copy'

    input:
    path sc_object      
    path sample_sheet

    output:
    path "celltype_de_results/*.tsv",      emit: results
    path "fig5d_epithelial_volcano.pdf",   emit: volcano
    path "celltype_de_provenance.json",    emit: provenance
    path "versions.yml",                   emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/celltype_pseudobulk_de.R \\
        --sc_object ${sc_object} \\
        --samples   ${sample_sheet} \\
        --condition tumor,normal \\
        --fdr       ${params.fdr_threshold} \\
        --lfc       ${params.lfc_threshold} \\
        --threads   ${task.cpus} \\
        --out_dir   .

    cat <<-VER > versions.yml
    "${task.process}":
        DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        Seurat: \$(Rscript -e 'cat(as.character(packageVersion("Seurat")))')
    VER
    """
}
