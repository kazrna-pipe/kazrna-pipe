// modules/tximport.nf
// tximport 1.32.0 - aggregate Salmon transcript-level estimates to gene level.

process TXIMPORT {
    label      'process_low'
    publishDir "${params.outdir}/quant/tximport", mode: 'copy'

    input:
    path quant_dirs   // collected salmon/<sample>/quant.sf trees
    path tx2gene
    path sample_sheet

    output:
    path "salmon.gene_counts.tsv",  emit: counts
    path "salmon.gene_tpm.tsv",     emit: tpm
    path "salmon.gene_lengths.tsv", emit: lengths
    path "versions.yml",            emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/tximport_aggregate.R \\
        --quant_root  . \\
        --tx2gene     ${tx2gene} \\
        --samples     ${sample_sheet} \\
        --out_prefix  salmon.gene

    cat <<-VER > versions.yml
    "${task.process}":
        tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
        R: \$(R --version | head -1 | grep -oP '[0-9.]+' | head -1)
    VER
    """
}
