// modules/bayesprism.nf
//


process BAYESPRISM {
    label      'process_high_memory'
    publishDir "${params.outdir}/deconvolution/bayesprism", mode: 'copy'

    input:
    path bulk_counts
    path sc_reference     // .rds with cell-type-labelled Seurat object

    output:
    path "bayesprism_proportions.tsv", emit: proportions
    path "bayesprism_cellexp_Z.rds",   emit: cellexp
    path "provenance.json",            emit: provenance
    path "versions.yml",               emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/bayesprism_deconv.R \\
        --bulk_counts  ${bulk_counts} \\
        --sc_reference ${sc_reference} \\
        --n_cores      ${task.cpus} \\
        --n_iter       ${params.bayesprism_iterations} \\
        --n_signatures ${params.bayesprism_n_signatures} \\
        --outdir       .

    cat <<-VER > versions.yml
    "${task.process}":
        BayesPrism: \$(Rscript -e 'cat(as.character(packageVersion("BayesPrism")))')
    VER
    """
}
