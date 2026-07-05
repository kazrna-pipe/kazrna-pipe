// modules/bayesprism.nf

process BAYESPRISM {
    label      'process_high_memory'
    publishDir "${params.outdir}/deconvolution/bayesprism", mode: 'copy'

    input:
    path bulk_counts
    path sc_reference     // .rds with cell-type-labelled Seurat object
    path sample_sheet

    output:
    path "bayesprism_theta.tsv",       emit: theta
    path "bayesprism_provenance.json", emit: provenance
    path "versions.yml",               emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/bayesprism_deconv.R \\
        --bulk         ${bulk_counts} \\
        --sc_reference ${sc_reference} \\
        --samples      ${sample_sheet} \\
        --threads      ${task.cpus} \\
        --out_prefix   bayesprism

    cat <<-VER > versions.yml
    "${task.process}":
        BayesPrism: \$(Rscript -e 'cat(as.character(packageVersion("BayesPrism")))')
    VER
    """
}
