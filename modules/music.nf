// modules/music.nf

process MUSIC {
    label      'process_medium'
    publishDir "${params.outdir}/deconvolution/music", mode: 'copy'

    input:
    path bulk_counts
    path sc_reference

    output:
    path "music_proportions.tsv",    emit: proportions
    path "music_provenance.json",    emit: provenance
    path "versions.yml",             emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/music_deconv.R \\
        --bulk         ${bulk_counts} \\
        --sc_reference ${sc_reference} \\
        --threads      ${task.cpus} \\
        --out_prefix   music

    cat <<-VER > versions.yml
    "${task.process}":
        MuSiC: \$(Rscript -e 'cat(as.character(packageVersion("MuSiC")))')
    VER
    """
}
