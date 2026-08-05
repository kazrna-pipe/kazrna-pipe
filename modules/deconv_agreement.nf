// modules/deconv_agreement.nf

process DECONV_AGREEMENT {
    label      'process_low'
    publishDir "${params.outdir}/deconvolution/agreement", mode: 'copy'

    input:
    path bayesprism_theta
    path music_proportions
    path cibersortx_fractions

    output:
    path "deconv_correlation.tsv", emit: corr
    path "fig5b_deconv_corr.pdf",  emit: plot
    path "deconv_agreement_provenance.json", emit: provenance
    path "versions.yml",           emit: versions

    script:
    def have_cibersort = cibersortx_fractions ? "--cibersortx ${cibersortx_fractions}" : ""
    """
    Rscript ${projectDir}/scripts/R/deconv_agreement.R \\
        --bayesprism ${bayesprism_theta} \\
        --music      ${music_proportions} \\
        ${have_cibersort} \\
        --out_dir    .

    cat <<-VER > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | grep -oP '[0-9.]+' | head -1)
    VER
    """
}
