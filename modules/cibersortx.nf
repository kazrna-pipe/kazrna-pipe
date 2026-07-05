// modules/cibersortx.nf
// CIBERSORTx requires a registered token and the official Docker image, which
// cannot be redistributed. This process is skipped unless params.cibersortx_token
// is set. The published Singularity image is built from the upstream Docker
// image at build time on the user's machine (see containers/Singularity.cibersortx).

process CIBERSORTX {
    label      'process_medium'
    publishDir "${params.outdir}/deconvolution/cibersortx", mode: 'copy'
    when:      params.cibersortx_token

    input:
    path bulk_mixture       // TSV: genes x samples (TPM)
    path signature_matrix   // produced from sc reference

    output:
    path "CIBERSORTx_Adjusted.txt", emit: fractions
    path "cibersortx_provenance.json", emit: provenance
    path "versions.yml",               emit: versions

    script:
    """
    /src/CIBERSORTxFractions \\
        --username  ${params.cibersortx_email} \\
        --token     ${params.cibersortx_token} \\
        --mixture   ${bulk_mixture} \\
        --sigmatrix ${signature_matrix} \\
        --rmbatchSmode TRUE \\
        --outdir .

    python ${projectDir}/scripts/python/write_provenance.py \\
        --tool     cibersortx \\
        --image    \${SINGULARITY_CONTAINER:-cibersortx-local} \\
        --out      cibersortx_provenance.json

    cat <<-VER > versions.yml
    "${task.process}":
        cibersortx_image: \${CIBERSORTX_IMAGE_DIGEST:-unknown}
    VER
    """
}
