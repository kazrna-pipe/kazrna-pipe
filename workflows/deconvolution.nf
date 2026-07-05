/*
 * Cross-modality deconvolution.
 *
 * Three deconvolution methods (CIBERSORTx, BayesPrism, MuSiC) run in parallel
 * against the GSE160269 single-cell reference. BayesPrism is the primary
 * method based on the highest cross-method agreement reported in the manuscript.
 *
 * Outputs Figure 4 panels and the 842 epithelial-specific DEGs.
 */

include { BAYESPRISM      } from '../modules/bayesprism.nf'
include { CIBERSORTX      } from '../modules/cibersortx.nf'
include { MUSIC           } from '../modules/music.nf'
include { DECONV_AGREEMENT } from '../modules/deconv_agreement.nf'
include { CELLTYPE_DE     } from '../modules/celltype_de.nf'

workflow DECONVOLUTION {

    take:
    bulk_counts     // gene x sample matrix
    sc_reference    // single-cell reference object (pseudobulk-ready)

    main:

    ch_versions = Channel.empty()
    ch_proportions = Channel.empty()

    // BayesPrism is always run.
    BAYESPRISM(bulk_counts, sc_reference)
    ch_proportions = ch_proportions.mix(BAYESPRISM.out.proportions.map { tuple('bayesprism', it) })
    ch_versions = ch_versions.mix(BAYESPRISM.out.versions)

    // MuSiC is always run.
    MUSIC(bulk_counts, sc_reference)
    ch_proportions = ch_proportions.mix(MUSIC.out.proportions.map { tuple('music', it) })
    ch_versions = ch_versions.mix(MUSIC.out.versions)

    // CIBERSORTx only if a Stanford token is available (registration required).
    if (params.cibersortx_token) {
        CIBERSORTX(bulk_counts, sc_reference, params.cibersortx_token)
        ch_proportions = ch_proportions.mix(CIBERSORTX.out.proportions.map { tuple('cibersortx', it) })
        ch_versions = ch_versions.mix(CIBERSORTX.out.versions)
    } else {
        log.warn "No CIBERSORTx token provided (params.cibersortx_token). " +
                 "Deconvolution will proceed with BayesPrism and MuSiC only. " +
                 "Register at https://cibersortx.stanford.edu/ to enable."
    }

    // Cross-method agreement (Figure 4A) and consensus proportions
    DECONV_AGREEMENT(ch_proportions.collect())
    ch_versions = ch_versions.mix(DECONV_AGREEMENT.out.versions)

    // Cell-type-specific DE using BayesPrism proportions (Figure 4C, 4D, and the
    // 842 epithelial-specific DEGs reported in the manuscript)
    CELLTYPE_DE(bulk_counts, BAYESPRISM.out.proportions)
    ch_versions = ch_versions.mix(CELLTYPE_DE.out.versions)

    emit:
    proportions = BAYESPRISM.out.proportions
    celltype_degs = CELLTYPE_DE.out.degs
    versions    = ch_versions
}
