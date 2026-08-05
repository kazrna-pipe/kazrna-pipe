/*
 * Cross-modality deconvolution.
 *
 * Three deconvolution methods (CIBERSORTx, BayesPrism, MuSiC) run in parallel
 * against the GSE160269 single-cell reference. BayesPrism is the primary
 * method based on the highest cross-method agreement.
 *
 * CIBERSORTx needs a registered Stanford token and a non-redistributable
 * image, so it is optional. DECONV_AGREEMENT always takes three inputs; when
 * CIBERSORTx is not run, an empty-list placeholder is passed and the module's
 * `have_cibersort` guard drops the flag from the R command line.
 */

include { BAYESPRISM       } from '../modules/bayesprism.nf'
include { CIBERSORTX       } from '../modules/cibersortx.nf'
include { MUSIC            } from '../modules/music.nf'
include { DECONV_AGREEMENT } from '../modules/deconv_agreement.nf'
include { CELLTYPE_DE      } from '../modules/celltype_de.nf'

workflow DECONVOLUTION {

    take:
    bulk_counts     // gene x sample matrix
    sc_reference    // single-cell reference object

    main:

    ch_versions = Channel.empty()

    // BayesPrism is always run.
    BAYESPRISM(bulk_counts, sc_reference)
    ch_versions = ch_versions.mix(BAYESPRISM.out.versions)

    // MuSiC is always run.
    MUSIC(bulk_counts, sc_reference)
    ch_versions = ch_versions.mix(MUSIC.out.versions)

    // CIBERSORTx only if a Stanford token is available (registration required).
    // The token itself is read from params inside the module, so it is not
    // passed as a channel argument.
    if (params.cibersortx_token) {
        CIBERSORTX(bulk_counts, sc_reference)
        ch_cibersortx = CIBERSORTX.out.fractions
        ch_versions = ch_versions.mix(CIBERSORTX.out.versions)
    } else {
        ch_cibersortx = Channel.value([])
        log.warn "No CIBERSORTx token provided (params.cibersortx_token). " +
                 "Deconvolution will proceed with BayesPrism and MuSiC only. " +
                 "Register at https://cibersortx.stanford.edu/ to enable."
    }

    // Cross-method agreement (Figure 5B) and consensus proportions
    DECONV_AGREEMENT(
        BAYESPRISM.out.proportions,
        MUSIC.out.proportions,
        ch_cibersortx
    )
    ch_versions = ch_versions.mix(DECONV_AGREEMENT.out.versions)

    // Cell-type-specific DE using BayesPrism proportions (Figure 5C, 5D, and the
    // epithelial-specific DEGs)
    CELLTYPE_DE(bulk_counts, BAYESPRISM.out.proportions)
    ch_versions = ch_versions.mix(CELLTYPE_DE.out.versions)

    emit:
    proportions   = BAYESPRISM.out.proportions
    celltype_degs = CELLTYPE_DE.out.degs
    versions      = ch_versions
}
