/*
 * Single-cell RNA-seq workflow - GPU path.
 *
 * Implements rapids-singlecell. Mirrors the CPU workflow stage-for-stage so
 * that ARI / NMI / ASW between the two paths is a like-for-like comparison.
 * Reports per-step wall-clock for manuscript Figure 3A.
 */

include { RAPIDS_LOAD       } from '../modules/rapids_singlecell.nf'
include { RAPIDS_QC         } from '../modules/rapids_singlecell.nf'
include { RAPIDS_NORMALIZE  } from '../modules/rapids_singlecell.nf'
include { RAPIDS_HVG        } from '../modules/rapids_singlecell.nf'
include { RAPIDS_PCA        } from '../modules/rapids_singlecell.nf'
include { RAPIDS_HARMONY    } from '../modules/rapids_singlecell.nf'
include { RAPIDS_NEIGHBORS  } from '../modules/rapids_singlecell.nf'
include { RAPIDS_LEIDEN     } from '../modules/rapids_singlecell.nf'
include { RAPIDS_UMAP       } from '../modules/rapids_singlecell.nf'
include { RAPIDS_MARKERS    } from '../modules/rapids_singlecell.nf'
include { CLUSTERING_AGREEMENT } from '../modules/clustering_agreement.nf'

workflow SCRNASEQ_GPU {

    take:
    ch_samples         // (sample_id, matrix_dir)

    main:

    ch_versions = Channel.empty()

    RAPIDS_LOAD(ch_samples)
    RAPIDS_QC(RAPIDS_LOAD.out.h5ad)
    RAPIDS_NORMALIZE(RAPIDS_QC.out.h5ad)
    RAPIDS_HVG(RAPIDS_NORMALIZE.out.h5ad)
    RAPIDS_PCA(RAPIDS_HVG.out.h5ad)
    RAPIDS_HARMONY(RAPIDS_PCA.out.h5ad)
    RAPIDS_NEIGHBORS(RAPIDS_HARMONY.out.h5ad)
    RAPIDS_LEIDEN(RAPIDS_NEIGHBORS.out.h5ad, params.sc_leiden_resolutions)
    RAPIDS_UMAP(RAPIDS_LEIDEN.out.h5ad)
    RAPIDS_MARKERS(RAPIDS_UMAP.out.h5ad)

    // Compare GPU labels to CPU labels at every Leiden resolution.
    // This produces Figure 3C and the JSON used by reproducibility checks.
    // The CPU labels are picked up from the published cache; if missing,
    // the comparison is skipped.
    if (file("${params.outdir}/sc/cpu/labels.tsv").exists()) {
        CLUSTERING_AGREEMENT(
            RAPIDS_LEIDEN.out.labels,
            file("${params.outdir}/sc/cpu/labels.tsv")
        )
        ch_versions = ch_versions.mix(CLUSTERING_AGREEMENT.out.versions)
    }

    ch_versions = ch_versions
        .mix(RAPIDS_LOAD.out.versions)
        .mix(RAPIDS_HARMONY.out.versions)
        .mix(RAPIDS_LEIDEN.out.versions)

    emit:
    h5ad     = RAPIDS_MARKERS.out.h5ad
    timing   = RAPIDS_MARKERS.out.timing
    versions = ch_versions
}
