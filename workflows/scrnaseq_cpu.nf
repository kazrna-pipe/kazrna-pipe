/*
 * Single-cell RNA-seq workflow - CPU path.
 *
 * Implements the Seurat v5 reference. Identical input matrices, normalisation
 * targets and HVG selection are shared with the GPU path so that clustering
 * agreement (ARI/NMI/ASW) is a meaningful comparison.
 */

include { SEURAT_LOAD       } from '../modules/seurat.nf'
include { SEURAT_QC         } from '../modules/seurat.nf'
include { SEURAT_SCTRANSFORM } from '../modules/seurat.nf'
include { SEURAT_HARMONY    } from '../modules/seurat.nf'
include { SEURAT_NEIGHBORS  } from '../modules/seurat.nf'
include { SEURAT_LEIDEN     } from '../modules/seurat.nf'
include { SEURAT_UMAP       } from '../modules/seurat.nf'
include { SEURAT_MARKERS    } from '../modules/seurat.nf'
include { SEURAT_ANNOTATE   } from '../modules/seurat.nf'

workflow SCRNASEQ_CPU {

    take:
    ch_samples         // (sample_id, matrix_dir)

    main:

    ch_versions = Channel.empty()

    SEURAT_LOAD(ch_samples)
    SEURAT_QC(SEURAT_LOAD.out.rds)
    SEURAT_SCTRANSFORM(SEURAT_QC.out.rds)

    // Merge all patients then integrate with Harmony for batch correction
    ch_merged = SEURAT_SCTRANSFORM.out.rds.collect()
    SEURAT_HARMONY(ch_merged)

    // Neighbors graph + parallel Leiden over multiple resolutions
    SEURAT_NEIGHBORS(SEURAT_HARMONY.out.rds)
    SEURAT_LEIDEN(SEURAT_NEIGHBORS.out.rds, params.sc_leiden_resolutions)
    SEURAT_UMAP(SEURAT_LEIDEN.out.rds)
    SEURAT_MARKERS(SEURAT_UMAP.out.rds)
    SEURAT_ANNOTATE(SEURAT_MARKERS.out.rds, file("${projectDir}/data/celltype_markers.yaml"))

    ch_versions = ch_versions
        .mix(SEURAT_LOAD.out.versions)
        .mix(SEURAT_QC.out.versions)
        .mix(SEURAT_SCTRANSFORM.out.versions)
        .mix(SEURAT_HARMONY.out.versions)
        .mix(SEURAT_LEIDEN.out.versions)
        .mix(SEURAT_MARKERS.out.versions)

    emit:
    rds       = SEURAT_ANNOTATE.out.rds       // final annotated Seurat object
    reference = SEURAT_ANNOTATE.out.reference // pseudobulk reference for deconvolution
    timing    = SEURAT_ANNOTATE.out.timing    // per-step wall-clock for Fig 3A
    versions  = ch_versions
}
