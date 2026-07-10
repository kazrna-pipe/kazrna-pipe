/*
 * Single-cell RNA-seq workflow (CPU path, Seurat).
 *
 * Loads each 10x matrix, runs the Seurat clustering and annotation pipeline,
 * and exposes the cell-type-labelled object used as the deconvolution reference.
 */

include { SEURAT_WORKFLOW } from '../modules/seurat.nf'

workflow SCRNASEQ_CPU {

    take:
    sc_samples      // (meta, matrix)
    marker_yaml

    main:
    SEURAT_WORKFLOW(sc_samples, marker_yaml)

    emit:
    reference = SEURAT_WORKFLOW.out.seurat
    celltypes = SEURAT_WORKFLOW.out.celltypes
    clusters  = SEURAT_WORKFLOW.out.clusters
    versions  = SEURAT_WORKFLOW.out.versions
}
