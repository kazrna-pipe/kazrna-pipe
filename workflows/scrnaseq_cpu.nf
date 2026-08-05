/*
 * Single-cell RNA-seq workflow (CPU path, Seurat).
 *
 * Loads every 10x matrix, runs the Seurat clustering pipeline, and exposes the
 * integrated object used as the deconvolution reference.
 *
 * The module runs once over all samples rather than once per sample: the
 * script merges the samples and Harmony-integrates across sample_id, so a
 * per-sample invocation would leave nothing for the integration step to do.
 */

include { SEURAT_WORKFLOW } from '../modules/seurat.nf'

workflow SCRNASEQ_CPU {

    take:
    sc_samples
    marker_yaml

    main:
    ch_h5  = sc_samples.map { meta, h5 -> h5 }.collect()
    ch_ids = sc_samples.map { meta, h5 -> meta.id }.collect().map { it.join(',') }

    SEURAT_WORKFLOW(ch_h5, ch_ids, marker_yaml)

    emit:
    reference = SEURAT_WORKFLOW.out.seurat
    clusters  = SEURAT_WORKFLOW.out.clusters
    markers   = SEURAT_WORKFLOW.out.markers
    umap      = SEURAT_WORKFLOW.out.umap
    timings   = SEURAT_WORKFLOW.out.timings
    versions  = SEURAT_WORKFLOW.out.versions
}
