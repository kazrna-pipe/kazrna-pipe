/*
 * Single-cell RNA-seq workflow (GPU path, rapids-singlecell).
 *
 * GPU-accelerated equivalent of the CPU path, used for the CPU-vs-GPU
 * benchmarking comparison.
 *
 * The module runs once over all samples rather than once per sample: the
 * script merges the samples and Harmony-integrates across sample_id, so a
 * per-sample invocation would leave nothing for the integration step to do.
 */

include { RAPIDS_SINGLECELL } from '../modules/rapids_singlecell.nf'

workflow SCRNASEQ_GPU {

    take:
    sc_samples
    marker_yaml

    main:
    ch_h5  = sc_samples.map { meta, h5 -> h5 }.collect()
    ch_ids = sc_samples.map { meta, h5 -> meta.id }.collect().map { it.join(',') }

    RAPIDS_SINGLECELL(ch_h5, ch_ids, marker_yaml)

    emit:
    adata     = RAPIDS_SINGLECELL.out.adata
    clusters  = RAPIDS_SINGLECELL.out.clusters
    markers   = RAPIDS_SINGLECELL.out.markers
    umap      = RAPIDS_SINGLECELL.out.umap
    timings   = RAPIDS_SINGLECELL.out.timings
    versions  = RAPIDS_SINGLECELL.out.versions
}
