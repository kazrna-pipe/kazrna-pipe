/*
 * Single-cell RNA-seq workflow (GPU path, rapids-singlecell).
 *
 * GPU-accelerated equivalent of the CPU path, used for the CPU-vs-GPU
 * benchmarking comparison (Figure 3).
 */

include { RAPIDS_SINGLECELL } from '../modules/rapids_singlecell.nf'

workflow SCRNASEQ_GPU {

    take:
    sc_samples      // (meta, matrix)
    marker_yaml

    main:
    RAPIDS_SINGLECELL(sc_samples, marker_yaml)

    emit:
    adata     = RAPIDS_SINGLECELL.out.adata
    celltypes = RAPIDS_SINGLECELL.out.celltypes
    clusters  = RAPIDS_SINGLECELL.out.clusters
    versions  = RAPIDS_SINGLECELL.out.versions
}
