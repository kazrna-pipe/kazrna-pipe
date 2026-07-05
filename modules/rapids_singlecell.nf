// modules/rapids_singlecell.nf
// Wraps scripts/python/rapids_singlecell_workflow.py - GPU single-cell pipeline.

process RAPIDS_SINGLECELL {
    tag        "${meta.id}"
    label      'process_gpu'
    publishDir "${params.outdir}/scrnaseq/gpu/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(matrix_dir)
    path  marker_yaml

    output:
    tuple val(meta), path("adata.h5ad"),              emit: adata
    tuple val(meta), path("clusters.tsv"),            emit: clusters
    tuple val(meta), path("celltype_assignments.tsv"),emit: celltypes
    tuple val(meta), path("stage_timings.tsv"),       emit: timings
    tuple val(meta), path("gpu_memory_trace.tsv"),    emit: gpu_trace
    path "rapids_provenance.json",                    emit: provenance
    path "versions.yml",                              emit: versions

    script:
    """
    python ${projectDir}/scripts/python/rapids_singlecell_workflow.py \\
        --counts_dir    ${matrix_dir} \\
        --sample_id     ${meta.id} \\
        --markers       ${marker_yaml} \\
        --n_features    ${params.sc_n_hvgs} \\
        --n_pcs         ${params.sc_n_pcs} \\
        --resolutions   ${params.sc_resolutions} \\
        --out_dir       .

    cat <<-VER > versions.yml
    "${task.process}":
        rapids_singlecell: \$(python -c 'import rapids_singlecell as r; print(r.__version__)')
        cudf:              \$(python -c 'import cudf; print(cudf.__version__)')
        cuml:              \$(python -c 'import cuml; print(cuml.__version__)')
        scanpy:            \$(python -c 'import scanpy; print(scanpy.__version__)')
        anndata:           \$(python -c 'import anndata; print(anndata.__version__)')
    VER
    """
}
