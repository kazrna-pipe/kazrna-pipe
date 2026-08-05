// modules/rapids_singlecell.nf
// Wraps scripts/python/rapids_singlecell_workflow.py - GPU single-cell pipeline.
//
// Runs ONCE over all samples: the script merges them and Harmony-integrates
// across sample_id, so a per-sample invocation would make integration a no-op.

process RAPIDS_SINGLECELL {
    label      'process_gpu'
    publishDir "${params.outdir}/scrnaseq/gpu", mode: 'copy'

    input:
    path  h5_files
    val   sample_ids
    path  marker_yaml

    output:
    path "adata.h5ad",       emit: adata
    path "labels.tsv",       emit: clusters
    path "markers.tsv",      emit: markers
    path "umap_coords.tsv",  emit: umap
    path "timing.tsv",       emit: timings
    path "provenance.json",  emit: provenance
    path "versions.yml",     emit: versions

    script:
    """
    python ${projectDir}/scripts/python/rapids_singlecell_workflow.py \\
        --input_h5      ${h5_files.join(',')} \\
        --sample_ids    ${sample_ids} \\
        --markers       ${marker_yaml} \\
        --hvg_n         ${params.sc_hvg_n} \\
        --n_pcs         ${params.sc_n_pcs} \\
        --resolutions   ${params.sc_leiden_resolutions} \\
        --normalization ${params.sc_normalization} \\
        --mt_threshold  ${params.sc_mt_threshold} \\
        --min_features  ${params.sc_min_genes} \\
        --threads       ${task.cpus} \\
        --outdir        .

    cat <<-VER > versions.yml
    "${task.process}":
        rapids_singlecell: \$(python -c 'import rapids_singlecell as r; print(r.__version__)' 2>/dev/null || echo "not available")
        scanpy:            \$(python -c 'import scanpy; print(scanpy.__version__)')
        anndata:           \$(python -c 'import anndata; print(anndata.__version__)')
    VER
    """
}
