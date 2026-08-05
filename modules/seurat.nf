// modules/seurat.nf
// Wraps scripts/R/seurat_workflow.R - CPU single-cell pipeline.
//
// Runs ONCE over all samples: the script merges them and Harmony-integrates
// across sample_id, so a per-sample invocation would make integration a no-op.

process SEURAT_WORKFLOW {
    label      'process_high_memory'
    publishDir "${params.outdir}/scrnaseq/cpu", mode: 'copy'

    input:
    path  h5_files
    val   sample_ids
    path  marker_yaml

    output:
    path "seurat_object.rds", emit: seurat
    path "labels.tsv",        emit: clusters
    path "markers.tsv",       emit: markers
    path "umap_coords.tsv",   emit: umap
    path "timing.tsv",        emit: timings
    path "provenance.json",   emit: provenance
    path "versions.yml",      emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/seurat_workflow.R \\
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
        Seurat:       \$(Rscript -e 'cat(as.character(packageVersion("Seurat")))')
        SeuratObject: \$(Rscript -e 'cat(as.character(packageVersion("SeuratObject")))')
        harmony:      \$(Rscript -e 'cat(as.character(packageVersion("harmony")))')
    VER
    """
}
