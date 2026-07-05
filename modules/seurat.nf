// modules/seurat.nf
// Wraps scripts/R/seurat_workflow.R - CPU single-cell pipeline.

process SEURAT_WORKFLOW {
    tag        "${meta.id}"
    label      'process_high_memory'
    publishDir "${params.outdir}/scrnaseq/cpu/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(matrix_dir)
    path  marker_yaml

    output:
    tuple val(meta), path("seurat_object.rds"),       emit: seurat
    tuple val(meta), path("clusters.tsv"),            emit: clusters
    tuple val(meta), path("celltype_assignments.tsv"),emit: celltypes
    tuple val(meta), path("stage_timings.tsv"),       emit: timings
    path "seurat_provenance.json",                    emit: provenance
    path "seurat_session_info.txt",                   emit: session
    path "versions.yml",                              emit: versions

    script:
    """
    Rscript ${projectDir}/scripts/R/seurat_workflow.R \\
        --counts_dir    ${matrix_dir} \\
        --sample_id     ${meta.id} \\
        --markers       ${marker_yaml} \\
        --n_features    ${params.sc_n_hvgs} \\
        --n_pcs         ${params.sc_n_pcs} \\
        --resolutions   ${params.sc_resolutions} \\
        --threads       ${task.cpus} \\
        --out_dir       .

    cat <<-VER > versions.yml
    "${task.process}":
        Seurat:           \$(Rscript -e 'cat(as.character(packageVersion("Seurat")))')
        SeuratObject:     \$(Rscript -e 'cat(as.character(packageVersion("SeuratObject")))')
        glmGamPoi:        \$(Rscript -e 'cat(as.character(packageVersion("glmGamPoi")))')
    VER
    """
}
