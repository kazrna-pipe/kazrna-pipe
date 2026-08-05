// modules/clustering_agreement.nf
// Wraps scripts/python/clustering_agreement.py - Figure 3C (ARI/NMI/ASW
// between CPU Seurat and GPU rapids-singlecell at matched resolutions).
//
// Both paths emit a file called labels.tsv, so they are staged under distinct
// names to stop Nextflow silently renaming one of them.

process CLUSTERING_AGREEMENT {
    label      'process_low'
    publishDir "${params.outdir}/scrnaseq/agreement", mode: 'copy'

    input:
    path(cpu_clusters, stageAs: 'cpu_labels.tsv')
    path(gpu_clusters, stageAs: 'gpu_labels.tsv')

    output:
    path "ari_nmi_asw_by_resolution.tsv", emit: metrics
    path "ari_nmi_asw_by_resolution.pdf", emit: plot
    path "provenance.json",               emit: provenance
    path "versions.yml",                  emit: versions

    script:
    """
    python ${projectDir}/scripts/python/clustering_agreement.py \\
        --cpu_labels cpu_labels.tsv \\
        --gpu_labels gpu_labels.tsv \\
        --outdir     . \\
        --self_check

    cat <<-VER > versions.yml
    "${task.process}":
        scikit-learn: \$(python -c 'import sklearn; print(sklearn.__version__)')
        scanpy:       \$(python -c 'import scanpy; print(scanpy.__version__)')
    VER
    """
}
