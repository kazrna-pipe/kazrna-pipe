// modules/clustering_agreement.nf
// Wraps scripts/python/clustering_agreement.py - Figure 3C (ARI/NMI/ASW
// between CPU Seurat and GPU rapids-singlecell at matched resolutions).

process CLUSTERING_AGREEMENT {
    label      'process_low'
    publishDir "${params.outdir}/scrnaseq/agreement", mode: 'copy'

    input:
    path cpu_clusters    // collected clusters.tsv from Seurat
    path gpu_clusters    // collected clusters.tsv from rapids

    output:
    path "agreement_metrics.tsv",      emit: metrics
    path "fig3c_ari_curve.pdf",        emit: plot
    path "agreement_provenance.json",  emit: provenance
    path "versions.yml",               emit: versions

    script:
    """
    python ${projectDir}/scripts/python/clustering_agreement.py \\
        --cpu_dir   . \\
        --gpu_dir   . \\
        --out_dir   . \\
        --self_check 0.99

    cat <<-VER > versions.yml
    "${task.process}":
        scikit-learn: \$(python -c 'import sklearn; print(sklearn.__version__)')
        scanpy:       \$(python -c 'import scanpy; print(scanpy.__version__)')
    VER
    """
}
