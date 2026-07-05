// modules/software_versions.nf
// Aggregates the per-process versions.yml emissions into one report that is
// always shipped alongside the results. Cross-referenced in the manuscript's
// Methods section and the GitHub release notes.

process SOFTWARE_VERSIONS {
    label      'process_single'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path versions_files

    output:
    path "software_versions.yml",      emit: yml
    path "software_versions.tsv",      emit: tsv
    path "software_versions.html",     emit: html

    script:
    """
    python ${projectDir}/scripts/python/collect_versions.py \\
        --inputs ${versions_files} \\
        --out_yml  software_versions.yml \\
        --out_tsv  software_versions.tsv \\
        --out_html software_versions.html
    """
}
