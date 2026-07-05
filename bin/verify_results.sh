#!/usr/bin/env bash
# bin/verify_results.sh
#
# Sanity-checks that every expected output of a successful test-profile run is
# present and non-empty. Used by the GitHub Actions CI to gate merges.
#
# Usage:
#   bash bin/verify_results.sh <outdir>
#
# Exit codes:
#   0  all expected outputs present and non-empty
#   1  one or more outputs missing or empty
#   2  argument error

set -uo pipefail

OUT="${1:-results_ci}"
if [[ -z "${OUT}" ]]; then
    echo "Usage: $0 <outdir>" >&2; exit 2
fi
if [[ ! -d "${OUT}" ]]; then
    echo "Output dir '${OUT}' does not exist" >&2; exit 1
fi

declare -a expected=(
    # Per-aligner gene counts
    "quant/featurecounts"             # one .counts.tsv per sample
    "quant/tximport/salmon.gene_counts.tsv"
    # DE results
    "de/deseq2/star"
    "de/deseq2/hisat2"
    "de/deseq2/salmon"
    "de/edger/star"
    "de/limma_voom/star"
    # Concordance
    "concordance/fig2a_spearman_matrix.tsv"
    "concordance/fig2b_deg_overlap.tsv"
    # scRNA-seq
    "scrnaseq/cpu"
    "scrnaseq/agreement/agreement_metrics.tsv"
    # Deconvolution
    "deconvolution/bayesprism/bayesprism_theta.tsv"
    "deconvolution/music/music_proportions.tsv"
    "deconvolution/agreement/deconv_correlation.tsv"
    # Software versions
    "pipeline_info/software_versions.yml"
    "pipeline_info/software_versions.tsv"
)

missing=0
for item in "${expected[@]}"; do
    path="${OUT}/${item}"
    if [[ -d "${path}" ]]; then
        # Directory must be non-empty
        if [[ -z "$(ls -A "${path}" 2>/dev/null)" ]]; then
            echo "[MISSING] empty dir: ${path}"
            missing=$((missing + 1))
        else
            echo "[OK] dir: ${path} ($(find "${path}" -type f | wc -l) files)"
        fi
    elif [[ -f "${path}" ]]; then
        if [[ ! -s "${path}" ]]; then
            echo "[MISSING] empty file: ${path}"
            missing=$((missing + 1))
        else
            echo "[OK] file: ${path} ($(wc -c < "${path}") bytes)"
        fi
    else
        echo "[MISSING] absent: ${path}"
        missing=$((missing + 1))
    fi
done

# Spot-check: software_versions.yml must list at least the major tools we expect.
versions="${OUT}/pipeline_info/software_versions.yml"
if [[ -f "${versions}" ]]; then
    for tool in STAR hisat2 salmon DESeq2 edgeR limma Seurat; do
        if ! grep -qi "${tool}" "${versions}"; then
            echo "[MISSING] ${tool} not in software_versions.yml"
            missing=$((missing + 1))
        fi
    done
fi

if [[ ${missing} -gt 0 ]]; then
    echo ""
    echo "VERIFY FAILED: ${missing} expected outputs missing or empty"
    exit 1
fi

echo ""
echo "VERIFY OK: all ${#expected[@]} expected outputs present and non-empty"
exit 0
