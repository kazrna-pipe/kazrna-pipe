#!/usr/bin/env bash
# containers/build.sh
#
# Builds the three KazRNA-Pipe images, verifies each one contains the packages
# the scripts import, and writes containers/MANIFEST.md with real digests,
# build dates and sizes.
#
#   bash containers/build.sh                 # build locally, record local IDs
#   bash containers/build.sh --push ghcr.io/ORG   # also push, record digests


set -euo pipefail

VERSION="1.1.0"
REGISTRY=""
PUSH=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --push)    PUSH=1; REGISTRY="$2"; shift 2 ;;
        --version) VERSION="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")"
BUILD_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

IMAGES=("bulk" "seurat" "py")

echo "=============================================================="
echo " Building KazRNA-Pipe containers v${VERSION}"
echo " Build date: ${BUILD_DATE}"
echo "=============================================================="

for img in "${IMAGES[@]}"; do
    tag="kazrna-${img}:${VERSION}"
    echo
    echo "--- Building ${tag}"
    docker build \
        --file "Dockerfile.${img}" \
        --tag "${tag}" \
        --label "org.opencontainers.image.version=${VERSION}" \
        --label "org.opencontainers.image.created=${BUILD_DATE}" \
        .
done

# --- Verify the images can actually run the scripts ------------------------
echo
echo "--- Verifying package availability"
docker run --rm "kazrna-bulk:${VERSION}" \
    Rscript -e "library(optparse); library(DESeq2); library(BayesPrism); cat('bulk OK\n')"
docker run --rm "kazrna-seurat:${VERSION}" \
    Rscript -e "library(optparse); library(Seurat); library(hdf5r); library(harmony); cat('seurat OK\n')"
docker run --rm "kazrna-py:${VERSION}" \
    python -c "import scanpy, sklearn, yaml; print('py OK')"

# --- Push and collect digests ----------------------------------------------
declare -A DIGESTS
if [[ ${PUSH} -eq 1 ]]; then
    echo
    echo "--- Pushing to ${REGISTRY}"
    for img in "${IMAGES[@]}"; do
        local_tag="kazrna-${img}:${VERSION}"
        remote_tag="${REGISTRY}/kazrna-${img}:${VERSION}"
        docker tag "${local_tag}" "${remote_tag}"
        docker push "${remote_tag}"
        DIGESTS[$img]=$(docker inspect --format='{{index .RepoDigests 0}}' "${remote_tag}" | cut -d@ -f2)
        echo "    ${img}: ${DIGESTS[$img]}"
    done
else
    for img in "${IMAGES[@]}"; do
        DIGESTS[$img]="(local build - push to obtain a registry digest)"
    done
fi

# --- Write the manifest -----------------------------------------------------
echo
echo "--- Writing MANIFEST.md"
{
    echo "# Container manifest"
    echo
    echo "Images used by KazRNA-Pipe v${VERSION}. Every process in the workflow"
    echo "runs inside one of these; none executes on the host."
    echo
    echo "Built: \`${BUILD_DATE}\`"
    if [[ ${PUSH} -eq 1 ]]; then
        echo "Registry: \`${REGISTRY}\`"
    else
        echo
        echo "> **These are local builds.** Push them to a registry before"
        echo "> submission: only a pushed image has an immutable digest that a"
        echo "> reader can pull and verify."
    fi
    echo
    echo "| Image | Tag | Digest | Size | Base | Used by |"
    echo "|---|---|---|---|---|---|"

    size_bulk=$(docker image inspect "kazrna-bulk:${VERSION}" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')
    size_seur=$(docker image inspect "kazrna-seurat:${VERSION}" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')
    size_py=$(docker image inspect "kazrna-py:${VERSION}" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')

    echo "| kazrna-bulk | \`${VERSION}\` | \`${DIGESTS[bulk]}\` | ${size_bulk} | bioconductor/bioconductor_docker:RELEASE_3_19 | DESEQ2, EDGER, LIMMA_VOOM, CONCORDANCE, CLUSTERPROFILER, TXIMPORT, MERGE_COUNTS, CELLTYPE_DE, BAYESPRISM, MUSIC, DECONV_AGREEMENT |"
    echo "| kazrna-seurat | \`${VERSION}\` | \`${DIGESTS[seurat]}\` | ${size_seur} | satijalab/seurat:5.0.0 | SEURAT_WORKFLOW |"
    echo "| kazrna-py | \`${VERSION}\` | \`${DIGESTS[py]}\` | ${size_py} | python:3.11-slim | RAPIDS_SINGLECELL, CLUSTERING_AGREEMENT, SOFTWARE_VERSIONS |"
    echo
    echo "## Third-party images"
    echo
    echo "Used unmodified from the BioContainers registry:"
    echo
    echo "| Tool | Image | Used by |"
    echo "|---|---|---|"
    echo "| FastQC 0.12.1 | \`quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0\` | FASTQC |"
    echo "| Trim Galore 0.6.10 | \`quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0\` | TRIMGALORE |"
    echo "| STAR 2.7.11b | \`quay.io/biocontainers/star:2.7.11b--h43eeafb_1\` | STAR_ALIGN |"
    echo "| HISAT2 2.2.1 | \`quay.io/biocontainers/hisat2:2.2.1--h87f3376_5\` | HISAT2_ALIGN |"
    echo "| SAMtools 1.20 | \`quay.io/biocontainers/samtools:1.20--h50ea8bc_0\` | SAMTOOLS_SORT |"
    echo "| Salmon 1.10.2 | \`quay.io/biocontainers/salmon:1.10.2--hecfa306_0\` | SALMON_QUANT |"
    echo "| Subread 2.0.6 | \`quay.io/biocontainers/subread:2.0.6--he4a0461_2\` | FEATURECOUNTS |"
    echo
    echo "## CIBERSORTx"
    echo
    echo "CIBERSORTx requires a Stanford-registered token and its image cannot"
    echo "be redistributed. The process is skipped unless"
    echo "\`params.cibersortx_token\` is set."
    echo
    echo "## Reproducing"
    echo
    echo '```bash'
    echo "bash containers/build.sh --version ${VERSION}"
    echo "bash containers/build.sh --version ${VERSION} --push ghcr.io/ORG"
    echo '```'
    echo
    echo "Dockerfiles: \`Dockerfile.bulk\`, \`Dockerfile.seurat\`, \`Dockerfile.py\`."
    echo "Each ends with a verification layer that fails the build if any"
    echo "required package is absent, so an image that builds is an image the"
    echo "pipeline can run in."
} > MANIFEST.md

echo
echo "=============================================================="
echo " Done. Wrote containers/MANIFEST.md"
if [[ ${PUSH} -eq 0 ]]; then
    echo
    echo " Next: update conf/containers.config to use kazrna-*:${VERSION},"
    echo " run the pipeline, then re-run with --push to record real digests."
fi
echo "=============================================================="
