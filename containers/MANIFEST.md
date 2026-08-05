# Container manifest

Images used by KazRNA-Pipe v1.1.0. Every process in the workflow
runs inside one of these; none executes on the host.

Built: `2026-08-05T20:33:50Z`
Registry: `ghcr.io/kazrna-pipe`

| Image | Tag | Digest | Size | Base | Used by |
|---|---|---|---|---|---|
| kazrna-bulk | `1.1.0` | `sha256:449731ddf753317a5026161e7e05eb73d2e74186b425ca3e826204db390ab427` | 5.4 GB | bioconductor/bioconductor_docker:RELEASE_3_19 | DESEQ2, EDGER, LIMMA_VOOM, CONCORDANCE, CLUSTERPROFILER, TXIMPORT, MERGE_COUNTS, CELLTYPE_DE, BAYESPRISM, MUSIC, DECONV_AGREEMENT |
| kazrna-seurat | `1.1.0` | `sha256:034eefc888af1f3efe22bb65821040cead37c8995a16ed731652b0765ab7cf2b` | 4.4 GB | satijalab/seurat:5.0.0 | SEURAT_WORKFLOW |
| kazrna-py | `1.1.0` | `sha256:17bc730681ca90d649f1bf25840f4672ded381747338085e0e6645a1c63f6832` | 1.2 GB | python:3.11-slim | RAPIDS_SINGLECELL, CLUSTERING_AGREEMENT, SOFTWARE_VERSIONS |

## Third-party images

Used unmodified from the BioContainers registry:

| Tool | Image | Used by |
|---|---|---|
| FastQC 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | FASTQC |
| Trim Galore 0.6.10 | `quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0` | TRIMGALORE |
| STAR 2.7.11b | `quay.io/biocontainers/star:2.7.11b--h43eeafb_1` | STAR_ALIGN |
| HISAT2 2.2.1 | `quay.io/biocontainers/hisat2:2.2.1--h87f3376_5` | HISAT2_ALIGN |
| SAMtools 1.20 | `quay.io/biocontainers/samtools:1.20--h50ea8bc_0` | SAMTOOLS_SORT |
| Salmon 1.10.2 | `quay.io/biocontainers/salmon:1.10.2--hecfa306_0` | SALMON_QUANT |
| Subread 2.0.6 | `quay.io/biocontainers/subread:2.0.6--he4a0461_2` | FEATURECOUNTS |

## CIBERSORTx

CIBERSORTx requires a Stanford-registered token and its image cannot
be redistributed. The process is skipped unless
`params.cibersortx_token` is set.

## Reproducing

```bash
bash containers/build.sh --version 1.1.0
bash containers/build.sh --version 1.1.0 --push ghcr.io/ORG
```

Dockerfiles: `Dockerfile.bulk`, `Dockerfile.seurat`, `Dockerfile.py`.
Each ends with a verification layer that fails the build if any
required package is absent, so an image that builds is an image the
pipeline can run in.
