# Container manifest

Images used by KazRNA-Pipe v1.2.0. Every process in the workflow runs inside
one of these; none executes on the host.

Registry: `ghcr.io/kazrna-pipe`

Images are referenced throughout `conf/containers.config` by immutable digest
rather than by tag. A tag can be repointed at different content at any time; a
digest is a hash of the content and cannot. A run of this pipeline therefore
either obtains exactly the images below or fails.

| Image | Tag | Digest | Built (UTC) | Size | Base | Used by |
|---|---|---|---|---|---|---|
| kazrna-bulk | `1.1.0` | `sha256:449731ddf753317a5026161e7e05eb73d2e74186b425ca3e826204db390ab427` | 2026-08-05T06:18:27Z | 5,850,336,567 B (5.45 GiB) | bioconductor/bioconductor_docker:RELEASE_3_19 | DESEQ2, EDGER, LIMMA_VOOM, CONCORDANCE, CLUSTERPROFILER, TXIMPORT, MERGE_COUNTS, CELLTYPE_DE, BAYESPRISM, MUSIC, DECONV_AGREEMENT |
| kazrna-seurat | `1.1.0` | `sha256:034eefc888af1f3efe22bb65821040cead37c8995a16ed731652b0765ab7cf2b` | 2026-08-05T05:49:58Z | 4,757,210,117 B (4.43 GiB) | satijalab/seurat:5.0.0 | SEURAT_WORKFLOW |
| kazrna-py | `1.1.0` | `sha256:17bc730681ca90d649f1bf25840f4672ded381747338085e0e6645a1c63f6832` | 2026-08-05T05:28:32Z | 1,264,981,178 B (1.18 GiB) | python:3.11-slim | RAPIDS_SINGLECELL, CLUSTERING_AGREEMENT, SOFTWARE_VERSIONS |
| hisat2_samtools (third party) | `6be64e12472a7b75` | `sha256:c34a62b3e0c61c75a10f869eaf75061a464d2177c32bb697ac864db4039a2ff4` | 2025-06-02T12:27:33Z | 543,838,687 B (0.51 GiB) | Seqera Wave (hisat2 2.2.1 + samtools 1.20) | HISAT2_ALIGN |

## Archive

The software release these images accompany is archived at Zenodo:

- All versions: [10.5281/zenodo.21832631](https://doi.org/10.5281/zenodo.21832631)
- v1.1.0: [10.5281/zenodo.21832632](https://doi.org/10.5281/zenodo.21832632)
- v1.2.0: [NEW_DOI]

To verify an image against this manifest:

```bash
docker pull ghcr.io/kazrna-pipe/kazrna-bulk@sha256:449731ddf753317a5026161e7e05eb73d2e74186b425ca3e826204db390ab427
docker image inspect --format='{{index .RepoDigests 0}}' \
  ghcr.io/kazrna-pipe/kazrna-bulk@sha256:449731ddf753317a5026161e7e05eb73d2e74186b425ca3e826204db390ab427
```

The packages present in each image are verified at build time: every Dockerfile
ends with a layer that fails the build if a required package is absent, so an
image that builds is an image the pipeline can run in.

## Third-party images

Used unmodified from the BioContainers registry:

| Tool | Image | Used by |
|---|---|---|
| FastQC 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | FASTQC |
| Trim Galore 0.6.10 | `quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0` | TRIMGALORE |
| STAR 2.7.11a | `quay.io/biocontainers/star:2.7.11a--h0033a41_0` | STAR_ALIGN |
| HISAT2 2.2.1 | `quay.io/biocontainers/hisat2:2.2.1--h87f3376_5` | HISAT2_ALIGN |
| SAMtools 1.20 | `quay.io/biocontainers/samtools:1.20--h50ea8bc_0` | SAMTOOLS_SORT |
| Salmon 1.10.2 | `quay.io/biocontainers/salmon:1.10.2--hecfa306_0` | SALMON_QUANT |
| Subread 2.0.6 | `quay.io/biocontainers/subread:2.0.6--he4a0461_2` | FEATURECOUNTS |

The versions actually loaded at runtime are recorded per run in
`results/<run>/pipeline_info/software_versions.yml`, which is generated from
each process rather than transcribed by hand.

## CIBERSORTx

CIBERSORTx requires a Stanford-registered token and its image cannot be
redistributed. The process is skipped unless `params.cibersortx_token` is set,
in which case deconvolution proceeds with BayesPrism and MuSiC.

## Reproducing

```bash
bash containers/build.sh --version 1.1.0
bash containers/build.sh --version 1.1.0 --push ghcr.io/kazrna-pipe
```

Dockerfiles: `Dockerfile.bulk`, `Dockerfile.seurat`, `Dockerfile.py`.
