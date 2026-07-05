# Container manifest

This file is the canonical record of every container image used by
**kazrna-pipe v1.0.0**. Each image is pinned by its **OCI digest** (the
content-addressable SHA-256 of the manifest). Tags are mutable and are listed
only as a convenience.


| Image                   | Tag      | OCI digest (sha256)               | Built       | Size  |
|-------------------------|----------|-----------------------------------|-------------|-------|
| `kazrna/bulk`           | `1.0.0`  | TBA — populated at v1.0.0 release | TBA         | TBA   |
| `kazrna/seurat`         | `1.0.0`  | TBA — populated at v1.0.0 release | TBA         | TBA   |
| `kazrna/rapids`         | `1.0.0`  | TBA — populated at v1.0.0 release | TBA         | TBA   |

## Base images (pinned by digest once the build host is decided)

| Upstream image                                          | Digest |
|---------------------------------------------------------|--------|
| `mambaorg/micromamba:1.5.8-jammy`                       | TBA    |
| `rocker/r-ver:4.4.1`                                    | TBA    |
| `nvcr.io/nvidia/rapidsai/notebooks:24.06-cuda12.2-py3.11` | TBA  |

> **How to fill these in.** After running `singularity build` for each image,
> capture the digest with:
>
> ```bash
> for sif in kazrna-bulk_1.0.0.sif kazrna-seurat_1.0.0.sif kazrna-rapids_1.0.0.sif ; do
>     echo -n "${sif}: "
>     sha256sum "${sif}"
> done
> ```
>
> Paste each `sha256:<hex>` value into the OCI digest column above.

## Building

```bash
cd containers
singularity build kazrna-bulk_1.0.0.sif    Singularity.bulk
singularity build kazrna-seurat_1.0.0.sif  Singularity.seurat
singularity build kazrna-rapids_1.0.0.sif  Singularity.rapids
```

## CIBERSORTx

CIBERSORTx is **not** redistributed; its licence prohibits this. Users that
wish to enable that arm of the deconvolution comparison must:

1. Register at <https://cibersortx.stanford.edu/> and request a token.
2. Pull the official image: `docker pull cibersortx/fractions:latest`
3. Convert to Singularity:
   `singularity build cibersortx_fractions.sif docker-daemon://cibersortx/fractions:latest`
4. Set `params.cibersortx_token` and `params.cibersortx_email`.

The pipeline runs all other analyses without CIBERSORTx; the deconvolution
results figure (Fig 5B) will then show two methods rather than three.
