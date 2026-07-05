# Container manifest

This file is the canonical record of every container image used by **kazrna-pipe v1.0.0**. Each image is pinned by its **OCI digest** (the
content-addressable SHA-256 of the manifest). Tags are mutable and are listed
only as a convenience.

Pre-built images with their SHA-256 digests are archived at Zenodo
(DOI: [RESERVED_DOI]); the digest column is populated from that deposit.
Base images are already pinned by immutable OCI digest (linux/amd64) below,
so each image builds reproducibly from the recipes in this directory.

| Image           | Tag     | OCI digest (sha256)                      | Built | Size |
| --------------- | ------- | ---------------------------------------- | ----- | ---- |
| `kazrna/bulk`   | `1.0.0` | Deposited at Zenodo [RESERVED_DOI]       | —     | —    |
| `kazrna/seurat` | `1.0.0` | Deposited at Zenodo [RESERVED_DOI]       | —     | —    |
| `kazrna/rapids` | `1.0.0` | Deposited at Zenodo [RESERVED_DOI]       | —     | —    |

## Base images (pinned by digest)

| Upstream image                       | Digest                                                                   |
| ------------------------------------ | ------------------------------------------------------------------------ |
| `mambaorg/micromamba:1.5.8-jammy`    | `sha256:0d2870c1159dfb0c285e1c942c68385cd1110255cc9551838a29e48793e21af1` |
| `rocker/r-ver:4.4.1`                 | `sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312` |
| `rapidsai/base:24.06-cuda12.2-py3.11`| `sha256:8f0c5090edf797a31412c40d72ffd003c2a9a3a3e957529f0a044fc776c06a59` |

## Building the Singularity/Apptainer images

All three images are fully defined by the recipe files in this directory and
rebuild deterministically on any Linux host with Singularity ≥ 3.8 or
Apptainer ≥ 1.1. Each recipe pins its base image by the OCI digest listed above.

```
cd containers
singularity build --fakeroot kazrna-bulk_1.0.0.sif    Singularity.bulk
singularity build --fakeroot kazrna-seurat_1.0.0.sif  Singularity.seurat
singularity build --fakeroot kazrna-rapids_1.0.0.sif  Singularity.rapids

# If --fakeroot is unavailable, build with sudo, e.g.:
#   sudo singularity build kazrna-bulk_1.0.0.sif Singularity.bulk
# `apptainer` may be substituted for `singularity` on newer systems.
```

After building, record each image's SHA-256 digest and paste it into the
OCI digest column above:

```
for sif in kazrna-bulk_1.0.0.sif kazrna-seurat_1.0.0.sif kazrna-rapids_1.0.0.sif ; do
  echo -n "${sif}: "
  sha256sum "${sif}"
done
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
results figure (Fig 4B) will then show two methods rather than three.