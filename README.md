# KazRNA-Pipe

[![CI](https://github.com/kazrna-pipe/kazrna-pipe/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/kazrna-pipe/kazrna-pipe/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21832631.svg)](https://doi.org/10.5281/zenodo.21832631)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**GPU-Accelerated, Reproducible Nextflow Workflow for Integrated Bulk and Single-Cell Transcriptomic Profiling**

KazRNA-Pipe is a Nextflow DSL2 pipeline that processes bulk and single-cell RNA-seq data through a unified, containerised, GPU-aware workflow. It was developed for population-scale precision oncology research on Kazakhstan academic HPC infrastructure and applied to a Kazakhstani esophageal squamous cell carcinoma (ESCC) cohort (PRJNA608223, n = 22) and the largest public ESCC single-cell atlas (GSE160269).

---

## Quick links

- **Reproducibility guide:** [REPRODUCIBILITY.md](REPRODUCIBILITY.md) — start here if you are a reviewer.
- **Data provenance:** [data/README.md](data/README.md) — accessions, download commands, metadata snapshots.
- **Container digests:** [containers/MANIFEST.md](containers/MANIFEST.md) — every image pinned by SHA-256 digest.
- **Test fixture:** [test_data/README.md](test_data/README.md) — what the CI fixture is, and what it is not.
- **Issue tracker:** https://github.com/kazrna-pipe/kazrna-pipe/issues

---

## Contents

1. [What the pipeline does](#what-the-pipeline-does)
2. [Repository layout](#repository-layout)
3. [Installation](#installation)
4. [Running the test profile](#running-the-test-profile)
5. [How reproducibility is enforced](#how-reproducibility-is-enforced)
6. [Running on real data](#running-on-real-data)
7. [Pinned versions](#pinned-versions)
8. [Data availability and provenance](#data-availability-and-provenance)
9. [Citation](#citation)
10. [License](#license)

---

## What the pipeline does

KazRNA-Pipe is organised into four modules, each independently runnable:

| Module | Inputs | Tools | Outputs |
| --- | --- | --- | --- |
| **Bulk RNA-seq** | Paired-end FASTQ + samplesheet | FastQC, Trim Galore, STAR, HISAT2, Salmon, featureCounts, tximport | Gene × sample count matrices, QC reports |
| **Differential expression** | Count matrix + design | DESeq2, edgeR, limma-voom (parallel) | DEG tables, cross-method concordance, KEGG/GO enrichment |
| **scRNA-seq** | 10x matrices | Seurat v5 (CPU) **or** rapids-singlecell (GPU) | Cluster labels, UMAP, marker genes, per-stage timings |
| **Deconvolution** | Bulk counts + sc reference | CIBERSORTx, BayesPrism, MuSiC (parallel) | Cell-type proportions, per-compartment DEGs |

Plus an HPC benchmarking module (strong/weak scaling, GPU utilisation, energy via RAPL and nvidia-smi).

Both aligners are quantified by the **same** featureCounts process with the same
settings, so cross-aligner concordance reflects the aligners rather than
differing count parameters. The CPU and GPU single-cell paths share a
normalisation method and the same `leidenalg` build, so their comparison
isolates the compute device.

---

## Repository layout

```
kazrna-pipe/
├── README.md                  This file
├── REPRODUCIBILITY.md         Step-by-step reproduction
├── CHANGELOG.md               Version history, including results affected by each change
├── CITATION.cff               Machine-readable citation
├── nextflow.config            Top-level Nextflow config
├── main.nf                    Pipeline entry point
├── conf/                      Execution profiles and container digests
├── modules/                   Per-tool Nextflow process definitions
├── workflows/                 Module-level workflows
├── bin/                       Data download and reference preparation scripts
├── scripts/R/                 R analysis scripts
├── scripts/python/            Python scripts, including the determinism checker
├── tools/                     Development checks (module/script interface contract)
├── containers/                Dockerfiles, build script, MANIFEST.md with image digests
├── data/                      Sample sheets, accession lists, metadata snapshots
├── test_data/                 Synthetic fixture for CI (~5 MB) and its generators
├── docs/                      Architecture diagram, parameter reference, benchmarking method
└── .github/workflows/ci.yml   Continuous integration
```

---

## Installation

| Component | Required | Notes |
| --- | --- | --- |
| Nextflow | ≥ 24.04 | `curl -s https://get.nextflow.io \| bash` — 24.04 introduced the `resourceLimits` directive this pipeline uses |
| Java | 17+ | Required by Nextflow |
| Docker or Apptainer | either | Both are exercised by CI on every push |
| CUDA | 12.2 | GPU profile only. The GPU path falls back to Scanpy when `rapids_singlecell` cannot be imported, so it can be tested without a GPU |

```bash
git clone https://github.com/kazrna-pipe/kazrna-pipe
cd kazrna-pipe
```

---

## Running the test profile

```bash
nextflow run main.nf -profile test,docker
```

Completes in roughly **90 seconds** on a 2-core machine with no GPU.

The test profile runs on a **synthetic fixture**, not on a subset of real data.
`test_data/` contains 200 artificial genes on a 304 kb artificial contig, with
paired-end reads generated *in silico* by the scripts in that directory. Forty
of the 200 genes carry a defined fold change between two arbitrary sample
groups, so the differential-expression steps produce non-empty results and CI
can check that the planted signal is recovered. No patient-derived sequence is
redistributed in this repository. See [test_data/README.md](test_data/README.md).

Expected outputs under `results/<outdir>/`:

```
quant/merged/star.gene_counts.tsv          Merged STAR counts
quant/merged/hisat2.gene_counts.tsv        Merged HISAT2 counts
quant/tximport/salmon.gene_counts.tsv      Salmon counts via tximport
de/deseq2/star/deseq2_results.tsv          Differential expression, three methods
de/edger/star/edger_results.tsv
de/limma_voom/star/limma_voom_results.tsv
concordance/spearman_matrix.tsv            Cross-aligner correlation
concordance/method_summary.tsv             DEG counts per method and three-way core
concordance/venn_DESeq2_edgeR_limma.pdf
scrnaseq/cpu/labels.tsv                    Cluster labels at each resolution
scrnaseq/cpu/umap_coords.tsv
scrnaseq/cpu/timing.tsv                    Per-stage wall time
pipeline_info/software_versions.yml        Every tool version actually invoked
```

---

## How reproducibility is enforced

Four mechanisms, all checked automatically on every push:

**Interface contract.** `tools/check_interfaces.py` verifies that each Nextflow module passes only options its script defines, and declares only outputs its script writes. It runs as the first CI job.

**Container pinning.** Every process declares a container. The three analysis images are referenced by SHA-256 digest in `conf/containers.config`, and CI fails if any is pinned by tag instead. Digests are listed in [containers/MANIFEST.md](containers/MANIFEST.md).

**End-to-end execution under both engines.** CI runs the complete pipeline on the fixture under Docker and Apptainer, asserts the expected outputs exist and are non-empty, and checks that the planted differential signal is recovered.

**Determinism.** Two independent runs are compared numerically by `scripts/python/diff_results.py`. Counts and all differential-expression output are bit-reproducible. Salmon's effective length and TPM are excluded, because `--gcBias` and `--seqBias` fit their models on a read subsample; the exclusion and its reasoning are recorded in the script, and `NumReads` — the column every downstream analysis consumes — is still compared.

---

## Running on real data

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the full walkthrough.

```bash
# 1. Fetch the cohort
bash bin/sra_fetch.sh --accession PRJNA608223 --target data/raw/PRJNA608223

# 2. Prepare references (GRCh38.p14 + GENCODE v44)
bash bin/prepare_references.sh --target data/refs

# 3. Run
nextflow run main.nf \
    -profile slurm,singularity \
    --input    data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --refs_dir data/refs \
    --outdir   results/manuscript
```

Building the STAR index for GRCh38 requires roughly 32 GB of RAM, so reference preparation is not feasible on a small workstation.

---

## Pinned versions

Container digests in [containers/MANIFEST.md](containers/MANIFEST.md) are the authoritative record; a tag can be retagged, a digest cannot.

The versions actually loaded at runtime are written per run to `results/<run>/pipeline_info/software_versions.yml`, generated from the processes themselves rather than transcribed by hand. The list below is from a reference run:

```
Nextflow           24.04+
STAR               2.7.11a
HISAT2             2.2.1
SAMtools           1.20
Salmon             1.10.2
Subread            2.0.6
tximport           1.32.0
DESeq2             1.44.0
edgeR              4.2.2
limma              3.60.6
clusterProfiler    4.12.6
org.Hs.eg.db       3.19.1
Seurat             5.0.0
SeuratObject       5.0.0
harmony            2.0.5
leidenalg          0.10.2
R                  4.4.1
Python             3.11
```

Reference annotation: **GRCh38.p14 / GENCODE v44**.

---

## Data availability and provenance

See [data/README.md](data/README.md). Summary:

| Dataset | Source | Accession | Access |
| --- | --- | --- | --- |
| Kazakh ESCC bulk | NCBI SRA | PRJNA608223 (22 tumour samples) | Open |
| ESCC single-cell atlas | GEO | GSE160269 | Open |
| TCGA-ESCA | GDC | TCGA-ESCA (92 squamous tumours; 13 paired tumour/normal, of which 2 squamous) | Open (harmonised counts) |
| Reference | GENCODE | v44 / GRCh38.p14 | Open |

**PRJNA608223 contains tumour samples only; it has no adjacent normal tissue.** This constrains what can be claimed from it, and is discussed in `data/README.md`.

**Anonymised local data:** none. Both clinical cohorts are publicly archived under their original anonymisation, and no patient-level identifiers are stored in this repository.

---

## Citation

If you use KazRNA-Pipe, please cite the archived software:

> Zhumadillayeva, A. KazRNA-Pipe: a reproducible Nextflow workflow for
> integrated bulk and single-cell transcriptomic profiling. Zenodo.
> https://doi.org/10.5281/zenodo.21832631

That is the concept DOI and always resolves to the most recent release. To cite the exact version used in an analysis, use the version DOI: v1.1.0 is https://doi.org/10.5281/zenodo.21832632.

A machine-readable citation is in [`CITATION.cff`](CITATION.cff), so GitHub's "Cite this repository" button produces a formatted reference.

---

## License

Released under the MIT License — see [LICENSE](LICENSE). Third-party tools invoked by the pipeline retain their own licenses; see [`docs/THIRD_PARTY_LICENSES.md`](docs/THIRD_PARTY_LICENSES.md). CIBERSORTx requires academic registration with the Stanford host (https://cibersortx.stanford.edu/); the pipeline detects whether a token has been provided and proceeds with BayesPrism and MuSiC if not.

---

## Acknowledgements

Funded by the Committee of Science of the Ministry of Science and Higher Education of the Republic of Kazakhstan, Grant No. BR28713313, *"Development of an intelligent computational model for automating the solution of problems using supercomputer resources, with applications in bioinformatics and quantum computing."*

Computational resources were provided by the KazHPC cluster at L.N. Gumilyov Eurasian National University.
