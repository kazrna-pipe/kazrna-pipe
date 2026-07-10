# KazRNA-Pipe

**GPU-Accelerated, Reproducible Nextflow Workflow for Integrated Bulk and Single-Cell Transcriptomic Profiling**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/singularity-3.11-blue)](https://sylabs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

KazRNA-Pipe is a Nextflow DSL2 pipeline that processes bulk and single-cell RNA-seq data through a unified, containerised, GPU-aware workflow. It was developed for population-scale precision oncology research on Kazakhstan academic HPC infrastructure and was validated on a Kazakhstani esophageal squamous cell carcinoma (ESCC) cohort (PRJNA608223, n = 22) and on the largest public ESCC single-cell atlas (GSE160269, n = 208,659 cells).

---

## Quick links

- **Manuscript:** *KazRNA-Pipe: GPU-Accelerated, Reproducible Nextflow Workflow for Integrated Bulk and Single-Cell Transcriptomic Profiling.* BioMedInformatics 2026 (under review).
- **Reproducibility guide:** [REPRODUCIBILITY.md](REPRODUCIBILITY.md) - start here if you are a reviewer.
- **Data provenance:** [data/README.md](data/README.md) - exact accessions, download dates, query commands, SHA-256 manifests.
- **Tool versions and container digests:** [containers/MANIFEST.md](containers/MANIFEST.md) - fully pinned.
- **Issue tracker:** https://github.com/kazrna-pipe/kazrna-pipe/issues

---

## Contents

1. [What the pipeline does](#what-the-pipeline-does)
2. [Repository layout](#repository-layout)
3. [Installation](#installation)
4. [Running the test profile (12 min)](#running-the-test-profile-12-min)
5. [Reproducing the manuscript figures](#reproducing-the-manuscript-figures)
6. [Pinned versions](#pinned-versions)
7. [Data availability and provenance](#data-availability-and-provenance)
8. [Citation](#citation)
9. [License](#license)

---

## What the pipeline does

KazRNA-Pipe is organised into four modules, each independently runnable:

| Module                 | Inputs                       | Tools                                                            | Outputs                                       |
| ---------------------- | ---------------------------- | ---------------------------------------------------------------- | --------------------------------------------- |
| **Bulk RNA-seq**       | Paired-end FASTQ + samplesheet | FastQC, Trim Galore, STAR, HISAT2, Salmon, featureCounts, tximport | Count matrix, QC reports                      |
| **Differential expression** | Count matrix + design     | DESeq2, edgeR, limma-voom (parallel)                             | DEG tables, concordance summary, KEGG pathway results |
| **scRNA-seq**          | 10x filtered matrices        | Seurat v5 (CPU) **or** rapids-singlecell (GPU)                   | AnnData/.rds object, UMAP, cluster labels, marker genes |
| **Deconvolution**      | Bulk counts + sc reference   | CIBERSORTx, BayesPrism, MuSiC (parallel)                         | Cell-type proportions, per-compartment DEGs   |

Plus an HPC benchmarking module (strong/weak scaling, GPU utilisation, energy via RAPL + nvidia-smi).

A diagram is provided in [`docs/architecture.svg`](docs/architecture.svg).

---

## Repository layout

```
kazrna-pipe/
├── README.md                  This file
├── REPRODUCIBILITY.md         Step-by-step reproduction
├── CITATION.cff               Machine-readable citation
├── LICENSE                    MIT
├── nextflow.config            Top-level Nextflow config
├── main.nf                    Pipeline entry point
├── conf/                      Execution profiles (test, slurm, aws, docker)
├── modules/                   Per-tool Nextflow process definitions
├── workflows/                 Module-level workflows (bulk, scrnaseq, deconv, benchmark)
├── bin/                       Shell scripts: download data, prepare refs, run end-to-end
├── scripts/
│   ├── R/                     R analysis scripts (DESeq2, edgeR, limma-voom, Seurat, BayesPrism)
│   └── python/                Python scripts (rapids-singlecell, clustering agreement, plots)
├── containers/                Singularity definition files, MANIFEST.md with image digests
├── data/                      Sample sheets, accession lists, metadata snapshots, SHA-256 manifests
├── test_data/                 Subsampled FASTQ + sc matrices for CI (~1 GB)
├── docs/                      Architecture diagrams, parameter reference
├── assets/                    Logos, schematics
└── .github/workflows/ci.yml   Continuous integration (runs test profile)
```

---

## Installation

### Prerequisites

| Component   | Tested version | Notes                                        |
| ----------- | -------------- | -------------------------------------------- |
| Nextflow    | 23.10.0        | `curl -s https://get.nextflow.io \| bash`    |
| Java        | 17+ (OpenJDK)  | Required by Nextflow                         |
| Singularity | 3.11.0         | Or Docker 24+ via `-profile docker`          |
| Git LFS     | optional       | Only for cloning test data ≥ 100 MB         |
| CUDA        | 12.2           | GPU profile only; tested on A100, H100, L40  |

### Clone

```bash
git clone https://github.com/kazrna-pipe/kazrna-pipe
cd kazrna-pipe
```

---

## Running the test profile (12 min)

The test profile runs the **complete pipeline end-to-end** on a small subset of PRJNA608223 (2 tumour / 2 normal samples, downsampled to 2 million reads each) and a subset of GSE160269 (8 patients, ~25,000 cells). It is designed to complete in **≤ 12 minutes on a 4-core workstation with no GPU**.

```bash
nextflow run main.nf -profile test,singularity
```

Expected outputs land in `results/test/` and include:

- A multi-aligner concordance heatmap (`results/test/bulk/concordance.pdf`)
- A 3-way DEG Venn diagram (`results/test/de/venn.pdf`)
- A UMAP from the Seurat CPU path (`results/test/sc/umap_cpu.pdf`)
- A BayesPrism cell-type proportion barplot (`results/test/deconv/proportions.pdf`)
- A `pipeline_report.html` summary
- `software_versions.yml` listing every tool version actually used

This same profile runs automatically on every push via GitHub Actions (see [`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

---

## Reproducing the manuscript figures

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for a figure-by-figure walkthrough. In short:

```bash
# 1. Fetch the raw data (writes to data/raw/, ~ 280 GB)
bash bin/download_data.sh --target data/raw --verify-sha256

# 2. Prepare references (GRCh38.p14 + GENCODE v44, ~ 25 GB)
bash bin/prepare_references.sh --target data/refs

# 3. Run the full pipeline on the published cohort
nextflow run main.nf \
    -profile slurm,singularity \
    --input data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --outdir results/manuscript \
    --run_benchmark true

# 4. Generate the figures
python3 scripts/python/make_figures.py --indir results/manuscript --outdir figures/
python scripts/python/make_figures.py --indir results/manuscript --outdir figures/
```

Total wall-clock on the KazHPC reference configuration (16 nodes, 1× A100 per node): **~9 hours**.

---

## Pinned versions

All tool versions are pinned at two layers:

1. **Container digests** in [`containers/MANIFEST.md`](containers/MANIFEST.md) (e.g. `quay.io/biocontainers/star@sha256:f3a4...`). Tags can be retagged; digests cannot.
2. **Conda lockfiles** in [`containers/locks/`](containers/locks/) for users who prefer not to use containers.

Headline versions used in the manuscript:

```
Nextflow           23.10.0
Singularity        3.11.0
STAR               2.7.11a
HISAT2             2.2.1
Salmon             1.10.2
featureCounts      2.0.6
tximport           1.32.0
DESeq2             1.44.0
edgeR              4.2.1
limma              3.60.6
clusterProfiler    4.12.6
Seurat             5.1.0
Harmony            1.2.0
leidenalg          0.10.2
rapids-singlecell  0.10.10
BayesPrism         2.2
MuSiC              1.0.0
CIBERSORTx         v1.06 (Docker image fdc0c12a07c5)
```

Reference annotation: **GRCh38.p14 / GENCODE v44** (release date 2023-08-04).

---

## Data availability and provenance

See [data/README.md](data/README.md) for full details. Summary:

| Dataset       | Source                   | Accession                  | Access date  | Download command                                    |
| ------------- | ------------------------ | -------------------------- | ------------ | --------------------------------------------------- |
| Kazakh ESCC bulk | NCBI SRA              | PRJNA608223 (22 samples)   | 2025-11-03   | `bin/sra_fetch.sh PRJNA608223`                      |
| ESCC sc atlas | GEO                      | GSE160269                  | 2025-11-03   | `bin/geo_fetch.sh GSE160269`                        |
| TCGA-ESCA     | GDC via TCGAbiolinks     | TCGA-ESCA (n=80, SCC only) | 2025-11-04   | `Rscript scripts/R/fetch_tcga.R`                    |
| GTEx esophagus | GTEx portal v8          | dbGaP phs000424.v9.p2      | 2025-11-04   | `bash bin/fetch_gtex.sh`                            |
| GRCh38.p14    | Ensembl / GENCODE        | GENCODE v44                | 2025-11-04   | `bash bin/prepare_references.sh`                    |

Every download script writes a `<dataset>.sha256` manifest into `data/raw/<dataset>/`. Re-running the pipeline against a different manifest will produce a warning. The exact file lists, FASTQ-dump parameters, and request URLs are recorded in `data/README.md`.

**Anonymised local data:** None. Both clinical cohorts used in this work (PRJNA608223 and TCGA-ESCA) are publicly archived under their original anonymisation. No new patient-level identifiers are stored in this repository.

---

## Citation

If you use KazRNA-Pipe, please cite the manuscript:

```bibtex
@article{kazrna_pipe_2026,
  title   = {KazRNA-Pipe: GPU-Accelerated, Reproducible Nextflow Workflow for
             Integrated Bulk and Single-Cell Transcriptomic Profiling},
  author  = {TBA},
  journal = {BioMedInformatics},
  year    = {2026},
  doi     = {TBA}
}

```

A machine-readable citation is in [`CITATION.cff`](CITATION.cff). (To be Created)

---

## License

KazRNA-Pipe is released under the MIT License - see [LICENSE](LICENSE). Third-party tools invoked by the pipeline retain their own licenses; see [`docs/THIRD_PARTY_LICENSES.md`](docs/THIRD_PARTY_LICENSES.md) for the full list. CIBERSORTx in particular requires academic registration with the Stanford host (https://cibersortx.stanford.edu/); the pipeline auto-detects whether a token has been provided and falls back to BayesPrism + MuSiC if not.

---

## Acknowledgements

This work was funded by the Committee of Science of the Ministry of Science and Higher Education of the Republic of Kazakhstan, Grant No. BR28713313, *"Development of an intelligent computational model for automating the solution of problems using supercomputer resources, with applications in bioinformatics and quantum computing."*

Computational resources were provided by the KazHPC cluster at L.N. Gumilyov Eurasian National University.
