# Reproducibility guide

This document tells you exactly how to make every quantitative result, on three independently-tested environments (a SLURM cluster, an AWS EC2 instance, and a local Ubuntu workstation). If you encounter a discrepancy, please open an issue at https://github.com/kazrna-pipe/kazrna-pipe/issues with the output of `nextflow log` and the contents of `results/<run>/pipeline_info/software_versions.yml`.

---

## 0. What a successful reproduction looks like

After running every step in this guide, you will have:

| Manuscript figure | Output file (under `results/manuscript/`)                          | Approximate runtime on KazHPC |
|-------------------|--------------------------------------------------------------------|-------------------------------|
| Figure 2A         | `bulk/concordance/spearman_matrix.pdf` + `.tsv`                     | 6 h (all 22 samples × 3 aligners) |
| Figure 2B         | `de/venn_DESeq2_edgeR_limma.pdf` + per-method TSVs                  | 12 min                        |
| Figure 2C         | `de/tcga_vs_kazakh_jaccard.pdf` + `tcga_kazakh_shared_genes.tsv`    | 1.5 h (TCGA download + DESeq2) |
| Figure 2D         | `bulk/resource_per_sample.pdf` + `nextflow_trace.tsv`               | Generated alongside Fig 2A    |
| Figure 3A         | `sc/timing/per_step_walltime.pdf` + `per_step_walltime.tsv`         | 2.5 h (CPU+GPU full atlas)    |
| Figure 3B         | `sc/timing/scaling_vs_n_cells.pdf`                                  | 4 h (subsampled benchmark)    |
| Figure 3C         | `sc/agreement/leiden_ari_nmi_asw.pdf`                               | 30 min                        |
| Figure 3D         | `sc/umaps/umap_cpu.pdf`, `umap_gpu.pdf`                             | Generated alongside Fig 3A    |
| Figure 4A         | `deconv/method_agreement.pdf` + `cross_method_correlations.tsv`     | 8 h (CIBERSORTx is the slow one) |
| Figure 4B         | `deconv/sample_proportions.pdf` + `bayesprism_proportions.tsv`      | Generated alongside Fig 4A    |
| Figure 4C         | `deconv/celltype_specific_degs.pdf` + per-compartment DEG TSVs      | 1 h                           |
| Figure 4D         | `deconv/effect_size_amplification.pdf`                              | Generated alongside Fig 4C    |
| Figure 5          | `benchmark/strong_scaling.pdf`, `weak_scaling.pdf`                  | 18 h (full scaling sweep)     |

Total wall-clock on the reference KazHPC configuration: **~9 hours** for everything except the strong-scaling sweep, which is run separately and adds another 18 h.

---

## 1. Environment setup

### 1.1 Install Nextflow and a container runtime

```bash
# Nextflow
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/
nextflow -version    # must print 23.10.0 or higher

# Singularity (preferred) or Docker
singularity --version  # 3.11.0 used in the manuscript
```

### 1.2 Clone and check out the manuscript tag

```bash
git clone https://github.com/kazrna-pipe/kazrna-pipe.git
cd kazrna-pipe
git checkout v1.0.0
git rev-parse HEAD
```

### 1.3 Validate the test profile first

Before downloading 280 GB of data, confirm that the pipeline runs at all on your machine:

```bash
nextflow run main.nf -profile test,singularity
```

This must complete in ≤ 15 minutes and produce `results/test/pipeline_report.html`. If it does not, stop here and open an issue - there is no point trying the full reproduction until the test profile passes.

---

## 2. Fetching the input data

### 2.1 Bulk RNA-seq (PRJNA608223)

```bash
bash bin/sra_fetch.sh \
    --accession PRJNA608223 \
    --target data/raw/PRJNA608223 \
    --verify-sha256
```

This script:

1. Resolves PRJNA608223 to its 44 run accessions (`SRR11242833`…`SRR11242876`) via the NCBI EUtils API. The exact run list is also stored in `data/PRJNA608223_runs.txt`.
2. Downloads each run with `prefetch --max-size 100g`, then converts to paired-end FASTQ with `fasterq-dump --split-files --threads 8`.
3. Compresses with `pigz -p 8`.
4. Writes `data/raw/PRJNA608223/sha256.manifest`.
5. Compares against `data/PRJNA608223.sha256` (the manifest at the time of the manuscript, **frozen on 2025-11-03**). If any hash differs, the script exits with a non-zero code and a per-file diff.

Expected size on disk: 142 GB.

### 2.2 Single-cell atlas (GSE160269)

```bash
bash bin/geo_fetch.sh \
    --accession GSE160269 \
    --target data/raw/GSE160269 \
    --verify-sha256
```

This downloads the filtered count matrices (`*_filtered_feature_bc_matrix.tar.gz` for each of the 60 patients) directly from GEO FTP. The mirrored snapshot of the GEO `series_matrix.txt` from the download date is in `data/snapshots/GSE160269_series_matrix_2025-11-03.txt.gz` so that you can verify what GEO looked like at the time of the manuscript even if it is later updated.

Expected size on disk: 88 GB.

### 2.3 TCGA-ESCA

```bash
Rscript scripts/R/fetch_tcga.R \
    --project TCGA-ESCA \
    --histology "Squamous Cell Neoplasms" \
    --workflow "STAR - Counts" \
    --outdir data/raw/TCGA-ESCA
```

The exact `TCGAbiolinks::GDCquery` call is recorded in the script. The query was issued on **2025-11-04**; the GDC manifest returned at that time is mirrored at `data/snapshots/TCGA-ESCA_gdc_manifest_2025-11-04.tsv`. Re-running the query later may return more samples (the GDC is updated periodically); restrict your run to the manifest snapshot to reproduce the manuscript exactly.

### 2.4 Reference annotation

```bash
bash bin/prepare_references.sh --target data/refs
```

Downloads GRCh38.p14 primary assembly and GENCODE v44 annotation from EBI, builds STAR (sjdbOverhang 100), HISAT2, and Salmon (decoy-aware) indices. All download URLs and SHA-256 manifests are listed in `bin/prepare_references.sh` and the resulting index files are hashed into `data/refs/sha256.manifest`.

---

## 3. End-to-end pipeline run

### 3.1 On a SLURM cluster

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    -resume \
    --input          data/samples_PRJNA608223.csv \
    --sc_input       data/samples_GSE160269.csv \
    --tcga_input     data/samples_TCGA_ESCA.csv \
    --refs_dir       data/refs \
    --outdir         results/manuscript \
    --run_benchmark  true \
    -with-trace      results/manuscript/trace.tsv \
    -with-report     results/manuscript/report.html \
    -with-timeline   results/manuscript/timeline.html
```

`-resume` lets you re-enter the workflow without recomputing already-completed processes; this matters because the full run takes ~9 hours.

### 3.2 On AWS EC2

```bash
nextflow run main.nf \
    -profile aws,docker \
    --input data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --refs_dir s3://kazrna-pipe-refs/GRCh38.p14_GENCODE_v44 \
    --outdir s3://kazrna-pipe-results/manuscript
```

The `aws` profile is documented in `conf/aws.config`. A `c5.24xlarge` instance with one attached G5 GPU completed the full pipeline in 11h 14m at a total cost of $94.20 (us-east-2, on-demand pricing, November 2025).

### 3.3 On a local workstation with one GPU

```bash
nextflow run main.nf \
    -profile local,singularity \
    --max_memory 64.GB \
    --max_cpus 16 \
    --input data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --outdir results/manuscript
```

A workstation with 64 GB RAM and a single A100 will run the bulk module fine but cannot fit the full GSE160269 single-cell run in memory. Use `--sc_subsample 60000` to verify the GPU code path on a feasible subset.

---

## 4. Generating the figures

```bash
Rscript scripts/R/figure_generation.R \
    --indir results/manuscript \
    --outdir figures/

python scripts/python/make_figures.py \
    --indir results/manuscript \
    --outdir figures/
```

Each figure script writes a `<figure>.provenance.json` alongside the PDF that records:

- The input file paths and their SHA-256 hashes.
- The script version (`git describe`).
- The R / Python / library versions actually loaded at runtime.
- The exact `argparse` arguments.

If you regenerate a figure and the JSON differs from `figures/expected/<figure>.provenance.json`, that is your first diagnostic.

---

## 5. Verifying the results

Numerical reproducibility is checked against `results/expected/` (committed to the repository):

```bash
bash bin/verify_results.sh --results results/manuscript --expected results/expected
```

The script compares:

1. **Bit-identical** for: STAR count matrices, HISAT2 count matrices, Salmon quant.sf, DESeq2 result tables, edgeR result tables, limma-voom result tables.
2. **Numerically close** (relative error < 1e-6) for: BayesPrism cell-type proportions (Markov chain stochastic).
3. **ARI ≥ 0.99** between rapids-singlecell cluster labels across runs (cuGraph parallel Leiden has non-deterministic node visiting order). This was documented in Section 4.4 of the manuscript.

Any deviation outside these tolerances is a reproducibility failure and should be reported.

---

## 6. Known caveats and platform-specific notes

- **CIBERSORTx** requires an academic token from https://cibersortx.stanford.edu/. The pipeline auto-detects whether `params.cibersortx_token` is set; if not, deconvolution proceeds with BayesPrism and MuSiC only and Figure 4A is generated as a 2-method comparison.
- **rapids-singlecell** versions ≤ 0.10.9 had a known bug in PCA initialisation on H100. The pipeline pins 0.10.10 specifically to avoid this. If you build a custom container, do not downgrade.
- **GPU clustering** uses cuGraph's parallel Leiden, whose community-merging order is non-deterministic across runs. Cluster *IDs* will differ between runs but cluster *boundaries* are preserved (ARI > 0.99). All downstream analyses use ARI / NMI / ASW, which are invariant to label permutation.
- **STAR index loading** is the dominant I/O step at scale. The strong-scaling plateau beyond 16 nodes in Figure 5 is reproducible across all three test environments.
- **GENCODE v44 vs v45**: re-running with v45 changes 12 DEGs in the core 2,104 set. We recommend pinning to v44 for an exact reproduction of the manuscript.

---

## 7. Reporting problems

Open an issue at https://github.com/kazrna-pipe/kazrna-pipe/issues including:

1. Output of `nextflow -version`, `singularity --version`, `uname -a`.
2. The contents of `results/<run>/pipeline_info/software_versions.yml`.
3. The Nextflow `.nextflow.log` (or the relevant tail).
4. Which step failed and the contents of `work/<hash>/.command.log`.

We commit to responding to reproducibility issues within five working days.
