# Reproducibility guide

This document describes how to reproduce the results of KazRNA-Pipe, and states
precisely what is and is not expected to reproduce exactly.

If you encounter a discrepancy, please open an issue at
https://github.com/kazrna-pipe/kazrna-pipe/issues with the contents of
`results/<run>/pipeline_info/software_versions.yml` and the relevant
`work/<hash>/.command.err`.

---

## 0. Start here: verify the pipeline runs

Before downloading any real data, confirm the pipeline runs on your machine.
This takes about 90 seconds and requires no GPU:

```bash
git clone https://github.com/kazrna-pipe/kazrna-pipe.git
cd kazrna-pipe
git checkout v1.1.0

nextflow run main.nf -profile test,docker --outdir results/check
```

Substitute `-profile test,singularity` if you use Apptainer. Both are exercised
by continuous integration on every push, so both are expected to work.

This runs on a **synthetic fixture** — 200 artificial genes, reads generated
*in silico*, no patient-derived sequence. Its purpose is to exercise every
process in the workflow quickly, and to demonstrate that the pipeline completes
end to end. See [test_data/README.md](test_data/README.md).

If this does not complete, stop here: there is no point attempting a full
reproduction until it does.

### What continuous integration already verifies

Every push to `main` runs:

1. **Interface check** — every module's options and declared outputs match the script it invokes (`tools/check_interfaces.py`).
2. **End-to-end run under Docker**, asserting the expected outputs exist, are non-empty, and that the fixture's planted differential signal is recovered.
3. **End-to-end run under Apptainer**, the same assertions.
4. **Determinism** — two independent runs compared numerically.

The result of the most recent run is shown by the badge on the repository front
page, and the full logs and output artifacts are retained for 30 days.

---

## 1. Environment

| Component | Required |
| --- | --- |
| Nextflow | ≥ 24.04 (the `resourceLimits` directive is used) |
| Java | 17+ |
| Docker or Apptainer | either |

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

No separate step is needed to obtain the software environment: every process
declares a container, and the three analysis images are pinned by SHA-256
digest in `conf/containers.config`. A run either obtains exactly those images
or fails.

---

## 2. Obtaining the input data

### 2.1 Bulk RNA-seq (PRJNA608223)

```bash
bash bin/sra_fetch.sh --accession PRJNA608223 --target data/raw/PRJNA608223
```

The 22 run accessions are listed in `data/samples_PRJNA608223.csv`, which is
also the samplesheet the pipeline reads. Each run is fetched with `prefetch`,
converted with `fasterq-dump --split-files`, and compressed.

**This cohort contains tumour samples only.** It includes no adjacent normal
tissue. See `data/README.md` for what that means for the analyses that can be
performed with it.

### 2.2 Single-cell atlas (GSE160269)

```bash
bash bin/geo_fetch.sh --accession GSE160269 --target data/raw/GSE160269
```

A snapshot of the GEO series matrix from the download date is kept under
`data/snapshots/`, so that the state of the record at the time of analysis can
be inspected even if GEO is later updated.

### 2.3 TCGA-ESCA

```bash
bash bin/fetch_tcga_esca.sh data/tcga_esca
```

This queries the GDC API directly; the filters are written into the script, so
the cohort definition is part of the repository rather than a web session that
cannot be replayed. It retrieves GDC-harmonised STAR gene counts (open access)
together with the sample metadata needed to distinguish tumour from adjacent
normal.

The underlying TCGA sequence is controlled-access (dbGaP phs000178). The TCGA
arm of the analysis therefore uses GDC's harmonised counts rather than passing
through this pipeline; this is stated so that the two arms are not assumed to
share a processing path.

### 2.4 Reference annotation

```bash
bash bin/prepare_references.sh --target data/refs
```

Downloads GRCh38.p14 and GENCODE v44 and builds the STAR, HISAT2 and Salmon
(decoy-aware) indices. **Building the STAR index requires roughly 32 GB of
RAM.**

---

## 3. Running the pipeline

### On a SLURM cluster

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    -resume \
    --input    data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --refs_dir data/refs \
    --outdir   results/manuscript
```

`-resume` re-enters the workflow without recomputing completed processes.

### On AWS

```bash
nextflow run main.nf -profile aws,docker \
    --input data/samples_PRJNA608223.csv \
    --sc_input data/samples_GSE160269.csv \
    --refs_dir s3://<your-bucket>/GRCh38.p14_GENCODE_v44 \
    --outdir   s3://<your-bucket>/results
```

### On a single workstation

```bash
nextflow run main.nf -profile local,docker \
    --max_memory 64.GB --max_cpus 16 \
    --input data/samples_PRJNA608223.csv \
    --outdir results/manuscript
```

A workstation with 64 GB will run the bulk module but cannot hold the full
GSE160269 single-cell analysis in memory. Use `--sc_subsample` to exercise the
single-cell path on a feasible subset.

---

## 4. What reproduces exactly, and what does not

This is the part worth reading carefully. Claiming more exactness than the
software delivers would be worse than claiming less.

### Bit-reproducible between runs

- Gene-level counts from STAR and HISAT2 (`featureCounts` output and the merged matrices)
- Salmon `NumReads`
- tximport gene-level counts
- All DESeq2, edgeR and limma-voom result tables
- Cross-method concordance output
- Single-cell cluster labels and UMAP coordinates

### Not bit-reproducible, by design

- **Salmon `EffectiveLength` and `TPM`.** `--gcBias` and `--seqBias` fit their models on a subsample of reads, so the fitted effective length varies slightly between runs, and TPM varies with it. `NumReads` is unaffected, and `NumReads` is what every downstream analysis in this pipeline consumes. Dropping the bias correction would make the check pass at the cost of the quantification, which is the wrong trade.
- **Alignment record order within BAM files.** With more than one thread, reads with equal coordinates are written in thread-completion order. The derived counts are identical.
- **Timestamps** in reports, logs, provenance JSON and serialised R objects.

### Verifying determinism yourself

```bash
nextflow run main.nf -profile test,docker --outdir results/run1 -work-dir w1
nextflow run main.nf -profile test,docker --outdir results/run2 -work-dir w2

docker run --rm -v "$PWD":/w -w /w \
  ghcr.io/kazrna-pipe/kazrna-py@sha256:17bc730681ca90d649f1bf25840f4672ded381747338085e0e6645a1c63f6832 \
  python3 scripts/python/diff_results.py results/run1 results/run2 --tolerance 1e-8
```

The exclusion list, and the reasoning behind each entry, is documented in the
script's own docstring. The check fails if it finds nothing to compare, so it
cannot pass vacuously on two empty runs.

---

## 5. Provenance recorded by each run

Every run writes into `results/<run>/pipeline_info/`:

| File | Contents |
| --- | --- |
| `software_versions.yml` | Every tool version actually invoked, collected from the processes themselves |
| `trace.txt` | Per-task resource usage, exit status and duration |
| `report.html` | Execution report |
| `timeline.html` | Task timeline |
| `dag.svg` | The workflow graph as executed |

Individual analysis steps additionally write `*_provenance.json` recording their
input paths and checksums, the parameters they were given, the resolved design
formula where applicable, and the package versions loaded at runtime.

---

## 6. Obtaining the exact software used

Each release is archived at Zenodo:

- All versions: https://doi.org/10.5281/zenodo.21832631
- v1.1.0: https://doi.org/10.5281/zenodo.21832632

To reproduce a specific analysis, use the version recorded in that run's
`pipeline_info/` rather than the current `main` branch, which will have moved
on.

---

## 7. Known caveats

- **CIBERSORTx** requires an academic token from https://cibersortx.stanford.edu/. Without `params.cibersortx_token`, deconvolution proceeds with BayesPrism and MuSiC.
- **The GPU single-cell path** falls back to Scanpy when `rapids_singlecell` cannot be imported, so the code path can be exercised without a GPU. Timings from that fallback are not GPU timings.
- **Enrichment analysis** requires identifiers that map to Entrez. On the synthetic fixture nothing maps, and the step records `status: skipped` in its provenance file rather than failing — this is the intended behaviour, and the same applies to any dataset whose identifiers are not Ensembl gene IDs.
- **GENCODE version.** Results are annotation-dependent. Pin to v44 to match the reference run.

---

## 8. Reporting problems

Open an issue including:

1. `nextflow -version` and your container engine version
2. `results/<run>/pipeline_info/software_versions.yml`
3. The failing task's `work/<hash>/.command.err` and `.command.sh`

`.command.sh` is usually the most informative: it shows the command line after
variable interpolation, so a missing parameter or an unexpected argument is
immediately visible.
