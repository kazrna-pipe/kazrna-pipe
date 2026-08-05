# Datasets

None of the primary sequencing data is redistributed in this repository. This
file records what was used, where it comes from, and how to obtain it.

## Primary cohort — PRJNA608223

| | |
|---|---|
| Accession | [PRJNA608223](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA608223) |
| Samples used | 22 |
| Composition | **Tumour only. This cohort contains no adjacent-normal samples.** |
| Manifest | `data/samples_PRJNA608223.csv` |

The absence of matched normals in this cohort is a design constraint of the
study, not an oversight, and it determines what can and cannot be claimed from
it. Any tumour-versus-normal contrast involving these samples necessarily draws
its non-tumour samples from a different study, which confounds tumour status
with study of origin. See the Limitations section of the manuscript.

Download:

```bash
# requires sra-tools
while IFS=, read -r run rest; do
    [ "$run" = "run_accession" ] && continue
    fasterq-dump --split-files --outdir data/fastq "$run"
done < data/samples_PRJNA608223.csv
```

## Single-cell reference — GSE160269

| | |
|---|---|
| Accession | [GSE160269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160269) |
| Role | Cell-type reference for deconvolution and for benchmarking estimated proportions |

## Comparison cohort — TCGA-ESCA

Obtained through the GDC Data Portal. Used where a within-study
tumour/adjacent-normal contrast is required, since matched pairs exist there.

## Reference genome and annotation

| | |
|---|---|
| Genome | GRCh38.p14 |
| Annotation | GENCODE v44 |
| Indices | Built locally; see `bin/build_references.sh` |

## Checksums

After downloading, verify against the recorded checksums:

```bash
sha256sum -c data/PRJNA608223.sha256
```

If that file is absent from your clone, generate it once from a verified
download and commit it, so subsequent users can check their own copies.

## What is in this directory

| File | Contents |
|---|---|
| `samples_PRJNA608223.csv` | Run accessions and sample metadata for the primary cohort |
| `celltype_markers.yaml` | Marker gene sets used for single-cell cluster annotation |
| `benchmark_samplesheet.csv` | Input manifest for the scaling benchmark |
| `PRJNA608223.sha256` | Checksums for the downloaded FASTQ files |

For the synthetic fixture used by CI, see `test_data/README.md`.
