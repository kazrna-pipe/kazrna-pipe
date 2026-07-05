# Parameter reference

This page enumerates every `--param` accepted by `main.nf`. Default values are
those in `nextflow.config`; profile-specific overrides are noted.

## Input

| Parameter            | Type    | Default                              | Description                                                                                       |
|----------------------|---------|--------------------------------------|---------------------------------------------------------------------------------------------------|
| `--samplesheet`      | path    | (required)                           | CSV with `sample_id, subject_id, condition, fastq_1, fastq_2, strandedness, single_end, ...`      |
| `--sc_samplesheet`   | path    | (optional)                           | CSV for single-cell input: `sample_id, condition, matrix_h5, celltype_annotation, ...`            |
| `--workflows`        | string  | `bulk,scrnaseq,deconvolution`        | Comma-separated subset of `{bulk, scrnaseq, deconvolution, benchmark}`                            |

## References

| Parameter            | Type | Default                                                  | Description                            |
|----------------------|------|----------------------------------------------------------|----------------------------------------|
| `--star_index`       | path | `data/refs/star_index/`                                  | STAR index directory                    |
| `--hisat2_index`     | path | `data/refs/hisat2_index/GRCh38_v44`                      | HISAT2 index prefix                     |
| `--salmon_index`     | path | `data/refs/salmon_index/`                                | Salmon decoy-aware index directory      |
| `--gtf`              | path | `data/refs/gencode.v44.primary_assembly.annotation.gtf.gz` | GTF for featureCounts and STAR sjdb   |
| `--tx2gene`          | path | `data/refs/tx2gene.v44.tsv`                              | Transcript→gene map for tximport         |

## DE thresholds

| Parameter            | Type   | Default | Description                                  |
|----------------------|--------|---------|----------------------------------------------|
| `--fdr_threshold`    | double | 0.05    | Adjusted-p cutoff for DE and enrichment       |
| `--lfc_threshold`    | double | 1.0     | |log2 fold-change| cutoff for DE              |

## Single-cell

| Parameter            | Type   | Default      | Description                                                |
|----------------------|--------|--------------|------------------------------------------------------------|
| `--sc_n_hvgs`        | int    | 3000         | Number of highly variable features                          |
| `--sc_n_pcs`         | int    | 50           | Number of PCs retained                                      |
| `--sc_resolutions`   | string | `0.4,0.6,0.8,1.0,1.2,1.5` | Comma-separated Leiden resolutions                |
| `--run_gpu`          | bool   | true         | Enable the rapids-singlecell GPU back-end                  |
| `--marker_yaml`      | path   | `data/celltype_markers.yaml` | Cell-type marker panel                       |

## Deconvolution

| Parameter            | Type   | Default | Description                                                       |
|----------------------|--------|---------|-------------------------------------------------------------------|
| `--cibersortx_token` | string | null    | CIBERSORTx token (required to enable that back-end)               |
| `--cibersortx_email` | string | null    | Email associated with the token                                   |
| `--skip_deconvolution` | bool | false   | Skip the entire deconvolution sub-workflow                        |

## Benchmark

| Parameter            | Type   | Default              | Description                                                  |
|----------------------|--------|----------------------|--------------------------------------------------------------|
| `--benchmark_only`   | bool   | false                | Run only the benchmark sub-workflow                          |
| `--benchmark_nodes`  | string | `1,2,4,8,16`         | Comma-separated node counts                                  |
| `--benchmark_modes`  | string | `strong,weak`        | Comma-separated scaling modes                                |
| `--skip_benchmark`   | bool   | true (test profile)  | Skip the benchmark sub-workflow                              |

## Output

| Parameter   | Type | Default     | Description                |
|-------------|------|-------------|----------------------------|
| `--outdir`  | path | `results/`  | Top-level output directory |
