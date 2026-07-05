# Test data

This directory contains a minimal smoke-test dataset that drives the `test`
Nextflow profile and the GitHub Actions CI workflow.

The data are aggressive subsamples designed to:

- exercise every process node in the workflow graph;
- complete on a 2-vCPU, 7 GB GitHub-hosted runner in **under 12 minutes**; and
- produce deterministic numerical outputs that the second-run reproducibility
  check can compare against.

## Contents

| File                                    | Source                         | Size    | Description                                            |
|-----------------------------------------|--------------------------------|---------|--------------------------------------------------------|
| `bulk_samplesheet.csv`                  | (this repo)                    | 1 KB    | 4 samples (2 tumour, 2 normal), paired-end             |
| `sc_samplesheet.csv`                    | (this repo)                    | 1 KB    | 2 samples, ~1000 cells each                            |
| `bulk_fastq/`                           | PRJNA608223 chr22 subset       | 720 MB  | 4 samples × paired-end, chr22-only reads               |
| `sc_h5/`                                | GSE160269 chr22 subset         | 110 MB  | Two 10x H5 files filtered to chr22 features            |
| `refs/star_index_chr22/`                | GENCODE v44 chr22              | 84 MB   | STAR index for chr22 only                              |
| `refs/hisat2_index_chr22/`              | GENCODE v44 chr22              | 22 MB   | HISAT2 index                                           |
| `refs/salmon_index_chr22/`              | GENCODE v44 chr22 transcripts  | 31 MB   | Salmon decoy-aware index                               |
| `refs/gencode.v44.chr22.gtf.gz`         | GENCODE v44                    | 1.8 MB  | chr22 annotation only                                  |
| `refs/tx2gene.chr22.tsv`                | derived                        | 250 KB  | tx → gene mapping                                      |

Subsampling was done with `seqtk sample -s 2025 <orig.fastq.gz> 250000` and
`samtools view -b -L chr22.bed <orig.bam>` followed by re-conversion to FASTQ.
The seed (`2025`) is the same one used pipeline-wide for stochastic stages.

## Verifying integrity

```bash
cd test_data && sha256sum -c SHA256SUMS
```

Any mismatch indicates the test data has been altered. CI will refuse to run
in that state. Regenerate the manifest with:

```bash
find . -type f ! -name 'SHA256SUMS' ! -name 'README.md' -print0 \
  | xargs -0 sha256sum > SHA256SUMS
```

## Why not ship the full PRJNA608223?

The full PRJNA608223 cohort is ~170 GB and is hosted by the European
Nucleotide Archive at <https://www.ebi.ac.uk/ena/browser/view/PRJNA608223>.
Redistributing it inside the test bundle would (i) violate ENA's terms of use,
(ii) blow past GitHub's 2 GB per-file LFS quota, and (iii) be slower to
download for CI than the chr22-only subset. The full cohort is fetched at
analysis time with `bin/sra_fetch.sh`.
