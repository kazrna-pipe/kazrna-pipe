# Test data

A synthetic fixture that drives the `test` Nextflow profile and the GitHub
Actions CI workflow.

Every file
below is generated from scratch by the scripts in this directory. Sample
labels such as `test_T1` and `test_N1` denote the two arbitrary groups the
planted fold changes are assigned between; they are not tumour or normal
tissue and are unrelated to PRJNA608223, which contains tumour samples only.

The fixture is designed to exercise every process node in the workflow graph,
complete on a 2-vCPU GitHub-hosted runner, and produce deterministic outputs
that the second-run reproducibility check can compare.

## Contents

| Path | Origin | Size | Description |
|---|---|---|---|
| `bulk_samplesheet.csv` | this repo | 1 KB | 4 samples, paired-end, two groups |
| `sc_samplesheet.csv` | this repo | 1 KB | 2 single-cell samples |
| `bulk_fastq/` | `make_reads.py` | 1.7 MB | 4 samples x paired-end, simulated |
| `sc_h5/` | `make_sc_h5.py` | 560 KB | Two 10x-format H5 matrices |
| `refs/mini.fa` | `make_reference.py` | - | One 304 kb contig, `chrTest` |
| `refs/gencode.v44.annotation.gtf` | `make_reference.py` | - | 200 synthetic genes |
| `refs/mini.transcripts.fa` | `make_transcriptome.py` | - | Transcript sequences |
| `refs/{star,hisat2,salmon}_index/` | `build_test_fixtures.sh` | 9 MB total | Built in the pinned containers |
| `refs/tx2gene.tsv` | `make_reference.py` | - | Transcript to gene mapping |

The GTF is named `gencode.v44.annotation.gtf` so that module paths match the
real run. Its contents are synthetic and share nothing with GENCODE v44.

200 genes is the smallest reference at which DESeq2 behaves as it does on real
data: `estimateDispersionsFit()` needs enough points to fit a mean-dispersion
trend. Forty genes carry a defined fold change between the two groups, so CI
can check that the planted signal is recovered.

## Regenerating

```bash
bash test_data/build_test_fixtures.sh
```

All generators use a fixed seed (`2025`), so the FASTA, GTF and FASTQ files
regenerate byte for byte. The aligner indices are built inside the pinned
containers; they will match provided those digests are unchanged.

## Verifying integrity

```bash
cd test_data && sha256sum -c SHA256SUMS
```

`SHA256SUMS` covers the files as shipped. Regenerate it with the final step of
`build_test_fixtures.sh`.

## Why the real cohort is not shipped

PRJNA608223 is roughly 170 GB and is hosted by the European Nucleotide Archive
at <https://www.ebi.ac.uk/ena/browser/view/PRJNA608223>. It is fetched at
analysis time with `bin/sra_fetch.sh`. Using a synthetic fixture for CI also
keeps the repository free of controlled or attributable sequence.
