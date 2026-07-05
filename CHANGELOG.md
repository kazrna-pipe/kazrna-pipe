# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-04

First public, FAIR-compliant release accompanying the manuscript submission to
*BioMedInformatics* (2026).

### Added
- Nextflow DSL2 workflow covering bulk RNA-seq, single-cell RNA-seq (CPU and
  GPU), and bulk deconvolution, packaged as four Singularity images pinned by
  OCI digest.
- Three-aligner / three-method DE concordance analysis (STAR/HISAT2/Salmon ×
  DESeq2/edgeR/limma-voom).
- GPU-accelerated single-cell pipeline via rapids-singlecell 0.10.10 with a
  CPU-equivalence check against Seurat v5.1.0 (ARI ≥ 0.94 at matched
  resolutions).
- BayesPrism + MuSiC + optional CIBERSORTx deconvolution with cross-method
  correlation check (Fig 5B).
- Cell-type-resolved DE via pseudo-bulk DESeq2 (Fig 5D).
- HPC benchmarking sub-workflow with RAPL CPU and nvidia-smi GPU energy
  sampling (Fig 4).
- Reproducibility infrastructure:
  - Frozen SHA-256 manifests for every input dataset and reference file.
  - Frozen GDC query JSON for the TCGA-ESCA cohort.
  - Conda lockfiles and an R `install.packages()` script pinned by version.
  - GitHub Actions CI that runs the full test profile and verifies a
    bit-/numerically-equivalent second run.

### Notes
- CIBERSORTx is **not** redistributed (its licence forbids this). Users wishing
  to enable that arm of the deconvolution comparison must obtain a token from
  the CIBERSORTx project; the pipeline runs without it.
- The frozen manifests in `data/snapshots/` were captured 2025-11-03 to
  2025-11-04. If upstream archives re-version the records, `bin/sra_fetch.sh`
  and `bin/geo_fetch.sh` will halt and report the discrepancy.
