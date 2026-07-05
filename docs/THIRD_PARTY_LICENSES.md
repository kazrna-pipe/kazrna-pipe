# Third-party licences

KazRNA-Pipe itself is released under the MIT licence (see `LICENSE`). The
tools and datasets it composes have their own licences, summarised below.
This file is informational; it does not modify upstream licence terms.

## Software

| Tool                  | Version    | Licence            | Upstream                                                                  |
|-----------------------|------------|--------------------|----------------------------------------------------------------------------|
| Nextflow              | 24.04.4    | Apache-2.0         | <https://github.com/nextflow-io/nextflow>                                  |
| STAR                  | 2.7.11a    | MIT                | <https://github.com/alexdobin/STAR>                                        |
| HISAT2                | 2.2.1      | GPL-3.0            | <https://github.com/DaehwanKimLab/hisat2>                                  |
| Salmon                | 1.10.2     | GPL-3.0            | <https://github.com/COMBINE-lab/salmon>                                    |
| Subread (featureCounts) | 2.0.6    | GPL-3.0            | <https://subread.sourceforge.net>                                          |
| FastQC                | 0.12.1     | GPL-3.0            | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>               |
| Trim Galore           | 0.6.10     | GPL-3.0            | <https://github.com/FelixKrueger/TrimGalore>                               |
| Cutadapt              | 4.9        | MIT                | <https://github.com/marcelm/cutadapt>                                      |
| Samtools              | 1.20       | MIT                | <https://www.htslib.org>                                                   |
| MultiQC               | 1.24.1     | GPL-3.0            | <https://multiqc.info>                                                     |
| R                     | 4.4.1      | GPL-2.0 / GPL-3.0  | <https://www.r-project.org>                                                |
| DESeq2                | 1.44.0     | LGPL-3.0           | <https://bioconductor.org/packages/DESeq2/>                                |
| edgeR                 | 4.2.1      | GPL-2.0+           | <https://bioconductor.org/packages/edgeR/>                                 |
| limma                 | 3.60.6     | GPL-2.0+           | <https://bioconductor.org/packages/limma/>                                 |
| tximport              | 1.32.0     | GPL-2.0+           | <https://bioconductor.org/packages/tximport/>                              |
| clusterProfiler       | 4.12.6     | Artistic-2.0       | <https://bioconductor.org/packages/clusterProfiler/>                       |
| Seurat                | 5.1.0      | MIT                | <https://github.com/satijalab/seurat>                                      |
| BayesPrism            | 2.2.0      | GPL-3.0            | <https://github.com/Danko-Lab/BayesPrism>                                  |
| MuSiC                 | 1.0.0      | GPL-3.0            | <https://github.com/xuranw/MuSiC>                                          |
| rapids-singlecell     | 0.10.10    | MIT                | <https://github.com/scverse/rapids_singlecell>                             |
| Scanpy                | 1.10.2     | BSD-3-Clause       | <https://github.com/scverse/scanpy>                                        |
| AnnData               | 0.10.8     | BSD-3-Clause       | <https://github.com/scverse/anndata>                                       |
| RAPIDS (cudf, cuml)   | 24.06      | Apache-2.0         | <https://github.com/rapidsai>                                              |
| CIBERSORTx            | (image)    | Restricted / token | <https://cibersortx.stanford.edu>  *(not redistributed; see MANIFEST.md)*  |

## Data

| Dataset       | Source                  | Terms of use                                                                                             |
|---------------|-------------------------|----------------------------------------------------------------------------------------------------------|
| PRJNA608223   | NCBI SRA / ENA          | Open. Cite Sharip et al. 2024, Front. Genet. 15:1249751, doi:10.3389/fgene.2024.1249751.                                       |
| GSE160269     | NCBI GEO                | Open. Cite Zhang et al. 2021, *Nat Commun*, 12:5291.                                                      |
| TCGA-ESCA     | NCI GDC                 | Open (de-identified RNA-seq counts). Cite the TCGA Research Network and the GDC Data Release used.        |
| GTEx v8       | GTEx Portal             | Open (summary data). Cite the GTEx Consortium, *Science* 2020.                                            |
| GENCODE v44   | EMBL-EBI                | Open (CC0).                                                                                               |
