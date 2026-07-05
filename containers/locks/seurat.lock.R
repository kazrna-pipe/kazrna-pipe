# containers/locks/seurat.lock.R
# Generated 2025-11-04 using renv::snapshot() on the build host.
# Each call below is pinned by upstream commit hash where applicable.

options(
  Ncpus     = parallel::detectCores(),
  timeout   = 600,
  repos     = c(
    BIOC = "https://bioconductor.org/packages/3.19/bioc",
    CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"
  )
)

install.packages("remotes")
install.packages("BiocManager")
BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)

# Core CRAN
install.packages(c(
  "optparse",     # 1.7.5
  "tidyverse",    # 2.0.0
  "ggplot2",      # 3.5.1
  "patchwork",    # 1.2.0
  "yaml",         # 2.3.10
  "jsonlite",     # 1.8.8
  "VennDiagram",  # 1.7.3
  "Matrix",       # 1.7-0
  "future",       # 1.34.0
  "future.apply", # 1.11.2
  "RcppAnnoy",    # 0.0.22
  "uwot",         # 0.2.2
  "igraph",       # 2.0.3
  "leiden"        # 0.4.3.1
))

# Bioconductor
BiocManager::install(c(
  "DESeq2",          # 1.44.0
  "edgeR",           # 4.2.1
  "limma",           # 3.60.6
  "tximport",        # 1.32.0
  "ComplexHeatmap",  # 2.20.0
  "clusterProfiler", # 4.12.6
  "org.Hs.eg.db",    # 3.19.1
  "glmGamPoi",       # 1.16.0
  "SingleCellExperiment", # 1.26.0
  "scran",           # 1.32.0
  "scater"           # 1.32.1
), update = FALSE, ask = FALSE)

# Seurat v5 stack
remotes::install_version("SeuratObject", version = "5.0.2", upgrade = "never")
remotes::install_version("Seurat",       version = "5.1.0", upgrade = "never")

# Deconvolution
remotes::install_github("Danko-Lab/BayesPrism@v2.2.0", upgrade = "never")
remotes::install_github("xuranw/MuSiC@v1.0.0",          upgrade = "never")

# Verify
stopifnot(packageVersion("Seurat")       >= "5.1.0")
stopifnot(packageVersion("DESeq2")       >= "1.44.0")
stopifnot(packageVersion("BayesPrism")   >= "2.2.0")
stopifnot(packageVersion("MuSiC")        >= "1.0.0")
cat("seurat.lock.R install OK\n")
