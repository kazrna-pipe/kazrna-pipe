#!/usr/bin/env Rscript
# scripts/R/music_deconv.R
#
# MuSiC bulk deconvolution using the GSE160269 sc atlas as reference.

suppressPackageStartupMessages({
    library(optparse)
    library(MuSiC)
    library(SingleCellExperiment)
    library(Biobase)
    library(Seurat)
    library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--bulk",         type = "character"),
    make_option("--sc_reference", type = "character"),
    make_option("--threads",      type = "integer", default = 4),
    make_option("--out_prefix",   type = "character", default = "music")
)))

stopifnot(!is.null(opt$bulk), !is.null(opt$sc_reference))

bulk <- as.matrix(read.table(opt$bulk, header = TRUE, sep = "\t",
                             row.names = 1, check.names = FALSE))
sc   <- readRDS(opt$sc_reference)        # Seurat object

# Convert Seurat -> SingleCellExperiment for MuSiC
sce <- as.SingleCellExperiment(sc)
celltype_col <- if ("celltype" %in% colnames(colData(sce))) "celltype" else "seurat_clusters"
sample_col   <- if ("orig.ident" %in% colnames(colData(sce))) "orig.ident" else "sample_id"

set.seed(2025)
res <- music_prop(bulk.mtx        = bulk,
                  sc.sce          = sce,
                  clusters        = celltype_col,
                  samples         = sample_col,
                  select.ct       = unique(colData(sce)[[celltype_col]]),
                  verbose         = TRUE)

write.table(round(res$Est.prop.weighted, 6),
            paste0(opt$out_prefix, "_proportions.tsv"),
            quote = FALSE, sep = "\t", col.names = NA)

prov <- list(
    tool      = "MuSiC",
    version   = as.character(packageVersion("MuSiC")),
    seed      = 2025,
    n_bulk    = ncol(bulk),
    n_sc      = ncol(sc),
    celltypes = unique(colData(sce)[[celltype_col]]),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(prov, paste0(opt$out_prefix, "_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)

writeLines(capture.output(sessionInfo()),
           paste0(opt$out_prefix, "_session_info.txt"))
