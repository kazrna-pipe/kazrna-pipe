#!/usr/bin/env Rscript
# scripts/R/tximport_aggregate.R
#
# Aggregate Salmon transcript-level quant.sf files to gene-level counts,
# TPM, and effective lengths using tximport.
#
# Usage:
#   Rscript tximport_aggregate.R --quant_root . --tx2gene tx2gene.tsv \
#                                --samples samples.csv --out_prefix salmon.gene

suppressPackageStartupMessages({
    library(optparse)
    library(tximport)
    library(readr)
    library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--quant_root", type = "character", help = "Root dir containing <sample>/quant.sf"),
    make_option("--tx2gene",    type = "character", help = "TSV: transcript_id, gene_id, gene_name"),
    make_option("--samples",    type = "character", help = "Sample sheet CSV with sample_id column"),
    make_option("--out_prefix", type = "character", default = "salmon.gene")
)))

stopifnot(!is.null(opt$quant_root), !is.null(opt$tx2gene), !is.null(opt$samples))

samples <- read_csv(opt$samples,  show_col_types = FALSE)
tx2gene <- read_tsv(opt$tx2gene,  show_col_types = FALSE,
                    col_names = c("tx_id", "gene_id", "gene_name"))

files <- file.path(opt$quant_root, samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing quant.sf files: ", paste(missing, collapse = ", "))

txi <- tximport(files,
                type     = "salmon",
                tx2gene  = tx2gene[, c("tx_id", "gene_id")],
                countsFromAbundance = "no")

write.table(round(txi$counts, 4),
            paste0(opt$out_prefix, "_counts.tsv"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(round(txi$abundance, 6),
            paste0(opt$out_prefix, "_tpm.tsv"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(round(txi$length, 2),
            paste0(opt$out_prefix, "_lengths.tsv"),
            quote = FALSE, sep = "\t", col.names = NA)

prov <- list(
    tool        = "tximport",
    version     = as.character(packageVersion("tximport")),
    r_version   = paste(R.version$major, R.version$minor, sep = "."),
    n_samples   = nrow(samples),
    n_genes     = nrow(txi$counts),
    quant_root  = normalizePath(opt$quant_root),
    tx2gene     = normalizePath(opt$tx2gene),
    samples     = normalizePath(opt$samples),
    timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(prov, paste0(opt$out_prefix, "_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)
