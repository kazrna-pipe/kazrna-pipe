#!/usr/bin/env Rscript
# scripts/R/enrichment_analysis.R
#
# KEGG and GO Biological Process over-representation on DE gene lists.

suppressPackageStartupMessages({
    library(optparse)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(readr)
    library(dplyr)
    library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--de_results", type = "character"),
    make_option("--organism",   type = "character", default = "hsa"),
    make_option("--fdr",        type = "double",    default = 0.05),
    make_option("--lfc",        type = "double",    default = 1.0),
    make_option("--out_dir",    type = "character", default = ".")
)))

de <- read_tsv(opt$de_results, show_col_types = FALSE)
# Be flexible about column names emitted by DESeq2/edgeR/limma wrappers.
gene_col <- intersect(colnames(de), c("gene_id", "gene", "ENSEMBL"))[1]
lfc_col  <- intersect(colnames(de),
                      c("log2FoldChange", "logFC"))[1]
padj_col <- intersect(colnames(de),
                      c("padj", "FDR", "adj.P.Val"))[1]
stopifnot(!is.na(gene_col), !is.na(lfc_col), !is.na(padj_col))

sig <- de |>
    dplyr::filter(.data[[padj_col]] < opt$fdr,
                  abs(.data[[lfc_col]]) >= opt$lfc) |>
    dplyr::pull(!!gene_col) |>
    sub("\\..*$", "", x = _)            # strip Ensembl version suffix

if (length(sig) < 10) {
    message("Too few significant genes (", length(sig), ") for enrichment")
    quit(status = 0)
}

entrez <- bitr(sig, fromType = "ENSEMBL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID

if (length(entrez) >= 10) {
    kegg <- enrichKEGG(gene = entrez, organism = opt$organism,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.25)
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
        write_tsv(as.data.frame(kegg),
                  file.path(opt$out_dir, "kegg_enrichment.tsv"))
    }

    go <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db,
                   ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.25,
                   readable = TRUE)
    if (!is.null(go) && nrow(go@result) > 0) {
        write_tsv(as.data.frame(go),
                  file.path(opt$out_dir, "go_bp_enrichment.tsv"))
    }
}

prov <- list(
    tool        = "clusterProfiler",
    version     = as.character(packageVersion("clusterProfiler")),
    orgdb       = as.character(packageVersion("org.Hs.eg.db")),
    de_results  = normalizePath(opt$de_results),
    n_sig       = length(sig),
    n_mapped    = length(entrez),
    fdr_cutoff  = opt$fdr,
    lfc_cutoff  = opt$lfc,
    timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(prov, file.path(opt$out_dir, "enrichment_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)
