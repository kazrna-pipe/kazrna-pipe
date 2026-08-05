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
    make_option("--min_genes",  type = "integer",   default = 10,
                help = "Minimum mapped genes required to attempt enrichment"),
    make_option("--out_dir",    type = "character", default = ".")
)))

# ---- Mutable run record ----------------------------------------------------
rec <- new.env(parent = emptyenv())
rec$status   <- "ok"
rec$notes    <- character(0)
rec$n_sig    <- 0L
rec$n_mapped <- 0L
rec$n_kegg   <- 0L
rec$n_go     <- 0L

note <- function(msg) {
    rec$notes <- c(rec$notes, msg)
    message(msg)
}

# ---- Analysis --------------------------------------------------------------
run_enrichment <- function() {

    de <- read_tsv(opt$de_results, show_col_types = FALSE)

    gene_col <- intersect(colnames(de), c("gene_id", "gene", "ENSEMBL"))[1]
    lfc_col  <- intersect(colnames(de), c("log2FoldChange", "logFC"))[1]
    padj_col <- intersect(colnames(de), c("padj", "FDR", "adj.P.Val"))[1]
    stopifnot(!is.na(gene_col), !is.na(lfc_col), !is.na(padj_col))

    sig <- de |>
        dplyr::filter(!is.na(.data[[padj_col]]),
                      .data[[padj_col]] < opt$fdr,
                      abs(.data[[lfc_col]]) >= opt$lfc) |>
        dplyr::pull(!!gene_col) |>
        sub("\\..*$", "", x = _) 

    rec$n_sig <- length(sig)
    message("Significant genes: ", rec$n_sig)

    if (rec$n_sig < opt$min_genes) {
        rec$status <- "skipped"
        note(sprintf("Only %d significant genes (minimum %d); enrichment not attempted.",
                     rec$n_sig, opt$min_genes))
        return(invisible(NULL))
    }
    entrez <- tryCatch(
        suppressWarnings(
            bitr(sig, fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID),
        error = function(e) {
            note(paste("ID mapping failed:", conditionMessage(e)))
            character(0)
        })
    rec$n_mapped <- length(entrez)
    message("Genes mapped to Entrez: ", rec$n_mapped, " of ", rec$n_sig)

    if (rec$n_mapped < opt$min_genes) {
        rec$status <- "skipped"
        note(sprintf(paste("Only %d of %d genes mapped to Entrez (minimum %d);",
                           "enrichment not attempted. This is expected when",
                           "identifiers are not Ensembl gene IDs."),
                     rec$n_mapped, rec$n_sig, opt$min_genes))
        return(invisible(NULL))
    }

    # KEGG queries a remote API; a network hiccup must not fail the whole run.
    kegg <- tryCatch(
        enrichKEGG(gene = entrez, organism = opt$organism,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.25),
        error = function(e) {
            note(paste("KEGG enrichment failed:", conditionMessage(e)))
            NULL
        })
    if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
        rec$n_kegg <- nrow(as.data.frame(kegg))
        write_tsv(as.data.frame(kegg), file.path(opt$out_dir, "kegg_enrichment.tsv"))
        message("KEGG terms: ", rec$n_kegg)
    } else {
        note("No KEGG terms passed the significance cutoffs.")
    }

    go <- tryCatch(
        enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.25, readable = TRUE),
        error = function(e) {
            note(paste("GO enrichment failed:", conditionMessage(e)))
            NULL
        })
    if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
        rec$n_go <- nrow(as.data.frame(go))
        write_tsv(as.data.frame(go), file.path(opt$out_dir, "go_bp_enrichment.tsv"))
        message("GO BP terms: ", rec$n_go)
    } else {
        note("No GO BP terms passed the significance cutoffs.")
    }

    invisible(NULL)
}

# ---- Run, then always record provenance ------------------------------------
tryCatch(run_enrichment(), error = function(e) {
    rec$status <- "error"
    note(paste("Enrichment aborted:", conditionMessage(e)))
})

write_json(list(
    tool         = "clusterProfiler",
    version      = as.character(packageVersion("clusterProfiler")),
    orgdb        = as.character(packageVersion("org.Hs.eg.db")),
    de_results   = normalizePath(opt$de_results, mustWork = FALSE),
    organism     = opt$organism,
    fdr_cutoff   = opt$fdr,
    lfc_cutoff   = opt$lfc,
    n_sig        = rec$n_sig,
    n_mapped     = rec$n_mapped,
    n_kegg_terms = rec$n_kegg,
    n_go_terms   = rec$n_go,
    status       = rec$status,
    notes        = rec$notes,
    timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version    = R.version.string
), file.path(opt$out_dir, "enrichment_provenance.json"),
   pretty = TRUE, auto_unbox = TRUE)

message("Wrote enrichment_provenance.json (status: ", rec$status, ")")

# A genuine analysis error should still fail the task; a skip should not.
if (rec$status == "error") quit(status = 1)