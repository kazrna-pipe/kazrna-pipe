#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# DESeq2 differential expression analysis for KazRNA-Pipe.
#
# Inputs
#   --counts   Tab-separated gene x sample raw count matrix (STAR featureCounts).
#   --meta     CSV with columns: sample_id, condition, age, sex.
#   --outdir   Output directory.
#
# Outputs
#   deseq2_results.tsv       Full DE table (all genes).
#   deseq2_significant.tsv   FDR < 0.05 & |log2FC| > 1 subset.
#   deseq2_normalized.tsv    Library-size-normalised counts (VST).
#   deseq2_session.txt       sessionInfo() output.
#
# This script is invoked by modules/deseq2.nf inside the container
#   quay.io/biocontainers/bioconductor-deseq2:1.44.0--r44hdfd78af_0
#
# Bioconductor DESeq2 1.44.0 - Love MI, Huber W, Anders S. Genome Biol 2014.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(tibble)
    library(readr)
    library(dplyr)
})

# ---- CLI -------------------------------------------------------------------
opt_list <- list(
    make_option("--counts",  type="character", help="Gene x sample count matrix (TSV)"),
    make_option("--meta",    type="character", help="Sample metadata CSV"),
    make_option("--outdir",  type="character", default="results/de/deseq2"),
    make_option("--fdr",     type="double",    default=0.05),
    make_option("--lfc",     type="double",    default=1.0),
    make_option("--contrast", type="character", default="condition,tumor,normal",
                help="Comma-separated contrast: factor,numerator,denominator"),
    make_option("--seed",    type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# ---- Load ------------------------------------------------------------------
message("Reading counts: ", opt$counts)
counts <- as.matrix(read.table(opt$counts, header=TRUE, row.names=1,
                               sep="\t", check.names=FALSE))
storage.mode(counts) <- "integer"

message("Reading metadata: ", opt$meta)
meta <- read_csv(opt$meta, show_col_types=FALSE) %>%
    as.data.frame() %>%
    { rownames(.) <- .$sample_id; . }

# Align columns of counts to rows of meta
common <- intersect(colnames(counts), rownames(meta))
if (length(common) < 4L) {
    stop("Fewer than 4 samples in common between counts and metadata; aborting.")
}
counts <- counts[, common, drop=FALSE]
meta   <- meta[common, , drop=FALSE]
meta$condition <- factor(meta$condition, levels=c("normal", "tumor"))

stopifnot(identical(colnames(counts), rownames(meta)))
message(sprintf("Analysing %d samples x %d genes.", ncol(counts), nrow(counts)))

# ---- DESeq2 ----------------------------------------------------------------
# Design includes age + sex as covariates per the manuscript Methods.
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = ~ age + sex + condition
)

# Independent filtering: drop genes with mean count < 10 (a conservative cut)
keep <- rowMeans(counts(dds)) >= 10
dds  <- dds[keep, ]
message(sprintf("Retained %d genes after expression filtering.", nrow(dds)))

dds <- DESeq(dds, parallel=FALSE)

# Extract the requested contrast
contrast_parts <- strsplit(opt$contrast, ",")[[1]]
if (length(contrast_parts) != 3L) stop("--contrast must be 'factor,num,denom'")
res <- results(dds, contrast=contrast_parts, alpha=opt$fdr)

# Shrink LFC for visualisation (apeglm). The unshrunken estimates are used
# for the DEG count reported in the manuscript Table 2.
res_shrunk <- lfcShrink(dds, contrast=contrast_parts, res=res, type="ashr")

# ---- Save ------------------------------------------------------------------
res_df <- res %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

res_shrunk_df <- res_shrunk %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    select(gene_id, log2FoldChange_shrunk = log2FoldChange, lfcSE_shrunk = lfcSE)

res_df <- left_join(res_df, res_shrunk_df, by="gene_id")

write_tsv(res_df, file.path(opt$outdir, "deseq2_results.tsv"))

sig <- res_df %>% filter(!is.na(padj), padj < opt$fdr, abs(log2FoldChange) > opt$lfc)
write_tsv(sig, file.path(opt$outdir, "deseq2_significant.tsv"))
message(sprintf("DEGs at FDR < %.3g and |log2FC| > %.2f: %d", opt$fdr, opt$lfc, nrow(sig)))

# VST-normalised counts for downstream visualisation
vst_counts <- assay(vst(dds, blind=FALSE))
vst_df <- as.data.frame(vst_counts) %>% rownames_to_column("gene_id")
write_tsv(vst_df, file.path(opt$outdir, "deseq2_normalized.tsv"))

# ---- Provenance ------------------------------------------------------------
sink(file.path(opt$outdir, "deseq2_session.txt"))
cat("KazRNA-Pipe DESeq2 module\n")
cat("=========================\n\n")
cat("Run datetime: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep="")
cat("Hostname:     ", Sys.info()["nodename"], "\n", sep="")
cat("Command:      ", paste(commandArgs(trailingOnly=FALSE), collapse=" "), "\n\n", sep="")
cat("Input SHA-256:\n")
cat(sprintf("  counts: %s\n", tools::md5sum(opt$counts)))
cat(sprintf("  meta:   %s\n\n", tools::md5sum(opt$meta)))
print(sessionInfo())
sink()

message("DESeq2 module complete. Output in: ", opt$outdir)
