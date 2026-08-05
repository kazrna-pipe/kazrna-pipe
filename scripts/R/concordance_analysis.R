#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Cross-method concordance analysis for KazRNA-Pipe.
#
# Generates Figure 2A (Spearman correlation matrix across STAR/HISAT2/Salmon)
# and Figure 2B (3-way Venn of DEGs across DESeq2/edgeR/limma-voom).
#
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(readr); library(dplyr); library(tibble); library(tidyr)
    library(ggplot2); library(VennDiagram); library(pheatmap)
    library(jsonlite)
})

opt_list <- list(
    make_option("--star_counts",   type="character"),
    make_option("--hisat2_counts", type="character"),
    make_option("--salmon_counts", type="character"),
    make_option("--deseq2_degs",   type="character"),
    make_option("--edger_degs",    type="character"),
    make_option("--limma_degs",    type="character"),
    make_option("--outdir",        type="character", default="results/bulk/concordance"),
    make_option("--fdr",           type="double", default=0.05),
    make_option("--lfc",           type="double", default=1.0)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

read_counts <- function(path, label) {
    m <- as.matrix(read.table(path, header=TRUE, row.names=1, sep="\t",
                              check.names=FALSE))
    log1p(m)
}

# ---- Figure 2A: Spearman correlation across the three aligners -------------
sm <- read_counts(opt$star_counts,   "star")
hm <- read_counts(opt$hisat2_counts, "hisat2")
qm <- read_counts(opt$salmon_counts, "salmon")

# Align genes and samples
common_genes   <- Reduce(intersect, list(rownames(sm), rownames(hm), rownames(qm)))
common_samples <- Reduce(intersect, list(colnames(sm), colnames(hm), colnames(qm)))
sm <- sm[common_genes, common_samples, drop=FALSE]
hm <- hm[common_genes, common_samples, drop=FALSE]
qm <- qm[common_genes, common_samples, drop=FALSE]

# Per-sample Spearman across aligners, averaged
sample_rho <- function(a, b) {
    sapply(seq_len(ncol(a)), function(i) cor(a[, i], b[, i], method="spearman"))
}
rho_star_hisat2 <- mean(sample_rho(sm, hm))
rho_star_salmon <- mean(sample_rho(sm, qm))
rho_hisat2_salmon <- mean(sample_rho(hm, qm))

cor_mat <- matrix(
    c(1.0,             rho_star_hisat2,  rho_star_salmon,
      rho_star_hisat2, 1.0,              rho_hisat2_salmon,
      rho_star_salmon, rho_hisat2_salmon, 1.0),
    nrow=3, dimnames=list(c("STAR","HISAT2","Salmon"), c("STAR","HISAT2","Salmon"))
)
write.table(cor_mat, file.path(opt$outdir, "spearman_matrix.tsv"),
            sep="\t", quote=FALSE)

pdf(file.path(opt$outdir, "spearman_matrix.pdf"), width=5, height=4.5)
pheatmap(cor_mat, display_numbers=TRUE, fontsize_number=11,
         number_format="%.3f",
         color=colorRampPalette(c("#ffffff","#2e75b6"))(50),
         breaks=seq(0.9, 1, length.out=51),
         cluster_rows=FALSE, cluster_cols=FALSE,
         main="Cross-aligner gene-level Spearman ρ")
dev.off()
message(sprintf("Figure 2A: STAR-HISAT2=%.3f STAR-Salmon=%.3f HISAT2-Salmon=%.3f",
                rho_star_hisat2, rho_star_salmon, rho_hisat2_salmon))

# ---- Figure 2B: 3-way Venn across DE methods -------------------------------

read_sig <- function(path, gene_col="gene_id") {
    df <- read_tsv(path, show_col_types=FALSE)

    if (!gene_col %in% colnames(df))
        stop("No '", gene_col, "' column in ", path)

    padj_col <- intersect(colnames(df), c("padj", "FDR", "adj.P.Val"))[1]
    lfc_col  <- intersect(colnames(df), c("log2FoldChange", "logFC"))[1]
    if (is.na(padj_col) || is.na(lfc_col))
        stop("Could not identify adjusted-p and log-fold-change columns in ", path,
             ". Columns present: ", paste(colnames(df), collapse=", "))

    n_all <- nrow(df)
    df <- df %>% filter(!is.na(.data[[padj_col]]),
                        .data[[padj_col]] < opt$fdr,
                        abs(.data[[lfc_col]]) > opt$lfc)
    message(sprintf("  %s: %d of %d genes pass FDR < %.3g and |LFC| > %.2f (%s, %s)",
                    basename(path), nrow(df), n_all, opt$fdr, opt$lfc,
                    padj_col, lfc_col))
    df[[gene_col]]
}

message("Applying significance filter to each DE result:")
deseq2_genes <- read_sig(opt$deseq2_degs)
edger_genes  <- read_sig(opt$edger_degs)
limma_genes  <- read_sig(opt$limma_degs)

# Triple intersection: genes called significant by all three methods.
core <- Reduce(intersect, list(deseq2_genes, edger_genes, limma_genes))
write_tsv(tibble(gene_id=core), file.path(opt$outdir, "core_degs_3way.tsv"))

message(sprintf("DEGs: DESeq2=%d  edgeR=%d  limma-voom=%d  core=%d",
                length(deseq2_genes), length(edger_genes),
                length(limma_genes), length(core)))

# Venn
venn.plot <- venn.diagram(
    x = list(DESeq2=deseq2_genes, edgeR=edger_genes, `limma-voom`=limma_genes),
    filename=NULL, output=FALSE,
    fill=c("#2e75b6","#e6a23c","#67c23a"), alpha=0.5, lwd=1,
    main="Three-way DEG concordance"
)
pdf(file.path(opt$outdir, "venn_DESeq2_edgeR_limma.pdf"), width=5, height=5)
grid::grid.draw(venn.plot)
dev.off()

# Concordance summary TSV
summary_df <- tibble(
    method  = c("DESeq2","edgeR","limma-voom","core_3way"),
    n_degs  = c(length(deseq2_genes), length(edger_genes), length(limma_genes), length(core)),
    pct_in_core = c(
        length(intersect(deseq2_genes, core)) / max(length(deseq2_genes), 1),
        length(intersect(edger_genes,  core)) / max(length(edger_genes),  1),
        length(intersect(limma_genes,  core)) / max(length(limma_genes),  1),
        1.0
    )
)
write_tsv(summary_df, file.path(opt$outdir, "method_summary.tsv"))

# ---- Provenance ------------------------------------------------------------
write_json(list(
    script      = "concordance_analysis.R",
    inputs      = list(star_counts   = opt$star_counts,
                       hisat2_counts = opt$hisat2_counts,
                       salmon_counts = opt$salmon_counts,
                       deseq2_degs   = opt$deseq2_degs,
                       edger_degs    = opt$edger_degs,
                       limma_degs    = opt$limma_degs),
    params      = list(fdr = opt$fdr, lfc = opt$lfc),
    n_common_genes   = length(common_genes),
    n_common_samples = length(common_samples),
    spearman    = list(star_hisat2   = rho_star_hisat2,
                       star_salmon   = rho_star_salmon,
                       hisat2_salmon = rho_hisat2_salmon),
    n_degs      = list(deseq2 = length(deseq2_genes),
                       edger  = length(edger_genes),
                       limma  = length(limma_genes),
                       core   = length(core)),
    timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version   = R.version.string
), file.path(opt$outdir, "concordance_provenance.json"),
   auto_unbox = TRUE, pretty = TRUE)

message("Concordance analysis complete. Output: ", opt$outdir)