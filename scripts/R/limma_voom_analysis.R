#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# limma-voom differential expression analysis for KazRNA-Pipe.
#
# Third leg of the cross-method concordance comparison (manuscript Fig 2B).
#
# Bioconductor limma 3.60.6 - Law CW, Chen Y, Shi W, Smyth GK. Genome Biol 2014.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(limma)
    library(edgeR)
    library(tibble)
    library(readr)
    library(dplyr)
})

opt_list <- list(
    make_option("--counts",  type="character"),
    make_option("--meta",    type="character"),
    make_option("--outdir",  type="character", default="results/de/limma_voom"),
    make_option("--fdr",     type="double",    default=0.05),
    make_option("--lfc",     type="double",    default=1.0),
    make_option("--seed",    type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

counts <- as.matrix(read.table(opt$counts, header=TRUE, row.names=1,
                               sep="\t", check.names=FALSE))
meta <- read_csv(opt$meta, show_col_types=FALSE) %>%
    as.data.frame() %>% { rownames(.) <- .$sample_id; . }

common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common, drop=FALSE]
meta   <- meta[common, , drop=FALSE]
meta$condition <- factor(meta$condition, levels=c("normal", "tumor"))

dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, group=meta$condition)
dge  <- dge[keep, , keep.lib.sizes=FALSE]
dge  <- calcNormFactors(dge, method="TMM")

design <- model.matrix(~ age + sex + condition, data=meta)
v <- voom(dge, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust=TRUE)

tt <- topTable(fit, coef="conditiontumor", number=Inf, sort.by="P") %>%
    as.data.frame() %>% rownames_to_column("gene_id")

write_tsv(tt, file.path(opt$outdir, "limma_voom_results.tsv"))
sig <- tt %>% filter(!is.na(adj.P.Val), adj.P.Val < opt$fdr, abs(logFC) > opt$lfc)
write_tsv(sig, file.path(opt$outdir, "limma_voom_significant.tsv"))
message(sprintf("limma-voom DEGs at FDR < %.3g and |logFC| > %.2f: %d", opt$fdr, opt$lfc, nrow(sig)))

sink(file.path(opt$outdir, "limma_voom_session.txt"))
cat("KazRNA-Pipe limma-voom module\n=============================\n\n")
cat("Run datetime: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep="")
print(sessionInfo()); sink()
message("limma-voom module complete.")
