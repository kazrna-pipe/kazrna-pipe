#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# edgeR differential expression analysis for KazRNA-Pipe.
#
# Parallels deseq2_analysis.R for cross-method concordance comparisons
# reported in Figure 2B.
#
# Bioconductor edgeR 4.2.1 - Robinson MD, McCarthy DJ, Smyth GK. Bioinformatics 2010.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(edgeR)
    library(tibble)
    library(readr)
    library(dplyr)
    library(jsonlite)
})

opt_list <- list(
    make_option("--counts",  type="character"),
    make_option("--meta",    type="character"),
    make_option("--outdir",  type="character", default="results/de/edger"),
    make_option("--fdr",     type="double",    default=0.05),
    make_option("--lfc",     type="double",    default=1.0),
    make_option("--covariates", type="character", default="",
                help="Comma-separated covariate columns; silently ignored if absent from --meta"),
    make_option("--seed",    type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

counts <- as.matrix(read.table(opt$counts, header=TRUE, row.names=1,
                               sep="\t", check.names=FALSE))
meta <- read_csv(opt$meta, show_col_types=FALSE, comment="#") %>%
    as.data.frame() %>% { rownames(.) <- .$sample_id; . }

common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common, drop=FALSE]
meta   <- meta[common, , drop=FALSE]
meta$condition <- factor(meta$condition, levels=c("normal", "tumor"))

# A single observed level gives a rank-deficient design: NA coefficients, exit 0.
observed <- levels(droplevels(meta$condition))
if (length(observed) < 2)
    stop("Contrast variable 'condition' has a single observed level ('",
         paste(observed, collapse=", "), "') across ", nrow(meta),
         " samples. Differential expression requires at least two levels. ",
         "Refusing to fit a rank-deficient model.")


# Build DGEList and filter weakly expressed genes
dge <- DGEList(counts=counts, group=meta$condition)
keep <- filterByExpr(dge, group=meta$condition)
dge  <- dge[keep, , keep.lib.sizes=FALSE]
dge  <- calcNormFactors(dge, method="TMM")


# ---- Design ----------------------------------------------------------------
covars <- if (is.null(opt$covariates) || opt$covariates == "") {
    character(0)
} else {
    trimws(strsplit(opt$covariates, ",")[[1]])
}
dropped <- setdiff(covars, colnames(meta))
if (length(dropped))
    message("Covariate(s) not present in metadata, dropped: ", paste(dropped, collapse=", "))
covars <- intersect(covars, colnames(meta))
design_formula <- as.formula(paste("~", paste(c(covars, "condition"), collapse=" + ")))
message("Design: ", deparse(design_formula))
design <- model.matrix(design_formula, data=meta)
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmQLFit(dge, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef="conditiontumor")

tt <- topTags(qlf, n=Inf, sort.by="PValue")$table %>%
    as.data.frame() %>% rownames_to_column("gene_id")

write_tsv(tt, file.path(opt$outdir, "edger_results.tsv"))
sig <- tt %>% filter(!is.na(FDR), FDR < opt$fdr, abs(logFC) > opt$lfc)
write_tsv(sig, file.path(opt$outdir, "edger_significant.tsv"))
message(sprintf("edgeR DEGs at FDR < %.3g and |logFC| > %.2f: %d", opt$fdr, opt$lfc, nrow(sig)))

sink(file.path(opt$outdir, "edger_session.txt"))
cat("KazRNA-Pipe edgeR module\n========================\n\n")
cat("Run datetime: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep="")
print(sessionInfo()); sink()
message("edgeR module complete.")

# ---- Machine-readable provenance -------------------------------------------
write_json(list(
    script     = "edger_analysis.R",
    inputs     = list(counts = opt$counts, meta = opt$meta),
    input_md5  = list(counts = unname(tools::md5sum(opt$counts)),
                      meta   = unname(tools::md5sum(opt$meta))),
    design     = deparse(design_formula),
    params     = list(fdr = opt$fdr, lfc = opt$lfc,
                      covariates = opt$covariates, seed = opt$seed),
    n_samples  = ncol(counts),
    n_genes    = nrow(counts),
    timestamp  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version  = R.version.string,
    packages   = list(edgeR = as.character(packageVersion("edgeR")))
), file.path(opt$outdir, "edger_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
