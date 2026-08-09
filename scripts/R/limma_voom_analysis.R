#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# limma-voom differential expression analysis for KazRNA-Pipe.
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
    library(jsonlite)
})

opt_list <- list(
    make_option("--counts",  type="character"),
    make_option("--meta",    type="character"),
    make_option("--outdir",  type="character", default="results/de/limma_voom"),
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

observed <- levels(droplevels(meta$condition))
if (length(observed) < 2)
    stop("Contrast variable 'condition' has a single observed level ('",
         paste(observed, collapse=", "), "') across ", nrow(meta),
         " samples. Differential expression requires at least two levels. ",
         "Refusing to fit a rank-deficient model.")


dge <- DGEList(counts=counts)
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

# ---- Machine-readable provenance -------------------------------------------
write_json(list(
    script     = "limma_voom_analysis.R",
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
    packages   = list(limma = as.character(packageVersion("limma")), edgeR = as.character(packageVersion("edgeR")))
), file.path(opt$outdir, "limma_voom_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
