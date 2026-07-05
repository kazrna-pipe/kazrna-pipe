#!/usr/bin/env Rscript
# scripts/R/deconv_agreement.R
#
# Compares cell-type proportions from BayesPrism, MuSiC, and (optionally)
# CIBERSORTx. Produces Figure 5B (correlation heatmap + scatter panel).

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ComplexHeatmap)
    library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--bayesprism", type = "character"),
    make_option("--music",      type = "character"),
    make_option("--cibersortx", type = "character", default = NULL),
    make_option("--out_dir",    type = "character", default = ".")
)))

read_prop <- function(path) {
    m <- as.matrix(read.table(path, header = TRUE, sep = "\t",
                              row.names = 1, check.names = FALSE))
    # Rows = samples, columns = cell types. Some tools write the transpose.
    if (sum(m, na.rm = TRUE) / nrow(m) < 0.5) m <- t(m)
    m
}

bp <- read_prop(opt$bayesprism)
mu <- read_prop(opt$music)

common_samples   <- Reduce(intersect, list(rownames(bp), rownames(mu)))
common_celltypes <- Reduce(intersect, list(colnames(bp), colnames(mu)))

if (!is.null(opt$cibersortx)) {
    cb <- read_prop(opt$cibersortx)
    common_samples   <- intersect(common_samples,   rownames(cb))
    common_celltypes <- intersect(common_celltypes, colnames(cb))
}

bp <- bp[common_samples, common_celltypes, drop = FALSE]
mu <- mu[common_samples, common_celltypes, drop = FALSE]

corr_df <- data.frame(
    celltype = common_celltypes,
    bayesprism_vs_music = sapply(common_celltypes, function(ct)
        cor(bp[, ct], mu[, ct], method = "pearson"))
)

if (!is.null(opt$cibersortx)) {
    cb <- cb[common_samples, common_celltypes, drop = FALSE]
    corr_df$bayesprism_vs_cibersortx <- sapply(common_celltypes, function(ct)
        cor(bp[, ct], cb[, ct], method = "pearson"))
    corr_df$music_vs_cibersortx      <- sapply(common_celltypes, function(ct)
        cor(mu[, ct], cb[, ct], method = "pearson"))
}

write.table(corr_df, file.path(opt$out_dir, "deconv_correlation.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Figure 5B
long <- as.data.frame(bp) |> tibble::rownames_to_column("sample") |>
    pivot_longer(-sample, names_to = "celltype", values_to = "bp") |>
    inner_join(
        as.data.frame(mu) |> tibble::rownames_to_column("sample") |>
            pivot_longer(-sample, names_to = "celltype", values_to = "mu"),
        by = c("sample", "celltype")
    )

p <- ggplot(long, aes(bp, mu)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "BayesPrism", y = "MuSiC",
         title = "Cross-method deconvolution agreement") +
    theme_bw(base_size = 9)
ggsave(file.path(opt$out_dir, "fig5b_deconv_corr.pdf"),
       p, width = 8, height = 6)

prov <- list(
    inputs = list(bayesprism = opt$bayesprism,
                  music      = opt$music,
                  cibersortx = opt$cibersortx),
    n_samples   = length(common_samples),
    n_celltypes = length(common_celltypes),
    timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(prov,
           file.path(opt$out_dir, "deconv_agreement_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)
