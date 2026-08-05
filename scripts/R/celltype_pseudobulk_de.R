#!/usr/bin/env Rscript
# scripts/R/celltype_pseudobulk_de.R
#

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(DESeq2)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--sc_object", type = "character"),
    make_option("--samples",   type = "character"),
    make_option("--condition", type = "character", default = "tumor,normal"),
    make_option("--fdr",       type = "double",    default = 0.05),
    make_option("--lfc",       type = "double",    default = 1.0),
    make_option("--threads",   type = "integer",   default = 4),
    make_option("--out_dir",   type = "character", default = ".")
)))

dir.create(file.path(opt$out_dir, "celltype_de_results"),
           showWarnings = FALSE, recursive = TRUE)

sc <- readRDS(opt$sc_object)
md <- sc@meta.data
required <- c("celltype", "sample_id", "condition")
if (!all(required %in% colnames(md))) {
    stop("sc object missing one of: ", paste(required, collapse = ", "))
}

contrast <- strsplit(opt$condition, ",")[[1]]
stopifnot(length(contrast) == 2)

n_de_per_celltype <- list()
all_volcano_data  <- list()

for (ct in unique(md$celltype)) {
    cells_ct <- rownames(md)[md$celltype == ct]
    if (length(cells_ct) < 50) next

    counts_ct <- GetAssayData(sc, assay = "RNA", layer = "counts")[, cells_ct]
    samples_ct <- md[cells_ct, "sample_id"]
    cond_ct    <- md[cells_ct, "condition"]

    pb <- sapply(split(seq_along(samples_ct), samples_ct),
                 function(idx) rowSums(counts_ct[, idx, drop = FALSE]))
    pb <- as.matrix(pb)

    samp_cond <- vapply(colnames(pb),
                        function(s) unique(cond_ct[samples_ct == s])[1],
                        character(1))
    if (length(unique(samp_cond)) < 2) next
    if (sum(samp_cond == contrast[1]) < 3 ||
        sum(samp_cond == contrast[2]) < 3) next

    coldata <- data.frame(sample = colnames(pb),
                          condition = factor(samp_cond, levels = rev(contrast)))
    pb <- pb[rowSums(pb >= 5) >= 3, , drop = FALSE]
    if (nrow(pb) < 500) next

    dds <- DESeqDataSetFromMatrix(countData = pb,
                                  colData   = coldata,
                                  design    = ~ condition)
    dds <- DESeq(dds, parallel = opt$threads > 1)
    res <- as.data.frame(results(dds, contrast = c("condition", contrast[1], contrast[2])))
    res$gene_id  <- rownames(res)
    res$celltype <- ct
    res <- res[!is.na(res$padj), ]

    sig <- subset(res, padj < opt$fdr & abs(log2FoldChange) >= opt$lfc)
    n_de_per_celltype[[ct]] <- nrow(sig)

    write.table(res,
                file.path(opt$out_dir, "celltype_de_results",
                          paste0(gsub("[^A-Za-z0-9_]", "_", ct), ".tsv")),
                quote = FALSE, sep = "\t", row.names = FALSE)

    if (ct %in% c("epithelial", "Epithelial", "basal_epithelial")) {
        all_volcano_data[[ct]] <- res
    }
}

# ---- Fig 5D: epithelial volcano ---------------------------------------------
ep <- bind_rows(all_volcano_data)
if (nrow(ep) > 0) {
    ep$sig <- with(ep, padj < opt$fdr & abs(log2FoldChange) >= opt$lfc)
    top <- ep |>
        filter(sig) |>
        arrange(padj) |>
        head(30)

    p <- ggplot(ep, aes(log2FoldChange, -log10(padj))) +
        geom_point(aes(colour = sig), size = 0.6, alpha = 0.6) +
        scale_colour_manual(values = c("TRUE" = "#c0392b", "FALSE" = "grey70")) +
        geom_text_repel(data = top, aes(label = gene_id),
                        size = 2.5, max.overlaps = 30) +
        geom_hline(yintercept = -log10(opt$fdr), linetype = 2) +
        geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = 2) +
        labs(x = paste0("log2 FC (", contrast[1], " vs ", contrast[2], ")"),
             y = "-log10 adj. p",
             title = "Epithelial-cell DE (pseudo-bulk)") +
        theme_bw(base_size = 9) +
        theme(legend.position = "none")
    ggsave(file.path(opt$out_dir, "fig5d_epithelial_volcano.pdf"),
           p, width = 5, height = 5)
}

prov <- list(
    tool        = "DESeq2_pseudobulk",
    version     = as.character(packageVersion("DESeq2")),
    seurat      = as.character(packageVersion("Seurat")),
    contrast    = contrast,
    fdr         = opt$fdr,
    lfc         = opt$lfc,
    n_celltypes_tested = length(n_de_per_celltype),
    n_de_per_celltype  = n_de_per_celltype,
    timestamp   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
write_json(prov, file.path(opt$out_dir, "celltype_de_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)
