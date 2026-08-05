#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Merge per-sample featureCounts output into one gene x sample count matrix.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(readr); library(dplyr); library(purrr); library(tibble)
})

opt_list <- list(
    make_option("--counts", type = "character",
                help = "Comma-separated list of per-sample featureCounts files"),
    make_option("--sample_ids", type = "character", default = NULL,
                help = "Comma-separated sample IDs, same order as --counts. If omitted, IDs are derived from the filenames."),
    make_option("--out", type = "character", default = "merged_counts.tsv",
                help = "Output TSV path [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

files <- trimws(strsplit(opt$counts, ",")[[1]])
files <- files[nzchar(files)]
if (!length(files)) stop("No input files given to --counts")

missing <- files[!file.exists(files)]
if (length(missing)) stop("File(s) not found: ", paste(missing, collapse = ", "))

if (!is.null(opt$sample_ids) && nzchar(opt$sample_ids)) {
    ids <- trimws(strsplit(opt$sample_ids, ",")[[1]])
    if (length(ids) != length(files))
        stop(sprintf("--counts has %d entries but --sample_ids has %d",
                     length(files), length(ids)))
} else {
    # Strip .counts.tsv / .featurecounts.txt and any remaining extension.
    ids <- sub("\\.(counts\\.tsv|featurecounts\\.txt)$", "", basename(files))
    ids <- sub("\\.[^.]+$", "", ids)
}

read_one <- function(path, sid) {
    df <- read_tsv(path, comment = "#", show_col_types = FALSE,
                   progress = FALSE)
    if (!"Geneid" %in% names(df))
        stop("No Geneid column in ", path, " - is this really featureCounts output?")
    # The count column is the last one; the five before it are annotation.
    count_col <- names(df)[ncol(df)]
    df %>%
        select(gene_id = Geneid, !!sid := all_of(count_col)) %>%
        mutate(!!sid := as.integer(.data[[sid]]))
}

message("Merging ", length(files), " featureCounts files")
mats <- map2(files, ids, read_one)

merged <- reduce(mats, full_join, by = "gene_id")
merged[is.na(merged)] <- 0L
merged <- merged %>% arrange(gene_id)

write_tsv(merged, opt$out)

message(sprintf("Wrote %s: %d genes x %d samples (%s)",
                opt$out, nrow(merged), length(ids), paste(ids, collapse = ", ")))
