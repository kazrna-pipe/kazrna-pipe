#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# fetch_tcga.R
#
# Fetches TCGA-ESCA squamous-histology samples (manuscript Section 3.1).
# Uses TCGAbiolinks to query the GDC and STAR - Counts workflow outputs
# for direct comparability with the in-house pipeline.
#
# The exact GDCquery call is preserved here verbatim because the GDC is
# updated periodically; re-running this in 2027 may return more samples.
# The manifest captured at the time of the manuscript is saved to
# data/snapshots/TCGA-ESCA_gdc_manifest_<date>.tsv for forensic reference.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(TCGAbiolinks)        # >= 2.32.0
    library(SummarizedExperiment)
    library(readr); library(dplyr); library(tibble); library(jsonlite)
})

opt_list <- list(
    make_option("--project",   type="character", default="TCGA-ESCA"),
    make_option("--histology", type="character",
                default="Squamous Cell Neoplasms",
                help="primary_diagnosis filter; restricts to SCC."),
    make_option("--workflow",  type="character", default="STAR - Counts"),
    make_option("--outdir",    type="character", default="data/raw/TCGA-ESCA")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "..", "..", "snapshots"),
           recursive=TRUE, showWarnings=FALSE)

message("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), "] Querying GDC for ",
        opt$project, " (", opt$histology, ")...")

query <- GDCquery(
    project          = opt$project,
    data.category    = "Transcriptome Profiling",
    data.type        = "Gene Expression Quantification",
    workflow.type    = opt$workflow,
    sample.type      = c("Primary Tumor", "Solid Tissue Normal"),
    experimental.strategy = "RNA-Seq"
)

# Filter on histology (squamous only) after the initial query because the
# field has to come from the clinical metadata, not the manifest.
manifest <- getResults(query) %>% as_tibble()
clinical <- GDCquery_clinic(opt$project, type="clinical") %>% as_tibble()

squamous_cases <- clinical %>%
    filter(grepl(opt$histology, primary_diagnosis, ignore.case=TRUE)) %>%
    pull(submitter_id)

manifest_filtered <- manifest %>% filter(cases.submitter_id %in% squamous_cases)
message(sprintf("Kept %d files from %d squamous cases (started with %d in the project).",
                nrow(manifest_filtered), length(unique(manifest_filtered$cases.submitter_id)),
                nrow(manifest)))

# ---- Snapshot the manifest for the audit trail ----------------------------
snapshot_path <- file.path(opt$outdir, "..", "..", "snapshots",
                           sprintf("%s_gdc_manifest_%s.tsv",
                                   opt$project, format(Sys.Date())))
write_tsv(manifest_filtered, snapshot_path)
message("Manifest snapshot: ", snapshot_path)

# Rebuild the query with only the squamous file IDs and download.
query2 <- query
query2$results[[1]] <- manifest_filtered

GDCdownload(query2, directory = opt$outdir, files.per.chunk = 25)
se <- GDCprepare(query2, directory = opt$outdir)

saveRDS(se, file.path(opt$outdir, "TCGA-ESCA_squamous_SE.rds"))

# Export count matrix and sample metadata as plain TSV for the pipeline.
counts <- assay(se, "unstranded")
write.table(counts, file.path(opt$outdir, "TCGA-ESCA_squamous_counts.tsv"),
            sep="\t", quote=FALSE)

meta <- colData(se) %>%
    as.data.frame() %>% rownames_to_column("sample_id") %>%
    select(sample_id, barcode, sample_type, tissue_type) %>%
    mutate(condition = ifelse(sample_type == "Primary Tumor", "tumor", "normal"))
write_csv(meta, file.path(opt$outdir, "TCGA-ESCA_squamous_metadata.csv"))

provenance <- list(
    script    = "fetch_tcga.R",
    datetime  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    arguments = opt,
    n_files   = nrow(manifest_filtered),
    n_cases   = length(unique(manifest_filtered$cases.submitter_id)),
    workflow  = opt$workflow,
    sessionInfo = capture.output(sessionInfo())
)
write_json(provenance, file.path(opt$outdir, "provenance.json"),
           pretty=TRUE, auto_unbox=TRUE)

message("TCGA-ESCA fetch complete.")
