#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Seurat v5 single-cell workflow for KazRNA-Pipe.
#
# CPU reference path. Mirrors the GPU rapids-singlecell workflow in
# scripts/python/rapids_singlecell_workflow.py at every stage so that
# clustering agreement (ARI, NMI, ASW reported in Fig 3C) is
# a like-for-like comparison.
#
# Stages: load -> QC -> normalize -> merge -> PCA -> Harmony -> neighbors ->
#         Leiden (multi-resolution) -> UMAP -> markers -> annotation.
#
# IMPORTANT: --normalization controls the variance-stabilisation method.
#   lognorm     : LogNormalize + vst HVG selection. Identical to the GPU path's
#                 default, and the ONLY setting valid for the CPU-vs-GPU
#                 comparison in Figure 3.
#   sctransform : SCTransform v2 regularised negative-binomial residuals. A
#                 different model; use for biological analysis if preferred,
#                 but do not mix it into the CPU/GPU benchmark.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)            # 5.1.0
    library(harmony)           # 1.2.0
    library(hdf5r)             # required by Read10X_h5()
    library(yaml)
    library(tibble); library(readr); library(dplyr); library(jsonlite)
})

opt_list <- list(
    make_option("--input_h5",      type="character", default=NULL,
                help="Comma-separated list of 10x .h5 files"),
    make_option("--sample_ids",    type="character", default=NULL,
                help="Comma-separated sample IDs, same order as --input_h5"),
    make_option("--input_dir",     type="character", default=NULL,
                help="Directory of per-sample 10x subdirectories (alternative to --input_h5)"),
    make_option("--markers",       type="character", default=NULL,
                help="YAML of cell-type marker genes for cluster annotation"),
    make_option("--outdir",        type="character", default="results/sc/cpu"),
    make_option("--normalization", type="character", default="lognorm",
                help="lognorm (matches GPU path) or sctransform [default %default]"),
    make_option("--mt_threshold",  type="double",    default=20),
    make_option("--min_features",  type="integer",   default=200),
    make_option("--hvg_n",         type="integer",   default=3000),
    make_option("--n_pcs",         type="integer",   default=50),
    make_option("--n_neighbors",   type="integer",   default=15),
    make_option("--resolutions",   type="character", default="0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0"),
    make_option("--threads",       type="integer",   default=1),
    make_option("--seed",          type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$input_h5) && is.null(opt$input_dir))
    stop("One of --input_h5 or --input_dir is required")
if (!is.null(opt$input_h5) && is.null(opt$sample_ids))
    stop("--input_h5 requires --sample_ids")
if (!opt$normalization %in% c("lognorm", "sctransform"))
    stop("--normalization must be 'lognorm' or 'sctransform'")

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)
resolutions <- as.numeric(strsplit(opt$resolutions, ",")[[1]])
primary_res <- if (1.0 %in% resolutions) 1.0 else resolutions[ceiling(length(resolutions)/2)]
cluster_key <- sprintf("leiden_res_%.2f", primary_res)

timing <- list()
.time <- function(label, expr) {
    t0 <- Sys.time()
    val <- force(expr)
    dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    timing[[label]] <<- dt
    message(sprintf("[stage] %-15s %7.1f s", label, dt))
    invisible(val)
}

# ---- 1. Load all samples ---------------------------------------------------
use_h5 <- !is.null(opt$input_h5)

if (use_h5) {
    files <- trimws(strsplit(opt$input_h5,   ",")[[1]])
    ids   <- trimws(strsplit(opt$sample_ids, ",")[[1]])
    if (length(files) != length(ids))
        stop(sprintf("--input_h5 has %d entries but --sample_ids has %d",
                     length(files), length(ids)))
    missing <- files[!file.exists(files)]
    if (length(missing)) stop("Input file(s) not found: ", paste(missing, collapse=", "))
} else {
    files <- list.dirs(opt$input_dir, recursive=FALSE)
    if (length(files) == 0L) stop("No sample subdirectories under ", opt$input_dir)
    ids   <- basename(files)
}

seurat_list <- .time("load", Map(function(f, s) {
    mat <- if (use_h5) Read10X_h5(f) else Read10X(data.dir = f)
    CreateSeuratObject(counts = mat, project = s,
                       min.cells = 3, min.features = opt$min_features)
}, files, ids))
names(seurat_list) <- ids

# ---- 2. Quality control ---------------------------------------------------
seurat_list <- .time("qc", lapply(seurat_list, function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern="^MT-")
    subset(x, subset = nFeature_RNA > opt$min_features & percent.mt < opt$mt_threshold)
}))

seurat_list <- .time("qc", lapply(seurat_list, function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern="^MT-")
    y <- subset(x, subset = nFeature_RNA > opt$min_features & percent.mt < opt$mt_threshold)
    if (ncol(y) == 0L)
        stop("QC removed every cell from sample '", Project(x), "'. ",
             "Median genes/cell before filtering: ", median(x$nFeature_RNA),
             "; --min_features is ", opt$min_features, ". Lower the threshold.")
    y
}))

# ---- 3. Normalization ------------------------------------------------------
# lognorm mirrors the GPU path exactly; sctransform is the alternative model.
seurat_list <- .time("normalize", lapply(seurat_list, function(x) {
    if (opt$normalization == "sctransform") {
        SCTransform(x, vst.flavor="v2", verbose=FALSE, return.only.var.genes=FALSE)
    } else {
        x <- NormalizeData(x, normalization.method="LogNormalize",
                           scale.factor=1e4, verbose=FALSE)
        FindVariableFeatures(x, selection.method="vst",
                             nfeatures=opt$hvg_n, verbose=FALSE)
    }
}))

# ---- 4. Merge + PCA -------------------------------------------------------
seurat <- .time("merge", merge(seurat_list[[1]], y=seurat_list[-1], merge.data=TRUE))

if (opt$normalization == "sctransform") {
    DefaultAssay(seurat) <- "SCT"
    VariableFeatures(seurat) <- SelectIntegrationFeatures(seurat_list, nfeatures=opt$hvg_n)
} else {
    DefaultAssay(seurat) <- "RNA"
    seurat <- JoinLayers(seurat)
    seurat <- FindVariableFeatures(seurat, selection.method="vst",
                                   nfeatures=opt$hvg_n, verbose=FALSE)
    seurat <- ScaleData(seurat, features=VariableFeatures(seurat), verbose=FALSE)
}

seurat <- .time("pca", RunPCA(seurat, npcs=opt$n_pcs, verbose=FALSE))

# ---- 5. Harmony integration ------------------------------------------------
seurat <- .time("harmony", RunHarmony(
    seurat, group.by.vars="orig.ident", reduction.use="pca",
    dims.use=1:opt$n_pcs, plot_convergence=FALSE
))

# ---- 6. Neighbors graph + Leiden at multiple resolutions ------------------
seurat <- .time("neighbors", FindNeighbors(seurat, reduction="harmony",
                                           dims=1:opt$n_pcs, k.param=opt$n_neighbors))

seurat <- .time("leiden", {
    for (r in resolutions) {
        seurat <- FindClusters(seurat, resolution=r, algorithm=4)  # 4 = Leiden
        seurat[[sprintf("leiden_res_%.2f", r)]] <- seurat$seurat_clusters
    }
    seurat
})

# ---- 7. UMAP --------------------------------------------------------------
seurat <- .time("umap", RunUMAP(seurat, reduction="harmony", dims=1:opt$n_pcs,
                                seed.use=opt$seed))

# ---- 8. Marker genes (Wilcoxon, at the primary resolution) ----------------
Idents(seurat) <- cluster_key
markers <- .time("markers", FindAllMarkers(
    seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE
))

# ---- 9. Cluster annotation from the marker YAML ---------------------------
celltypes <- .time("annotate", {
    clusters <- levels(Idents(seurat))
    if (is.null(opt$markers) || !file.exists(opt$markers)) {
        warning("No marker file at ", opt$markers, "; clusters left unassigned")
        tibble(cluster = clusters, cell_type = "unassigned", score = NA_real_)
    } else {
        marker_sets <- yaml::read_yaml(opt$markers)
        present <- lapply(marker_sets, function(g) intersect(g, rownames(seurat)))
        present <- present[lengths(present) > 0]

        if (length(present) == 0L) {
            warning("No marker genes found in the data; clusters left unassigned")
            tibble(cluster = clusters, cell_type = "unassigned", score = NA_real_)
        } else {
            seurat <- AddModuleScore(seurat, features = present,
                                     name = "ctscore", seed = opt$seed)
            score_cols <- paste0("ctscore", seq_along(present))
            mean_scores <- sapply(score_cols, function(col)
                tapply(seurat[[col]][[1]], Idents(seurat), mean))
            mean_scores <- matrix(mean_scores, nrow = length(clusters),
                                  dimnames = list(clusters, names(present)))
            tibble(
                cluster   = rownames(mean_scores),
                cell_type = colnames(mean_scores)[max.col(mean_scores, ties.method = "first")],
                score     = apply(mean_scores, 1, max)
            )
        }
    }
})

# ---- 10. Save --------------------------------------------------------------
saveRDS(seurat,  file.path(opt$outdir, "seurat_object.rds"))
write_tsv(markers,   file.path(opt$outdir, "markers.tsv"))
write_tsv(celltypes, file.path(opt$outdir, "celltype_assignments.tsv"))

labels_long <- lapply(resolutions, function(r) {
    col <- sprintf("leiden_res_%.2f", r)
    data.frame(
        cell_id    = colnames(seurat),
        resolution = r,
        cluster    = as.character(seurat[[col]][[1]]),
        stringsAsFactors = FALSE
    )
}) %>% bind_rows()
write_tsv(labels_long, file.path(opt$outdir, "labels.tsv"))

# UMAP coords
umap_df <- as.data.frame(Embeddings(seurat, "umap")) %>% rownames_to_column("cell_id")
write_tsv(umap_df, file.path(opt$outdir, "umap_coords.tsv"))

# Per-stage timing for Figure 3A
timing_df <- tibble(stage = names(timing), wall_seconds = unlist(timing), path = "cpu")
write_tsv(timing_df, file.path(opt$outdir, "timing.tsv"))

# Provenance JSON
provenance <- list(
    script         = "seurat_workflow.R",
    git_describe   = system("git describe --always --dirty 2>/dev/null || echo unknown",
                            intern=TRUE),
    datetime       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    hostname       = unname(Sys.info()["nodename"]),
    arguments      = opt,
    normalization  = opt$normalization,
    primary_cluster_key = cluster_key,
    n_cells        = ncol(seurat),
    n_genes        = nrow(seurat),
    n_samples      = length(ids),
    sample_ids     = ids,
    sessionInfo    = capture.output(sessionInfo())
)
write_json(provenance, file.path(opt$outdir, "provenance.json"),
           pretty=TRUE, auto_unbox=TRUE)

message("Seurat CPU workflow complete. Output in: ", opt$outdir)
