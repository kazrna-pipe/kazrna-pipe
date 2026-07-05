#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Seurat v5 single-cell workflow for KazRNA-Pipe.
#
# CPU reference path. Mirrors the GPU rapids-singlecell workflow in
# scripts/python/rapids_singlecell_workflow.py at every stage so that
# clustering agreement (ARI, NMI, ASW reported in manuscript Fig 3C) is
# a like-for-like comparison.
#
# Stages: load -> QC -> SCTransform -> merge -> Harmony -> PCA -> neighbors ->
#         Leiden (multi-resolution) -> UMAP -> markers -> manual annotation.
#
# Per-stage wall-clock is recorded for manuscript Figure 3A.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)            # 5.1.0
    library(harmony)           # 1.2.0
    library(leidenAlg)         # 0.10.2 (leidenAlg backend) - or igraph for fallback
    library(tibble); library(readr); library(dplyr); library(jsonlite)
})

opt_list <- list(
    make_option("--input_dir",   type="character",
                help="Directory containing one subdirectory per sample, each with 10x filtered matrices"),
    make_option("--outdir",      type="character", default="results/sc/cpu"),
    make_option("--mt_threshold", type="double",   default=20),
    make_option("--min_features", type="integer",  default=200),
    make_option("--hvg_n",       type="integer",   default=3000),
    make_option("--n_pcs",       type="integer",   default=50),
    make_option("--resolutions", type="character", default="0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0"),
    make_option("--seed",        type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)
resolutions <- as.numeric(strsplit(opt$resolutions, ",")[[1]])

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
sample_dirs <- list.dirs(opt$input_dir, recursive=FALSE)
if (length(sample_dirs) == 0L) stop("No sample subdirectories under ", opt$input_dir)

seurat_list <- .time("load", {
    lapply(sample_dirs, function(d) {
        sid <- basename(d)
        mat <- Read10X(data.dir=d)
        CreateSeuratObject(counts=mat, project=sid, min.cells=3, min.features=opt$min_features)
    })
})
names(seurat_list) <- basename(sample_dirs)

# ---- 2. Quality control ---------------------------------------------------
seurat_list <- .time("qc", lapply(seurat_list, function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern="^MT-")
    subset(x, subset = nFeature_RNA > opt$min_features & percent.mt < opt$mt_threshold)
}))

# ---- 3. SCTransform per sample --------------------------------------------
seurat_list <- .time("sctransform", lapply(seurat_list, function(x) {
    SCTransform(x, vst.flavor="v2", verbose=FALSE, return.only.var.genes=FALSE)
}))

# ---- 4. Merge + PCA -------------------------------------------------------
seurat <- .time("merge", merge(seurat_list[[1]], y=seurat_list[-1], merge.data=TRUE))
DefaultAssay(seurat) <- "SCT"
VariableFeatures(seurat) <- SelectIntegrationFeatures(seurat_list, nfeatures=opt$hvg_n)

seurat <- .time("pca", RunPCA(seurat, npcs=opt$n_pcs, verbose=FALSE))

# ---- 5. Harmony integration ------------------------------------------------
seurat <- .time("harmony", RunHarmony(
    seurat, group.by.vars="orig.ident", reduction.use="pca",
    dims.use=1:opt$n_pcs, plot_convergence=FALSE
))

# ---- 6. Neighbors graph + Leiden at multiple resolutions ------------------
seurat <- .time("neighbors", FindNeighbors(seurat, reduction="harmony", dims=1:opt$n_pcs))

seurat <- .time("leiden", {
    for (r in resolutions) {
        seurat <- FindClusters(seurat, resolution=r, algorithm=4)  # 4 = Leiden
        seurat[[sprintf("leiden_res_%.2f", r)]] <- seurat$seurat_clusters
    }
    seurat
})

# ---- 7. UMAP --------------------------------------------------------------
seurat <- .time("umap", RunUMAP(seurat, reduction="harmony", dims=1:opt$n_pcs, seed.use=opt$seed))

# ---- 8. Marker genes (Wilcoxon, default resolution = 1.0) -----------------
Idents(seurat) <- "leiden_res_1.00"
markers <- .time("markers", FindAllMarkers(
    seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE
))

# ---- 9. Save --------------------------------------------------------------
saveRDS(seurat,  file.path(opt$outdir, "seurat_object.rds"))
write_tsv(markers, file.path(opt$outdir, "markers.tsv"))

# Export labels at every resolution for downstream agreement analysis
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
    n_cells        = ncol(seurat),
    n_genes        = nrow(seurat),
    n_samples      = length(sample_dirs),
    sessionInfo    = capture.output(sessionInfo())
)
write_json(provenance, file.path(opt$outdir, "provenance.json"),
           pretty=TRUE, auto_unbox=TRUE)

message("Seurat CPU workflow complete. Output in: ", opt$outdir)
