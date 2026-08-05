#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# BayesPrism deconvolution for KazRNA-Pipe.
#
# Primary deconvolution method (selected on the basis of highest cross-method
# agreement). Produces cell-type proportions for
# Figure 4B and the cell-type-resolved DEGs used to identify the
# 842 epithelial-specific DEGs reported in Figure 4C.
#
# BayesPrism 2.2 - Chu T, Wang Z, Pe'er D, Danko CG. Nat Cancer 2022.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(BayesPrism)        # 2.2
    library(Seurat)            # 5.1.0
    library(Matrix)
    library(readr); library(dplyr); library(tibble); library(jsonlite)
})

opt_list <- list(
    make_option("--bulk_counts",     type="character",
                help="Gene x sample raw count matrix (TSV)"),
    make_option("--sc_reference",    type="character",
                help="Path to a Seurat .rds object containing the GSE160269 reference"),
    make_option("--celltype_column", type="character", default="cell_type",
                help="Metadata column in the Seurat object with cell-type labels"),
    make_option("--outdir",          type="character", default="results/deconv/bayesprism"),
    make_option("--n_iter",          type="integer",   default=5000L),
    make_option("--n_signatures",    type="integer",   default=50L),
    make_option("--n_cores",         type="integer",   default=8L),
    make_option("--seed",            type="integer",   default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# ---- 1. Load bulk counts --------------------------------------------------
bulk <- as.matrix(read.table(opt$bulk_counts, header=TRUE, row.names=1,
                             sep="\t", check.names=FALSE))
storage.mode(bulk) <- "integer"
message(sprintf("Loaded bulk matrix: %d genes x %d samples", nrow(bulk), ncol(bulk)))

# ---- 2. Load sc reference and build a per-cell-type expression matrix -----
ref <- readRDS(opt$sc_reference)
if (!opt$celltype_column %in% colnames(ref@meta.data)) {
    stop("Cell-type column '", opt$celltype_column, "' not present in reference metadata.")
}
cell_types <- as.character(ref@meta.data[[opt$celltype_column]])
sc_counts <- GetAssayData(ref, layer="counts")   # genes x cells, sparse

# BayesPrism wants cells x genes
sc_counts <- t(sc_counts)
message(sprintf("Reference matrix: %d cells x %d genes; %d cell types.",
                nrow(sc_counts), ncol(sc_counts), length(unique(cell_types))))

# ---- 3. Align gene names ---------------------------------------------------
shared_genes <- intersect(rownames(bulk), colnames(sc_counts))
if (length(shared_genes) < 5000L) {
    warning("Only ", length(shared_genes),
            " genes shared between bulk and sc reference; results may be unstable.")
}
bulk      <- bulk[shared_genes, , drop=FALSE]
sc_counts <- sc_counts[, shared_genes, drop=FALSE]

# ---- 4. Build the Prism object --------------------------------------------
# BayesPrism takes a malignant-vs-stromal split. For ESCC we treat
# 'Epithelial' as the malignant compartment if present.
malignant_label <- if ("Epithelial" %in% cell_types) "Epithelial" else NA

prism <- new.prism(
    reference   = sc_counts,
    mixture     = t(bulk),                   # samples x genes
    input.type  = "count.matrix",
    cell.type.labels = cell_types,
    cell.state.labels = cell_types,         # cell-type-level only; no sub-states
    key         = malignant_label,           # malignant compartment label
    outlier.cut = 0.01,
    outlier.fraction = 0.1
)

# ---- 5. Run inference -----------------------------------------------------
res <- run.prism(prism = prism, n.cores = opt$n_cores)

# ---- 6. Extract proportions -----------------------------------------------
theta <- get.fraction(bp = res, which.theta = "final", state.or.type = "type")
prop_df <- as.data.frame(theta) %>% rownames_to_column("sample_id")
write_tsv(prop_df, file.path(opt$outdir, "bayesprism_proportions.tsv"))

# Per-cell-type expression Z (Z.tg)
z_tg <- get.exp(bp = res, state.or.type = "type")
saveRDS(z_tg, file.path(opt$outdir, "bayesprism_cellexp_Z.rds"))

# ---- 7. Provenance --------------------------------------------------------
provenance <- list(
    script        = "bayesprism_deconv.R",
    datetime      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    hostname      = unname(Sys.info()["nodename"]),
    arguments     = opt,
    n_bulk        = ncol(bulk),
    n_genes_shared = length(shared_genes),
    n_cells_ref   = nrow(sc_counts),
    cell_types    = unique(cell_types),
    sessionInfo   = capture.output(sessionInfo())
)
write_json(provenance, file.path(opt$outdir, "provenance.json"),
           pretty=TRUE, auto_unbox=TRUE)

message("BayesPrism complete. Output: ", opt$outdir)
