#!/usr/bin/env python3
"""Generate tiny 10x-format filtered_feature_bc_matrix.h5 files + a celltype annotation TSV."""
import sys, numpy as np, scipy.sparse as sp, h5py
np.random.seed(2025)
outdir = sys.argv[1]
N_GENES, N_CELLS = 500, 300
GENE_IDS   = np.array([f"ENSG{i:08d}".encode() for i in range(N_GENES)])
GENE_NAMES = np.array([f"GENE{i}".encode() for i in range(N_GENES)])

def write_10x_h5(path, n_cells):
    X = sp.random(N_GENES, n_cells, density=0.15, format="csc", dtype=np.float32)
    X.data = np.ceil(X.data * 30).astype(np.int32)
    barcodes = np.array([f"CELL{i:05d}-1".encode() for i in range(n_cells)])
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data",    data=X.data.astype(np.int32))
        g.create_dataset("indices", data=X.indices.astype(np.int64))
        g.create_dataset("indptr",  data=X.indptr.astype(np.int64))
        g.create_dataset("shape",   data=np.array([N_GENES, n_cells], dtype=np.int32))
        g.create_dataset("barcodes", data=barcodes)
        fg = g.create_group("features")
        fg.create_dataset("id",           data=GENE_IDS)
        fg.create_dataset("name",         data=GENE_NAMES)
        fg.create_dataset("feature_type", data=np.array([b"Gene Expression"] * N_GENES))
        fg.create_dataset("genome",       data=np.array([b"chrTest"] * N_GENES))

for s, n in [("sc_sample_1", 300), ("sc_sample_2", 250)]:
    write_10x_h5(f"{outdir}/{s}_filtered_feature_bc_matrix.h5", n)

# shared cell-type annotation table (barcode -> celltype)
cts = ["Epithelial", "Fibroblast", "Tcell", "Myeloid", "Bcell", "Endothelial"]
with open(f"{outdir}/celltype_annotation_chr22.tsv", "w") as f:
    f.write("barcode\tcelltype\n")
    for i in range(300):
        f.write(f"CELL{i:05d}-1\t{cts[i % len(cts)]}\n")
print("sc .h5 + annotation written")
