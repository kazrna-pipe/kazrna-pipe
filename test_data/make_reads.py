#!/usr/bin/env python3
"""
Generate matched paired-end reads sampled from the annotated gene intervals.

Two properties matter for the downstream analysis:

1. Depth is drawn per gene per sample from a negative binomial, so gene-wise
   dispersion varies across the mean range. DESeq2's estimateDispersionsFit()
   needs that spread; with uniform depth it aborts.

2. A defined subset of genes is differentially expressed between the two
   groups, so the differential-expression steps produce non-empty results and
   the CI assertions test something. Without signal, a green run only proves
   the code does not crash.

Reads are generated on the forward strand of each gene with R2 as the reverse
complement, i.e. an UNSTRANDED library. conf/test.config must therefore set
strandedness = 'unstranded'; featureCounts -s 2 assigns nothing to any feature.

Deterministic: fixed seed, so regeneration reproduces the files byte for byte.
"""
import sys
import random
import gzip

random.seed(7)

fa, outdir = sys.argv[1], sys.argv[2]
refs_dir = fa.rsplit("/", 1)[0]

# ---- Reference ------------------------------------------------------------
seq = []
with open(fa) as f:
    for line in f:
        if not line.startswith(">"):
            seq.append(line.strip())
genome = "".join(seq)

genes = []
with open(f"{refs_dir}/gene_intervals.tsv") as f:
    next(f)
    for line in f:
        gid, start, end = line.rstrip("\n").split("\t")
        genes.append((gid, int(start), int(end)))

# ---- Design ---------------------------------------------------------------
SAMPLES = [("test_T1", "tumor"), ("test_T2", "tumor"),
           ("test_N1", "normal"), ("test_N2", "normal")]

N_DE = 40                 # genes with a real group difference
FOLD_CHANGE = 4.0         # up in tumour
BASE_MEAN = 60            # expected fragments per gene per sample
DISPERSION = 0.25         # negative-binomial dispersion
RL = 100                  # read length
FRAG = 300                # fragment length

de_genes = set(g[0] for g in genes[:N_DE])
comp = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(comp)[::-1]


def nb_sample(mean, disp):
    """Negative binomial via gamma-Poisson mixture, so dispersion varies."""
    shape = 1.0 / disp
    scale = mean * disp
    lam = random.gammavariate(shape, scale)
    # Knuth's Poisson, adequate at these means
    l, k, p = pow(2.718281828459045, -lam), 0, 1.0
    while True:
        k += 1
        p *= random.random()
        if p <= l:
            return k - 1


total_reads = {}
for sample, condition in SAMPLES:
    n_written = 0
    with gzip.open(f"{outdir}/{sample}_R1.fastq.gz", "wt") as r1, \
         gzip.open(f"{outdir}/{sample}_R2.fastq.gz", "wt") as r2:
        for gid, gs, ge in genes:
            mean = BASE_MEAN
            if gid in de_genes and condition == "tumor":
                mean *= FOLD_CHANGE
            depth = nb_sample(mean, DISPERSION)
            for _ in range(depth):
                if ge - gs < FRAG:
                    continue
                pos = random.randint(gs, ge - FRAG)
                frag = genome[pos:pos + FRAG]
                q = "I" * RL
                r1.write(f"@{sample}_{n_written}/1\n{frag[:RL]}\n+\n{q}\n")
                r2.write(f"@{sample}_{n_written}/2\n{revcomp(frag[-RL:])}\n+\n{q}\n")
                n_written += 1
    total_reads[sample] = n_written

print("reads written:")
for s, n in total_reads.items():
    print(f"  {s}: {n} pairs")
print(f"{len(de_genes)} of {len(genes)} genes are differentially expressed "
      f"({FOLD_CHANGE}x up in tumour)")
