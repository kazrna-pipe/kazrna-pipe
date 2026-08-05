#!/usr/bin/env python3
"""
Generate a synthetic reference: FASTA + GTF + tx2gene for N_GENES genes on one
contig.

200 genes is the smallest size at which DESeq2 behaves as it does on real
data: estimateDispersionsFit() needs enough points to fit a mean-dispersion
trend, and aborts below roughly ten. A fixture smaller than that exercises the
plumbing while skipping the statistics.

Deterministic: fixed seed, so regeneration reproduces the files byte for byte.
"""
import sys
import random

random.seed(2025)
out = sys.argv[1]

N_GENES = 200
GENE_LEN = 1200          # exonic length per gene
GAP = 300                # intergenic spacing
FLANK = 2000             # padding at each end of the contig

CONTIG_LEN = FLANK * 2 + N_GENES * (GENE_LEN + GAP)


def randseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


genome = randseq(CONTIG_LEN)

with open(f"{out}/mini.fa", "w") as f:
    f.write(">chrTest\n")
    for i in range(0, len(genome), 60):
        f.write(genome[i:i + 60] + "\n")

genes = []
pos = FLANK
for i in range(1, N_GENES + 1):
    gid = f"GENE{i:03d}"
    start = pos
    end = pos + GENE_LEN - 1
    genes.append((gid, start, end))
    pos = end + GAP + 1

with open(f"{out}/gencode.v44.annotation.gtf", "w") as f:
    for gid, start, end in genes:
        tid = gid + "T1"
        attr_g = f'gene_id "{gid}"; gene_name "{gid}";'
        attr_t = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gid}";'
        f.write(f"chrTest\ttest\tgene\t{start}\t{end}\t.\t+\t.\t{attr_g}\n")
        f.write(f"chrTest\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attr_t}\n")
        f.write(f"chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attr_t}\n")

with open(f"{out}/tx2gene.tsv", "w") as f:
    f.write("transcript_id\tgene_id\tgene_name\n")
    for gid, *_ in genes:
        f.write(f"{gid}T1\t{gid}\t{gid}\n")

# Written for make_reads.py so the two stay consistent.
with open(f"{out}/gene_intervals.tsv", "w") as f:
    f.write("gene_id\tstart\tend\n")
    for gid, start, end in genes:
        f.write(f"{gid}\t{start}\t{end}\n")

print(f"reference written: {out} ({N_GENES} genes, contig {CONTIG_LEN} bp)")
