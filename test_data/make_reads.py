#!/usr/bin/env python3
"""Generate matched paired-end reads sampled from transcript regions."""
import sys, random, gzip
random.seed(7)
fa, outdir = sys.argv[1], sys.argv[2]
seq = []
with open(fa) as f:
    for line in f:
        if not line.startswith(">"):
            seq.append(line.strip())
genome = "".join(seq)
GENES = [(1000, 2000), (4000, 5500), (8000, 9200)]
RL, N = 100, 2000
comp = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(comp)[::-1]
for sample in ["test_T1", "test_N1", "test_T2", "test_N2"]:
    with gzip.open(f"{outdir}/{sample}_R1.fastq.gz", "wt") as r1, \
         gzip.open(f"{outdir}/{sample}_R2.fastq.gz", "wt") as r2:
        for i in range(N):
            gs, ge = random.choice(GENES)
            pos = random.randint(gs, ge - 300)
            frag = genome[pos:pos+300]
            read1 = frag[:RL]
            read2 = revcomp(frag[-RL:])
            q = "I" * RL
            r1.write(f"@{sample}_{i}/1\n{read1}\n+\n{q}\n")
            r2.write(f"@{sample}_{i}/2\n{read2}\n+\n{q}\n")
print("reads written for 4 samples")
