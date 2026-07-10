#!/usr/bin/env python3
"""Generate a tiny synthetic reference: FASTA + GTF for 3 genes on one contig."""
import sys, random
random.seed(2025)
out = sys.argv[1]
GENES = [("GENE1", 1000, 2000), ("GENE2", 4000, 5500), ("GENE3", 8000, 9200)]
CONTIG_LEN = 10000

def randseq(n):
    return "".join(random.choice("ACGT") for _ in range(n))

# One contig with genic regions embedded
seq = list(randseq(CONTIG_LEN))
genome = "".join(seq)
with open(f"{out}/mini.fa", "w") as f:
    f.write(">chrTest\n")
    for i in range(0, len(genome), 60):
        f.write(genome[i:i+60] + "\n")

with open(f"{out}/gencode.v44.annotation.gtf", "w") as f:
    for gid, start, end in GENES:
        tid = gid + "T1"
        attr_g = f'gene_id "{gid}"; gene_name "{gid}";'
        attr_t = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gid}";'
        f.write(f"chrTest\ttest\tgene\t{start}\t{end}\t.\t+\t.\t{attr_g}\n")
        f.write(f"chrTest\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attr_t}\n")
        f.write(f"chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attr_t}\n")

with open(f"{out}/tx2gene.tsv", "w") as f:
    f.write("transcript_id\tgene_id\tgene_name\n")
    for gid, *_ in GENES:
        f.write(f"{gid}T1\t{gid}\t{gid}\n")
print("reference written:", out)
