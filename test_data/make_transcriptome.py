#!/usr/bin/env python3
"""Extract transcript sequences from mini.fa + GTF for Salmon indexing."""
import sys
refs = sys.argv[1]
# parse contig
seq = []
with open(f"{refs}/mini.fa") as f:
    for line in f:
        if not line.startswith(">"):
            seq.append(line.strip())
genome = "".join(seq)
# parse transcripts (single-exon)
txs = []
with open(f"{refs}/gencode.v44.annotation.gtf") as f:
    for line in f:
        p = line.split("\t")
        if len(p) > 2 and p[2] == "exon":
            start, end = int(p[3]), int(p[4])
            tid = [a for a in p[8].split(";") if "transcript_id" in a][0].split('"')[1]
            txs.append((tid, genome[start-1:end]))
with open(f"{refs}/mini.transcripts.fa", "w") as f:
    for tid, s in txs:
        f.write(f">{tid}\n")
        for i in range(0, len(s), 60):
            f.write(s[i:i+60] + "\n")
print("transcriptome written:", len(txs), "transcripts")
