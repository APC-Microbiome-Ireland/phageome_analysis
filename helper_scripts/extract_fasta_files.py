from Bio import SeqIO
import os
file = "../data/megaphage_contigs.fasta"
loci = ["SRS015941_NODE_2_length_315556_cov_47.6479", "SRS014470_NODE_1_length_389119_cov_8.88358"]

for record in SeqIO.parse(file, "fasta"):
    for locus in loci:
        if record.id == locus:
            SeqIO.write(record, "../data/" + locus + ".fasta", "fasta")
