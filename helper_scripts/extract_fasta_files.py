from Bio import SeqIO
import os
file = "../data/megaphage_contigs.fasta"
#loci = ["SRS015941_NODE_2_length_315556_cov_47.6479", "SRS014470_NODE_1_length_389119_cov_8.88358"]
loci = ["SRS078431_NODE_11_length_258293_cov_26.0675", "SRS017304_NODE_5_length_258288_cov_116.449"]
cluster = "VC_1419_6"

for record in SeqIO.parse(file, "fasta"):
    for locus in loci:
        if record.id == locus:
            SeqIO.write(record, "../data/" + cluster + "/" + locus + ".fasta", "fasta")
