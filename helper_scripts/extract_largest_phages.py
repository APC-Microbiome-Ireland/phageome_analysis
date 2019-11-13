from Bio import SeqIO
file = "../data/megaphage_contigs.fasta"
loci = ["SRS016200_NODE_2_length_551943_cov_20.3899", "SRS047100_NODE_1_length_449263_cov_16.0097", "SRS017007_NODE_1_length_547686_cov_137.51"]

for record in SeqIO.parse(file, "fasta"):
    for locus in loci:
        if record.id == locus:
            SeqIO.write(record, "../data/" + locus + ".fasta", "fasta")
