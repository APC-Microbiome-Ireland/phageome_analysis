from Bio import SeqIO
loci = ["SRS016200_NODE_2_length_551943_cov_20.3899", "SRS047100_NODE_1_length_449263_cov_16.0097", "SRS017007_NODE_1_length_547686_cov_137.51"]
fasta_file = '../data/megaphage_contigs_translations2.fna'

for locus in loci:
    jumbophages = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id.rsplit('_',1)[0] == locus:
            jumbophages.append(record)
            SeqIO.write(record, "../data/proteins/" + record.id + ".faa", "fasta")
