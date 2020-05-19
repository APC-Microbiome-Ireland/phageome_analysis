from Bio import SeqIO
import os
input_fasta = "../data/virsorter_positive.fasta"
output_fasta = "../data/jumbophage_contigs.fasta"
input = "../data/jumbophage_contigs_for_annot.txt"

headers = []
with open(input, "r") as f:
    for line in f:
        headers.append(line.split("\n")[0])

records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    if record.id in headers:
        records.append(record)

SeqIO.write(records, output_fasta, "fasta")
