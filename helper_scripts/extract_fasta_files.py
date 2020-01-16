from Bio import SeqIO
import os
file = "../data/megaphage_contigs.fasta"
# loci = ["SRS016200_NODE_2_length_551943_cov_20.3899", "SRS047100_NODE_1_length_449263_cov_16.0097", "SRS017007_NODE_1_length_547686_cov_137.51"]
loci = ["SRS013252_NODE_7_length_269911_cov_143.484", "SRS078677_NODE_7_length_263796_cov_15.3008", "SRS142483_NODE_8_length_267299_cov_53.3807"]
cluster = "VC_3069_2"
output_prefix = "../data/" + cluster + "/fasta/"

for record in SeqIO.parse(file, "fasta"):
    for locus in loci:
        if record.id == locus:
            if not os.path.exists(output_prefix):
                os.makedirs(output_prefix)
            SeqIO.write(record, output_prefix + locus + ".fasta", "fasta")
