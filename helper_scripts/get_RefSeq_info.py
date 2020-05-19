from Bio import SeqIO
from os import listdir

refseq_file = "../data/refseq_r99/refseq_r99_jumbophage_genbank.gb"

with open("../data/refseq_r99/refseq_r99_jumbophage_info.txt", "w") as f:
    f.write("name" + "\t" + "length" + "\t" + "taxonomy" + "\t" + "host" + "\t" + "circular" + "\t" + "reference" + "\t" + "environment" + "\n")
    for record in SeqIO.parse(refseq_file, "genbank"):
        f.write(record.id + "\t" + str(len(record)) + "\t" + "_".join(record.annotations.get("taxonomy")) + \
        "\t" + record.annotations.get("source") + "\t" + record.annotations.get("topology") + "\t" + record.annotations.get("references")[0].title + "\t")
        if record.features[0].qualifiers.get("isolation_source") is None or record.features[0].qualifiers.get("country") is None:
            f.write("NA" + "\n")
        else:
            f.write(record.features[0].qualifiers.get("isolation_source")[0] + ", " + record.features[0].qualifiers.get("country")[0] + "\n")
