from Bio import SeqIO
from os import listdir

mcp_identifiers = ["major capsid protein", "Major capsid protein"]

headers = []
translations = []

# Add Al-Shayeb jumbophage mcps
alshayeb_folder = "../data/alshayeb/major_capsid_protein_phages/"
alshayeb_files = [alshayeb_folder + f for f in listdir(alshayeb_folder)]

for file in alshayeb_files:
    count = 1
    for record in SeqIO.parse(file, "genbank"):
        if len(record) >= 200000:
            for i in range(len(record.features)):
                if record.features[i].type == "CDS":
                    if any(mcp in record.features[i].qualifiers.get("product")[0] for mcp in mcp_identifiers):
                        headers.append(file.split(alshayeb_folder)[1].split(".gbk")[0] + "_" + str(count))
                        translations.append(record.features[i].qualifiers.get("translation")[0])
                        count += 1

# Add RefSeq jumbophage mcps
refseq_file = "../data/refseq_r99/refseq_r99_jumbophage_genbank.gb"
refseq_jumbo_file = "../data/refseq_r99/refseq_r99_jumbophage_headers.txt"

refseq_jumbophage_headers = []
with open(refseq_jumbo_file, "r") as f:
    for line in f:
        refseq_jumbophage_headers.append(line.split("\n")[0])

for record in SeqIO.parse(refseq_file, "genbank"):
    count = 1;
    if record.id in refseq_jumbophage_headers:
        for i in range(len(record.features)):
            if record.features[i].type == "CDS":
                if any(mcp in record.features[i].qualifiers.get("product")[0] for mcp in mcp_identifiers):
                    headers.append(record.id + "_" + str(count))
                    translations.append(record.features[i].qualifiers.get("translation")[0])
                    count += 1

# Add IMG/VR jumbophage mcps
imgvr_file = "../data/imgvr_uvigs_translations.faa"
imgvr_mcp_file = "../data/imgvr_mcp_headers.txt"

imgvr_mcp_headers = []
with open(imgvr_mcp_file, "r") as f:
    for line in f:
        imgvr_mcp_headers.append(line.split("\n")[0])

for record in SeqIO.parse(imgvr_file, "fasta"):
    if record.id in imgvr_mcp_headers:
        headers.append(record.id)
        translations.append(record.seq.__str__())

# Add our jumbophage mcps
jumbophage_file = "../data/jumbophage_contigs_translations.faa"
jumbophage_mcp_file = "../data/jumbophage_mcp_headers.txt"

jumbo_mcp_headers = []
with open(jumbophage_mcp_file, "r") as f:
    for line in f:
        jumbo_mcp_headers.append(line.split("\n")[0])

for record in SeqIO.parse(jumbophage_file, "fasta"):
    if record.id in jumbo_mcp_headers:
        headers.append(record.id)
        translations.append(record.seq.__str__())

# Write final major capsid protein file
with open("../data/major_capsid_proteins.faa", "w") as f:
    for i in range(len(headers)):
        f.write(">" + headers[i] + "\n" + translations[i] + "\n")
