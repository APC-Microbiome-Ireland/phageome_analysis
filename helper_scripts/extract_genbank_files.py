from Bio import SeqIO
import os
file = "../data/circular_jumbo_phage_genomes_RAST.gbk"

loci = ["ERR589420_NODE_6_length_229992_cov_13.7426_182", "SRS019607_NODE_9_length_235695_cov_219.599_15",
"SRS044662_NODE_25_length_232905_cov_15.7753_251", "SRS078425_NODE_1_length_226348_cov_18.355_86"]

output_prefix = "../data/mcp/"
for record in SeqIO.parse(file, "genbank"):
    for locus in loci:
        prot = locus.split("_")[-1]
        locus = locus.split("_" + prot)[0]
        if record.id == locus:
            #record.features[int(prot)].type = "MCP"
            if not os.path.exists(output_prefix):
                os.makedirs(output_prefix)
            SeqIO.write(record, output_prefix + locus + ".gbk", "genbank")
