from Bio import SeqIO
import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna

folder = "../data/mcp/"

loci_rev = ["SRS078425_NODE_1_length_226348_cov_18.355", "SRS019607_NODE_9_length_235695_cov_219.599",
"SRS044662_NODE_25_length_232905_cov_15.7753"]

loci = ["SRS078425_NODE_1_length_226348_cov_18.355", "SRS019607_NODE_9_length_235695_cov_219.599",
"ERR589420_NODE_6_length_229992_cov_13.7426",
"SRS044662_NODE_25_length_232905_cov_15.7753"]


for locus in loci:
    file = folder + locus + ".gbk"
    for record in SeqIO.parse(file, "genbank"):
        if locus in loci_rev:
            order_genome = []
            order_list = []
            # Reorder genome
            genome = list(str(record.seq))
            order_genome.extend(range(len(genome)))
            order_genome.reverse()
            reordered_genome = [genome[n] for n in order_genome]
            for i,nucl in enumerate(reordered_genome):
                if nucl == "A":
                    reordered_genome[i] = "T"
                elif nucl == "C":
                    reordered_genome[i] = "G"
                elif nucl == "G":
                    reordered_genome[i] = "C"
                elif nucl == "T":
                    reordered_genome[i] = "A"

            reordered_genome = "".join(reordered_genome)
            record.seq = Seq(reordered_genome, ambiguous_dna)

            # Change position of genes
            last_end_pos = record.features[len(record.features)-1].location.end.position
            for i in range(1, len(record.features)):
                start_pos = last_end_pos - record.features[i].location.end.position
                end_pos = last_end_pos - record.features[i].location.start.position
                new_strand = record.features[i].location.strand * -1
                record.features[i].location = FeatureLocation(start_pos, end_pos, strand = new_strand)

            order_list.extend(range(1, len(record.features)))
            order_list.reverse()
            order_list.insert(0,0)
            record.features = [record.features[j] for j in order_list]

        order_list = [0]
        order_genome = []
        # Get positions of MCP and end
        for i in range(len(record.features)):
            if record.features[i].type == "MCP":
                order_from = i
                mcp_start_pos = record.features[order_from].location.start.position
                last_end_pos = record.features[0].location.end.position
                order_genome.extend(range(mcp_start_pos, last_end_pos))
                order_genome.extend(range(mcp_start_pos))

        order_list.extend(range(order_from, len(record.features)))
        order_list.extend(range(1, order_from))
        features = [record.features[j] for j in order_list]

        # Change locations
        for j in range(len(features)):
            if record.features[j].type != "source":
                start_pos = features[j].location.start.position
                end_pos = features[j].location.end.position
                if start_pos >= mcp_start_pos:
                    start_pos = start_pos - mcp_start_pos
                    end_pos = end_pos - mcp_start_pos
                else:
                    start_pos = start_pos + (last_end_pos - mcp_start_pos)
                    end_pos = end_pos + (last_end_pos - mcp_start_pos)
                features[j].location = FeatureLocation(start_pos, end_pos, strand = features[j].location.strand)
        record.features = features

        # Reorder genome
        genome = list(str(record.seq))
        reordered_genome = [genome[n] for n in order_genome]
        reordered_genome = "".join(reordered_genome)
        record.seq = Seq(reordered_genome, ambiguous_dna)

        SeqIO.write(record, folder + locus + "_reordered.gbk", "genbank")
