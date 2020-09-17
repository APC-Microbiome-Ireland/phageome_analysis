library(data.table)
library(reshape2)
library(dplyr)

########## Read in and process raw data ####
#Contigs
contig_ids = read.table("../data/virsorter_positive.ids")
contig_ids = unique(as.character(contig_ids$V1))

########## Read in and process contig annotations ####
contig_data = data.frame(
  row.names = contig_ids,
  name = contig_ids,
  size = as.numeric(sapply(strsplit(contig_ids, split = "_"), "[[", 5)),
  coverage = as.numeric(sapply(strsplit(contig_ids, split = "_"), "[[", 7)),
  circular = FALSE,
  gc_content = NA,
  hit_to_vir_refseq = FALSE,
  hit_to_crasslike = FALSE,
  demovir = NA,
  vir_refseq = NA,
  crasslike = NA,
  ribo_prot_count = 0,
  pvogs_count = 0,
  integrase = FALSE,
  sample = NA,
  country = NA
)

circular = read.table("../data/circular_contigs.ids", header=FALSE)
vir_refseq = read.table("../data/refseq_annot.txt", sep = "\t", header=FALSE)
crasslike = read.table("../data/crasslike_annot.txt", sep = "\t", header=FALSE)
demovir = read.table("../data/DemoVir_assignments.txt", sep = "\t", header=TRUE, row.names = 1)
gc_content = read.table("../data/virsorter_positive_gc_content.txt", header=FALSE, row.names = 1)
ribo_prot_count = read.table("../data/ribosomal_prot_counts.txt", header=FALSE)
rownames(ribo_prot_count) = as.character(ribo_prot_count$V2)
pvogs_count = read.table("../data/pvogs_counts_per_contig.txt", header=FALSE)
rownames(pvogs_count) = as.character(pvogs_count$V2)
vcontact = read.csv("../data/genome_by_genome_overview.csv")
integrase = read.table("../data/integrase_contigs.ids", header=FALSE)
contig_data[rownames(contig_data) %in% circular$V1,"circular"] = TRUE
contig_data[rownames(contig_data) %in% integrase$V1,"integrase"] = TRUE
contig_data$gc_content = as.numeric(as.character(gc_content[rownames(contig_data), "V2"]))
contig_data[rownames(contig_data) %in% vir_refseq$V1,"hit_to_vir_refseq"] = TRUE
contig_data[rownames(contig_data) %in% crasslike$V1,"hit_to_crasslike"] = TRUE
contig_data$demovir = demovir[rownames(contig_data),"Family"]
contig_data$demovir = as.character(contig_data$demovir)
contig_data[contig_data$demovir %in%
              c("Flaviviridae", "Iridoviridae", "Poxviridae", "Retroviridae", "Mimiviridae", "Pithoviridae",
                "Ascoviridae", "Nudiviridae", "Baculoviridae", "Marseilleviridae", "Picornaviridae",
                "Phycodnaviridae", "Herpesviridae", "Alloherpesviridae"), "demovir"] = "Other"
contig_data[contig_data$demovir %in% c(NA, "Unassigned"), "demovir"] = "Unassigned"
contig_data[contig_data$hit_to_crasslike, "demovir"] = "crAss-like"
contig_data$demovir = as.factor(contig_data$demovir)
contig_data$demovir = factor(contig_data$demovir, levels =
                               c("Bicaudaviridae", "Tectiviridae", "Inoviridae", "Microviridae",
                                 "Siphoviridae", "Myoviridae", "Podoviridae", "crAss-like",
                                 "Other", "Unassigned")) #Manually reordering levels. Check carefully!!!

for(i in levels(vir_refseq$V1))
{
  contig_data[i, "vir_refseq"] = as.character(vir_refseq[vir_refseq$V1 == i,][1,2]) #use first hit
}
contig_data$vir_refseq = as.factor(contig_data$vir_refseq)

for(i in levels(crasslike$V1))
{
  contig_data[i, "crasslike"] = as.character(crasslike[crasslike$V1 == i,][1,2]) #use first hit
}
contig_data$crasslike = as.factor(contig_data$crasslike)

contig_data$ribo_prot_count = ribo_prot_count[rownames(contig_data), "V1"]
contig_data$pvogs_count = pvogs_count[rownames(contig_data), "V1"]

contig_data <- left_join(contig_data, vcontact[,c("Genome", "VC.Subcluster")], by = c("name"="Genome")) %>%
  rename(vcontact_cluster = VC.Subcluster)
contig_data$sample = sapply(strsplit(as.character(contig_data$name), split = "_"), "[[", 1)

contig_data[is.na(contig_data$ribo_prot_count), "ribo_prot_count"] = 0
contig_data[is.na(contig_data$pvogs_count), "pvogs_count"] = 0


########## Host information from CRISPR and tRNA matches ####
#Read in and process CRISPR BLAST data
crispr_bitscore_cutoff = 45
trna_bitscore_cutoff = 65
crispr_evalue_cutoff = 1E-05
trna_evalue_cutoff = 1E-10

crispr_hits2 =
  read.table(file = '../data/hmp_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)

crispr_hits = crispr_hits2

colnames(crispr_hits) = c("qseqid", "sseqid", "bitscore", "pident", "length", "qlen", "evalue", "sstart", "send", "sname")
crispr_hits = crispr_hits[crispr_hits$evalue < crispr_evalue_cutoff,]
crispr_hits = crispr_hits[crispr_hits$bitscore >= crispr_bitscore_cutoff,]#Remove low score hits
crispr_hits$sname = gsub("UNVERIFIED: ", "", crispr_hits$sname)
crispr_hits$sname = gsub("TPA: ", "", crispr_hits$sname)
crispr_hits$sname = gsub("Candidatus ", "", crispr_hits$sname)
crispr_hits$sname = gsub("^complete chromosome ", "", crispr_hits$sname)
crispr_host_species <- gsub("\\]", "", gsub(".*\\[", "", crispr_hits$sname))
crispr_hits$crispr_host <- gsub(" .*", "", crispr_host_species)

crispr_summary <- crispr_hits %>% group_by(sseqid, crispr_host) %>%
  summarise(sum_bitscore = sum(bitscore)) %>%
  group_by(sseqid) %>%
  filter(sum_bitscore == max(sum_bitscore)) %>%
  select(-sum_bitscore)

contig_data <- left_join(contig_data, crispr_summary, by = c("name"="sseqid"))

#Dereplicate contigs according to the set used in read alignment step
contig_data <- contig_data[!duplicated(contig_data$name),]

#Cleanup
rm(integrase,circular,vir_refseq,crasslike,demovir,gc_content,ribo_prot_count,
   pvogs_count,vcontact,crispr_hits2)
contig_data <- data.table(contig_data)

saveRDS(contig_data,  file = "../data/contig_data.RDS")
