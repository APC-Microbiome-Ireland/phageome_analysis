library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)

########## Read in metadata and define parameters ####
metadata = read.csv("../../data/metadata_v2.csv", stringsAsFactors = FALSE)
rownames(metadata) = metadata$ID
metadata$Health[metadata$Health == "unhealthy" & metadata$Location == "China"] <- "rheumatoid arthritis"
metadata$Location_Health <- paste(metadata$Location, "-", metadata$Health)

# Body site colours
cols <- brewer.pal(11, "Spectral")[c(2, 4, 9, 10, 11)]
names(cols) <- unique(metadata$sample_type)



########## Read in and process raw data ####
#Contigs
contig_ids = read.table("../../data/virsorter_positive.ids")
contig_ids = unique(as.character(contig_ids$V1))

# Remove jumbophages that are not viral
jumbo_contigs <- read.delim("../../data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
contig_ids <- contig_ids[contig_ids %in% jumbo_contigs$name | as.numeric(sapply(strsplit(contig_ids, split = "_"), "[[", 5)) <= 2e5]

#Read alignment counts
vir_counts = read.table("../../data/bowtie2_read_counts.txt", header=TRUE, row.names = 1, sep = "\t", check.names = FALSE)
#vir_counts = fread("../../data/bowtie2_read_counts.txt", header=TRUE, sep = "\t", check.names = FALSE)
vir_coverage = read.table("../../data/breadth_cov_collated.tbl", header=TRUE, row.names = 1, check.names = FALSE)
#vir_coverage = fread("../../data/breadth_cov_collated.tbl", header=TRUE, sep = "\t", check.names = FALSE)

#Remove unwanted samples
vir_counts <- vir_counts[,names(vir_counts) %in%  metadata$ID]
vir_coverage <- vir_coverage[,names(vir_coverage) %in% metadata$ID]

#Remove unwanted contigs
vir_counts <- vir_counts[rownames(vir_counts) %in% c("Total_reads", contig_ids),]
vir_coverage <- vir_coverage[rownames(vir_coverage) %in% contig_ids,]

counts_total = vir_counts["Total_reads",]
vir_counts = vir_counts[-c(nrow(vir_counts), nrow(vir_counts)-1),]
vir_coverage = vir_coverage[rownames(vir_counts), colnames(vir_counts)]
vir_counts[vir_coverage < 0.75] = 0 #Apply breadth of coverage filter (<75% coverage are removed)
rm(vir_coverage)

vir_counts_unaligned = counts_total - apply(vir_counts, 2, sum)
vir_counts = rbind(vir_counts, vir_counts_unaligned)
vir_counts_prop = apply(vir_counts, 2, function(x) x/sum(x))
rm(vir_counts)
vir_counts_prop = vir_counts_prop[-nrow(vir_counts_prop),]

#Remove samples with missing total counts
vir_counts_prop = vir_counts_prop[,-which(apply(vir_counts_prop, 2, function(x) any(is.na(x))))]

#Remove rows with all zeroes
vir_counts_prop = vir_counts_prop[-which(apply(vir_counts_prop, 1, function(x) all(x == 0))),]

saveRDS(vir_counts_prop, file ="../../data/vir_counts_prop.RDS")




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
  vcontact_cluster = NA,
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
vcontact = read.table("../data/new_clusters_vcontact_w_taxa.txt", header=TRUE, sep = "\t")
rownames(vcontact) = as.character(vcontact$contig_ID)
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
contig_data$vcontact_cluster = vcontact[rownames(contig_data), "Cluster_id"]

country = read.table("../data/sample_countries.txt")
rownames(country) = country$V2
contig_data$sample = sapply(strsplit(as.character(contig_data$name), split = "_"), "[[", 1)

#contig_data[contig_data$sample %in% unlist(merged_samples), "sample"] =
#  as.character(sapply(contig_data[contig_data$sample %in% unlist(merged_samples), "sample"],
#                      function(x) merged_samples[[grep(x, merged_samples)]][1]))

contig_data$country = metadata[contig_data$sample,"Location"]
contig_data[which((is.na(contig_data$country))),"country"] =
  country[contig_data[which((is.na(contig_data$country))),"sample"],"V1"]

#Filter contigs based on number of ribosomal protein hits and pVOGs
contig_data[is.na(contig_data$ribo_prot_count), "ribo_prot_count"] = 0
contig_data[is.na(contig_data$pvogs_count), "pvogs_count"] = 0
#contig_data = contig_data[which(contig_data$ribo_prot_count <= 3),]
#contig_data = contig_data[which(contig_data$pvogs_count/(contig_data$size/1000) >= 0.3),]




########## Host information from CRISPR and tRNA matches ####
#Read in and process CRISPR BLAST data
crispr_bitscore_cutoff = 45
trna_bitscore_cutoff = 65
crispr_evalue_cutoff = 1E-05
trna_evalue_cutoff = 1E-10

#crispr_hits1 =
#  read.table(file = 'refseq89_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)
crispr_hits2 =
  read.table(file = '../data/hmp_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)

#crispr_hits = rbind(crispr_hits1, crispr_hits2)
crispr_hits = crispr_hits2

colnames(crispr_hits) = c("qseqid", "sseqid", "bitscore", "pident", "length", "qlen", "evalue", "sstart", "send", "sname")
crispr_hits = crispr_hits[crispr_hits$evalue < crispr_evalue_cutoff,]
crispr_hits = crispr_hits[crispr_hits$bitscore >= crispr_bitscore_cutoff,]#Remove low score hits
crispr_hits$sname = gsub("UNVERIFIED: ", "", crispr_hits$sname)
crispr_hits$sname = gsub("TPA: ", "", crispr_hits$sname)
crispr_hits$sname = gsub("Candidatus ", "", crispr_hits$sname)
crispr_hits$sname = gsub("^complete chromosome ", "", crispr_hits$sname)
crispr_hits$genus =
  sapply(strsplit(as.character(crispr_hits$sname), split = " "), "[[", 2) #Check manually!

crispr_summary <- crispr_hits %>% group_by(sseqid, genus) %>%
  summarise(sum_bitscore = sum(bitscore)) %>%
  group_by(sseqid) %>%
  filter(sum_bitscore == max(sum_bitscore)) %>%
  select(-sum_bitscore)

contig_data <- left_join(contig_data, crispr_summary, by = c("name"="sseqid")) %>%
  rename(crispr_host = genus)




########## Taxonomy from co-clustering with RefSeq ####
clusters_tax = list()
for(i in unique(as.character(contig_data$vcontact_cluster))) {
  tax_string = as.character(vcontact[which(vcontact$Cluster_id == i),"Taxonomy"])
  if(length(tax_string[!is.na(tax_string)]) > 0) {
    clusters_tax[[i]] = tax_string[!is.na(tax_string)]
  }
}
sink("clusters_taxonomy.txt")
print(clusters_tax)
sink()

#Dereplicate contigs according to the set used in read alignment step
contig_data <- contig_data[!duplicated(contig_data$name),]

#Cleanup
rm(integrase,circular,vir_refseq,crasslike,demovir,gc_content,ribo_prot_count,
   pvogs_count,vcontact,crispr_hits2)
contig_data <- data.table(contig_data)
saveRDS(contig_data,  file = "../data/contig_data.RDS")
write.table(contig_data, file = "../data/contig_data.txt")




########## Contig scatterplot ####
contig_scatter_demovir = ggplot(contig_data, aes(x = size, y = coverage, fill = demovir, shape = circular, size = circular)) +
  theme_classic() +
  geom_point(alpha = 0.7) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + annotation_logticks() +
  scale_fill_manual(values = c(brewer.pal(12, "Paired")[c(1,2,4,5,7,8,9,10,12)], "lightgray"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 3),
                                         keyheight = 0.3, default.unit = "cm", ncol = 2)) +
  scale_size_manual(values = c(1,3)) + scale_shape_manual(values = c(22,21)) +
  theme(legend.position = "bottom") +
  xlab("Size, bp") + ylab("Coverage, x") +
  facet_wrap(~integrase, nrow = 2)




########## Virome composition ####
# Make long format table
vir_counts_prop <- readRDS("../data/vir_counts_prop.RDS")
contig_data <- readRDS("../data/contig_data.RDS")
vir_counts_prop_melt = melt(vir_counts_prop)
rm(vir_counts_prop)
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, contig_data[,c("name", "demovir", "vcontact_cluster", "crispr_host", "integrase")], by = c("Var1"="name"))

# Aggregate counts by Demovir/vConTACT2
vir_counts_prop_melt <- vir_counts_prop_melt %>% filter(value != 0)

# Add metadata
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, metadata[,c("ID", "sample_type", "Location")], by = c("Var2"="ID"))

saveRDS(vir_counts_prop_melt, file = "../data/vir_counts_prop_melt.RDS")