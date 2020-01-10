library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(igraph)

########## Read in metadata and define parameters ####
metadata_contigs <- read.csv("../data/metadata_andrey.csv", stringsAsFactors = FALSE)
metadata = read.csv("../data/metadata_v2.csv", stringsAsFactors = FALSE)

# Body site colours
cols <- brewer.pal(11, "Spectral")[c(2, 4, 9, 10, 11)]
names(cols) <- unique(metadata$sample_type)

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

country = read.table("../data/sample_countries.txt")
rownames(country) = country$V2
contig_data$sample = sapply(strsplit(as.character(contig_data$name), split = "_"), "[[", 1)

#contig_data[contig_data$sample %in% unlist(merged_samples), "sample"] =
#  as.character(sapply(contig_data[contig_data$sample %in% unlist(merged_samples), "sample"],
#                      function(x) merged_samples[[grep(x, merged_samples)]][1]))

contig_data$country = country[contig_data$sample,"V1"]
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
crispr_host_species <- gsub("\\]", "", gsub(".*\\[", "", crispr_hits$sname))
crispr_hits$crispr_host <- gsub(" .*", "", crispr_host_species)

crispr_summary <- crispr_hits %>% group_by(sseqid, crispr_host) %>%
  summarise(sum_bitscore = sum(bitscore)) %>%
  group_by(sseqid) %>%
  filter(sum_bitscore == max(sum_bitscore)) %>%
  select(-sum_bitscore)

contig_data <- left_join(contig_data, crispr_summary, by = c("name"="sseqid"))




########## Taxonomy from co-clustering with RefSeq ####
# clusters_tax = list()
# for(i in unique(as.character(contig_data$vcontact_cluster))) {
#   tax_string = as.character(vcontact[which(vcontact$VC.Subcluster == i),"Taxonomy"])
#   if(length(tax_string[!is.na(tax_string)]) > 0) {
#     clusters_tax[[i]] = tax_string[!is.na(tax_string)]
#   }
# }
# sink("clusters_taxonomy.txt")
# print(clusters_tax)
# sink()

#Dereplicate contigs according to the set used in read alignment step
contig_data <- contig_data[!duplicated(contig_data$name),]

#Cleanup
rm(integrase,circular,vir_refseq,crasslike,demovir,gc_content,ribo_prot_count,
   pvogs_count,vcontact,crispr_hits2)
contig_data <- data.table(contig_data)

saveRDS(contig_data,  file = "../data/contig_data.RDS")
write.table(contig_data, file = "../data/contig_data.txt")

# ########## Contig scatterplot
# contig_scatter_demovir = ggplot(contig_data, aes(x = size, y = coverage, fill = demovir, shape = circular, size = circular)) +
#   theme_classic() +
#   geom_point(alpha = 0.7) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) + annotation_logticks() +
#   scale_fill_manual(values = c(brewer.pal(12, "Paired")[c(1,2,4,5,7,8,9,10,12)], "lightgray"),
#                     guide = guide_legend(override.aes = list(shape = 21, size = 3),
#                                          keyheight = 0.3, default.unit = "cm", ncol = 2)) +
#   scale_size_manual(values = c(1,3)) + scale_shape_manual(values = c(22,21)) +
#   theme(legend.position = "bottom") +
#   xlab("Size, bp") + ylab("Coverage, x") +
#   facet_wrap(~integrase, nrow = 2)

# Remove spurious jumbophages
#Contigs over 200kb
jumbophage_contigs = contig_data[which(contig_data$size >= 200000),]
#Circular genomes
jumbophage_contigs_circ = jumbophage_contigs[which(jumbophage_contigs$circular),]
#Extract sub-graph
vcontact_sig = fread("../data/cc_sig1.0_mcl1.5.ntw")
vcontact_sig_jumboph = vcontact_sig[V1 %in% jumbophage_contigs$name & V2 %in% jumbophage_contigs$name]
#Only contigs connected with a circular genome
vcontact_sig_jumboph = vcontact_sig_jumboph[V1 %in% jumbophage_contigs_circ$name | V2 %in% jumbophage_contigs_circ$name]
#Remove vertices connected with < 2 edges
vcontact_sig_jumboph = vcontact_sig_jumboph[V1 %in% names(which(table(vcontact_sig_jumboph$V1)>1)) & V2 %in% names(which(table(vcontact_sig_jumboph$V1)>1))]
#Update jumbophage contigs dataframe
jumbophage_contigs = jumbophage_contigs[jumbophage_contigs$name %in% unique(vcontact_sig_jumboph$V1),]
rm(jumbophage_contigs_circ)
# vcontact_igraph = igraph::graph.data.frame(vcontact_sig_jumboph, directed = FALSE, vertices = jumbophage_contigs)
# vcontact_igraph = simplify(vcontact_igraph)
# graph_layout = layout_with_fr(vcontact_igraph)
# 
# nclusters = length(levels(as.factor(V(vcontact_igraph)$vcontact_cluster)))
# rand_colours =
#   colorRampPalette(c(brewer.pal(8, "Spectral"),brewer.pal(8, "RdBu"), brewer.pal(8, "Dark2")))(nclusters)
# rand_colours = rand_colours[sample(1:nclusters, nclusters)]
# rand_colours = rand_colours[as.numeric(as.factor(V(vcontact_igraph)$vcontact_cluster))]
# 
# labels_clusters = paste(as.character(sapply(clusters_tax[V(vcontact_igraph)$vcontact_cluster], "[[", 1)),
#                         V(vcontact_igraph)$crispr_host)
# labels_clusters = gsub("NULL", "", labels_clusters)
# labels_clusters = gsub("NA", "", labels_clusters)
# labels_clusters[labels_clusters != " "] = paste(V(vcontact_igraph)$vcontact_cluster[labels_clusters != " "], labels_clusters[labels_clusters != " "])
# labels_clusters = gsub("  ", " ", labels_clusters)

########## Functional annotation of jumbophage candidates 

#### DIAMOND on NCBI NR
jumbophages_ncbi <- read.table(file="../data/megaphage_proteins_diamond.out", sep = "\t", stringsAsFactors = FALSE)
names(jumbophages_ncbi) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore") 

# Filter by best e-value
jumbophages_ncbi <- jumbophages_ncbi %>% group_by(qseqid) %>%
  filter(evalue == min(evalue))

# Read uniprotkb
uniprot_metadata <- read.table(file="../db/UNIPROTKB/uniprotkb_refseq.tab", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Remove duplicates
uniprot_metadata_char <- uniprot_metadata[uniprot_metadata$Protein.names != "Uncharacterized protein",]
uniprot_metadata_char <- uniprot_metadata_char[!duplicated(uniprot_metadata_char$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X),]
uniprot_metadata_unchar <- uniprot_metadata[uniprot_metadata$Protein.names == "Uncharacterized protein",]
uniprot_metadata_unchar <- uniprot_metadata_unchar[!uniprot_metadata_unchar$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X %in% uniprot_metadata_char$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X,]
uniprot_metadata <- rbind(uniprot_metadata_char, uniprot_metadata_unchar)

jumbophages_ncbi <- left_join(jumbophages_ncbi, uniprot_metadata, by = c("sseqid" = "yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X"))
jumbophages_ncbi <- jumbophages_ncbi[!is.na(jumbophages_ncbi$Entry),]
jumbophages_ncbi <- jumbophages_ncbi %>% select(qseqid, sseqid, Protein.names)

# Remove duplicates
jumbophages_ncbi <- jumbophages_ncbi[!duplicated(jumbophages_ncbi$qseqid),]

#### HMMER on pVOGs
filterBySmallestEvalue <- function(alignment_results) {
  alignment_results_filtered <- alignment_results %>% 
    group_by(query_name) %>%
    filter(evalue == min(evalue)) %>%
    ungroup()
  return(alignment_results_filtered)
}

jumbophages_pvogs <- read.table("../data/megaphage_pvogs_hmm_tblout.txt", stringsAsFactors = FALSE)

names(jumbophages_pvogs) <- c("target_name", "target_accession", "query_name", "query_accession", "evalue", 
                             "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                             "clu", "ov", "env", "dom", "rep", "inc", "description")
jumbophages_pvogs <- filterBySmallestEvalue(jumbophages_pvogs)

# Join descriptions
pvogs_metadata <- read.table("../db/AllvogHMMprofiles/VOGTable_singleprot.txt", stringsAsFactors = FALSE, sep = "\t", quote = "")
jumbophages_pvogs <- left_join(jumbophages_pvogs, pvogs_metadata, by = c("target_name"="V1")) %>%
  select(query_name, target_name, description = V7)

#### HMMER on Pfam
jumbophages_pfam <- read.table(file="../data/megaphage_proteins_pfam.out", sep = "", stringsAsFactors = F, header = FALSE, skip = 3)
names(jumbophages_pfam) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                            "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                            "clu", "ov", "env", "dom", "rep", "inc")
jumbophages_pfam <- filterBySmallestEvalue(jumbophages_pfam)

# Join descriptions
pfam_metadata <- data.frame(ids = unique(jumbophages_pfam$target_accession))
description <- rep(NA, nrow(pfam_metadata))
for(i in 1:nrow(pfam_metadata)) {
  tmp <- system(paste0("grep -A2 ", pfam_metadata$ids[i], "$ ../db/PFAM/Pfam-A.hmm.dat"), intern = TRUE)
  description[i] <- gsub("#=GF DE   ", "", tmp[grep("DE   ", tmp)])
}
pfam_metadata$description <- description
jumbophages_pfam <- left_join(jumbophages_pfam, pfam_metadata, by = c("target_accession"="ids")) %>%
  select(query_name, target_accession, description)

#### HMMER on TIGRFAM (with no eggnog annotation)
jumbophages_tigrfams <- read.table(file="../data/megaphage_proteins_tigrfams.out", sep = "", stringsAsFactors = FALSE, header = FALSE)
names(jumbophages_tigrfams) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                                "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                                "clu", "ov", "env", "dom", "rep", "inc")
jumbophages_tigrfams <- filterBySmallestEvalue(jumbophages_tigrfams)

# Join descriptions
tigrfams_metadata <- data.frame(ids = unique(jumbophages_tigrfams$target_name))
description <- rep(NA, nrow(tigrfams_metadata))
for(i in 1:nrow(tigrfams_metadata)) {
  tmp <- read.table(file = paste0("../db/TIGRFAMS/", tigrfams_metadata$ids[i], ".INFO"), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  description[i] <- gsub("DE  ", "", tmp$V1[grep("DE  ", tmp$V1)])
}
tigrfams_metadata$description <- description
jumbophages_tigrfams <- left_join(jumbophages_tigrfams, tigrfams_metadata, by = c("target_name"="ids")) %>%
  select(query_name, target_name, description)

# #### HHBLIT on largest phages
# hhblit_result <- data.frame()
# for (i in 1:length(largest_phages)) {
#   hhblit_result_tmp <- read.table(paste0("../data/proteins/", largest_phages[i], ".1E-5.tab"), stringsAsFactors = FALSE)
#   hhblit_result <- rbind(hhblit_result, hhblit_result_tmp)
# }
# names(hhblit_result) <- c("query", "target", "match", "tlen", "mismatch", "gap_open", "qstart", "qend", "tstart", "tend", "eval", "score")
# hhblit_result_filter <- hhblit_result %>% group_by(query) %>%
#   filter(eval == min(eval))
# 
# uni_lookup <- fread("../db/UNICLUST30/uniclust_uniprot_mapping.tsv", stringsAsFactors = FALSE)
# uniprot_id <- rep(NA, nrow(hhblit_result_filter))
# for (i in 1:length(uniprot_id)) {
#   uniprot_id[i] <- uni_lookup$V2[hhblit_result_filter$target[i] + 1]
# }
# write.table(uniprot_id, "../data/largest_phages_uniprot_ids.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# uniprot_id <- read.delim("../data/largest_phages_uniprot_ids.txt", stringsAsFactors = FALSE, header = FALSE)
# uniprot_lookup <- read.delim("../db/UNIPROTKB/uniprot_uniclust_largest_phages_lookup.tab", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# hhblit_result_filter$Entry <- uniprot_id$V1
# jumbophages_hhblit <- left_join(hhblit_result_filter, uniprot_lookup, by = "Entry") %>%
#   select(query, Entry, Protein.names)

#### Combine results (in order ncbi, pvogs, pfam, tigrfam)
header_names <- c("coding_region", "protein", "protein_description")
names(jumbophages_ncbi) <- header_names
names(jumbophages_pvogs) <- header_names
names(jumbophages_pfam) <- header_names
names(jumbophages_tigrfams) <- header_names
# names(jumbophages_hhblit) <- header_names
jumbophages_proteins <- rbind(as.data.frame(jumbophages_ncbi[jumbophages_ncbi$protein_description != "Uncharacterized protein",]),
                             as.data.frame(jumbophages_pvogs), as.data.frame(jumbophages_pfam), as.data.frame(jumbophages_tigrfams)
                             # ,as.data.frame(jumbophages_hhblit)
)

# Check duplications sum(duplicated(jumbophages_proteins$contig))
jumbophages_proteins <- jumbophages_proteins[!duplicated(jumbophages_proteins$coding_region),]

# Extract contig id
jumbophages_proteins$contig <- sub("_[^_]+$", "", jumbophages_proteins$coding_region)

# Modify prodigal GFF to include protein descriptions
prodigal_gff <- read.table("../data/megaphage_contigs_prodigal.gff", stringsAsFactors = FALSE)
prodigal_gff$coding_region <- paste0(prodigal_gff$V1, "_", gsub(".*_", "", gsub(";.*", "", prodigal_gff$V9)))
prodigal_gff <- left_join(prodigal_gff, jumbophages_proteins, by = "coding_region")

# Read protein family lookup
prot_lookup <- read.table("../db/ONTOLOGIES/protein_family_lookup.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
prot_lookup <- prot_lookup[order(prot_lookup$description_contains),] # Order so more unique lookups get selected

prodigal_gff$Name <- rep(NA, nrow(prodigal_gff))
prodigal_gff$Parent <- rep(NA, nrow(prodigal_gff))
for (i in 1:nrow(prot_lookup)) {
  tmp_inds <- grep(prot_lookup$description_contains[i], prodigal_gff$protein_description)
  prodigal_gff$Name[tmp_inds] <- prot_lookup$name[i]
  prodigal_gff$Parent[tmp_inds] <- prot_lookup$parent[i]
}

# Count how many structural proteins in each phage
struct_prot_count <- prodigal_gff %>% group_by(V1, Parent) %>%
  summarise(n = n()) %>%
  filter(Parent == "Structural proteins")

# Calculate prevalence of different types of structural proteins
struct_prot_prev <- prodigal_gff %>% 
  filter(Parent == "Structural proteins") %>% 
  group_by(Name) %>%
  summarise(n = n())

# Remove jumbophages without major capsid protein
jumbophage_contigs <- jumbophage_contigs[jumbophage_contigs$name %in% unique(prodigal_gff$V1[prodigal_gff$Name == "Major capsid protein"]),]
write.table(jumbophage_contigs, file = "../data/jumbophage_contigs.txt")

# # Add Name and Parent to attributes
# prodigal_gff$V9 <- ifelse(is.na(prodigal_gff$Name), 
#                           paste0(prodigal_gff$V9, "Name=hypothetical protein;"),
#                           paste0(prodigal_gff$V9, "Name=", gsub(" \\[.*", "", prodigal_gff$Name), ";"))
# prodigal_gff$V9 <- ifelse(is.na(prodigal_gff$Parent), 
#                           paste0(prodigal_gff$V9, "Parent=None;"),
#                           paste0(prodigal_gff$V9, "Parent=", prodigal_gff$Parent, ";"))
# prodigal_gff <- prodigal_gff[,c(1:9)]
# names(prodigal_gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# 
# # Add tRNA
# jumbophages_trna <- read.table("../data/jumbophage_trna_clean.tab", stringsAsFactors = FALSE)
# jumbophages_trna$strand <- ifelse(grepl("c", jumbophages_trna$V3), "-", "+")
# jumbophages_trna <- cbind(jumbophages_trna, map_df(.x = gsub("c", "", gsub("\\]", "", gsub("\\[", "", jumbophages_trna$V3))), 
#                                                  function(.x) {
#                                                    split <- strsplit(.x, ",")
#                                                    return(data.frame(start = as.numeric(split[[1]][1]), 
#                                                                      end = as.numeric(split[[1]][2])))
#                                                  }))
# jumbophages_trna$attributes <- paste0("Name=", jumbophages_trna$V2, ";")
# trna_gff <- jumbophages_trna %>% mutate(source = "ARAGORN_v1.2.36", type = "tRNA", phase = 0, score = format(as.numeric(V4), nsmall=1)) %>%
#   select(seqid = V6, source, type, start, end, score, strand, phase, attributes)
# 
# # Combine prodigal CDS and trna
# comb_gff <- rbind(prodigal_gff, trna_gff)
# 
# # Write GFF file for largest phages
# largest_phages <- jumbophage_contigs_meta$name[order(jumbophage_contigs_meta$size, decreasing = TRUE)][c(1:3)]
# for (i in 1:length(largest_phages)) {
#   largest_phage_gff <- comb_gff[comb_gff$seqid == largest_phages[i],]
#   write.table(largest_phage_gff, paste0("../data/", largest_phages[i], ".gff"), 
#               quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")                          
# }


########## Create read data ####
#Read alignment counts
vir_counts = read.table("../data/bowtie2_read_counts.txt", header=TRUE, row.names = 1, sep = "\t", check.names = FALSE)
#vir_counts = fread("../data/bowtie2_read_counts.txt", header=TRUE, sep = "\t", check.names = FALSE)
vir_coverage = read.table("../data/breadth_cov_collated.tbl", header=TRUE, row.names = 1, check.names = FALSE)
#vir_coverage = fread("../data/breadth_cov_collated.tbl", header=TRUE, sep = "\t", check.names = FALSE)

#Remove unwanted samples
vir_counts <- vir_counts[,names(vir_counts) %in%  metadata$ID]
vir_coverage <- vir_coverage[,names(vir_coverage) %in% metadata$ID]

#Remove unwanted contigs
vir_counts <- vir_counts[rownames(vir_counts) %in% c("Total_reads", contig_ids),]
vir_coverage <- vir_coverage[rownames(vir_coverage) %in% contig_ids,]

# Remove jumbophages that are not viral from vir_counts and vir_coverage
vir_counts <- vir_counts[rownames(vir_counts) %in% c("Total_reads", jumbophage_contigs$name, contig_data$name[contig_data$size < 200000]),]
vir_coverage <- vir_coverage[rownames(vir_coverage) %in% c(jumbophage_contigs$name, contig_data$name[contig_data$size < 200000]),]

counts_total = vir_counts["Total_reads",]
saveRDS(counts_total, "../data/counts_total.RDS")

vir_counts = vir_counts[-nrow(vir_counts),]
vir_coverage = vir_coverage[rownames(vir_counts), colnames(vir_counts)]
vir_counts[vir_coverage < 0.75] = 0 #Apply breadth of coverage filter (<75% coverage are removed)
# rm(vir_coverage)

vir_counts_unaligned = counts_total - apply(vir_counts, 2, sum)
vir_counts = rbind(vir_counts, vir_counts_unaligned)
vir_counts_prop = apply(vir_counts, 2, function(x) x/sum(x))
vir_counts_prop = vir_counts_prop[-nrow(vir_counts_prop),]
# rm(vir_counts)

#Remove samples with missing total counts
vir_counts_prop = vir_counts_prop[,!is.na(counts_total)]

#Remove rows with all zeroes
vir_counts_prop = vir_counts_prop[-which(apply(vir_counts_prop, 1, function(x) all(x == 0))),]

saveRDS(vir_counts_prop, file ="../data/vir_counts_prop.RDS")



########## Virome composition ####
# Make long format table
vir_counts_prop <- readRDS("../data/vir_counts_prop.RDS")
contig_data <- readRDS("../data/contig_data.RDS")
vir_counts_prop_melt = melt(vir_counts_prop)
# rm(vir_counts_prop)
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, contig_data[,c("name", "demovir", "vcontact_cluster", "crispr_host", "integrase")], by = c("Var1"="name"))
vir_counts_prop_melt <- vir_counts_prop_melt %>% filter(value != 0)

# Add metadata
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, metadata[,c("ID", "sample_type", "Location")], by = c("Var2"="ID"))
saveRDS(vir_counts_prop_melt, file = "../data/vir_counts_prop_melt.RDS")

# Aggregate counts by Demovir/vConTACT2
vir_counts_prop_melt <- as.data.table(vir_counts_prop_melt)
vir_counts_prop_melt_agg = vir_counts_prop_melt[, sum(value), by=.(Var2, demovir, vcontact_cluster, sample_type, Location)]
saveRDS(vir_counts_prop_melt_agg,  file = "../data/vir_counts_prop_melt_agg.RDS")

# Aggregate counts by CRISPR host
vir_counts_prop_melt_agg2 = vir_counts_prop_melt[, sum(value), by=.(Var2, crispr_host, sample_type, Location)]
saveRDS(vir_counts_prop_melt_agg2, file = "../data/vir_counts_prop_melt_agg2.RDS")

# Re-cast counts matrix by vConTACT2 clusters
vir_counts_prop_agg = dcast.data.table(vir_counts_prop_melt_agg, vcontact_cluster ~ Var2, value.var = "V1", fun.aggregate = sum)
vir_counts_prop_agg = as.data.frame(vir_counts_prop_agg)
vir_counts_prop_agg = vir_counts_prop_agg[-which(is.na(vir_counts_prop_agg[,1])),]
rownames(vir_counts_prop_agg) = vir_counts_prop_agg[,1]
vir_counts_prop_agg = vir_counts_prop_agg[,-1]
vir_counts_prop_agg = as.matrix(vir_counts_prop_agg)
saveRDS(vir_counts_prop_agg,  file = "../data/vir_counts_prop_agg.RDS")

# Re-cast counts matrix by phages
vir_counts_prop_agg_phage <- dcast.data.table(vir_counts_prop_melt, Var1 ~ Var2, value.var = "value", fun.aggregate = sum)
rownames(vir_counts_prop_agg_phage) <- unlist(vir_counts_prop_agg_phage[,1])
vir_counts_prop_agg_phage <- vir_counts_prop_agg_phage[,-1]
vir_counts_prop_agg_phage <- as.matrix(vir_counts_prop_agg_phage)
saveRDS(vir_counts_prop_agg_phage, file = "../data/vir_counts_prop_agg_phage.RDS")
