library(data.table)
library(reshape2)
library(dplyr)

# Read data
jumbophage_contigs <- read.table("../data/jumbophage_contigs.txt", stringsAsFactors = FALSE)
metadata = read.csv("../data/metadata_v2.csv", stringsAsFactors = FALSE)

#### HMMER 
filterBySmallestEvalue <- function(alignment_results) {
  alignment_results_filtered <- alignment_results %>% 
    group_by(query_name) %>%
    filter(evalue == min(evalue)) %>%
    ungroup() %>%
    group_by(query_name) %>%
    filter(evalue_domain == min(evalue_domain))
  return(alignment_results_filtered)
}

jumbophages_pvogs <- read.table("../data/jumbophage_pvogs_hmm_tblout.txt", stringsAsFactors = FALSE)
jumbophages_pvogs$contig <- sub("_[^_]+$", "", jumbophages_pvogs$V1)
jumbophages_pvogs <- jumbophages_pvogs[jumbophages_pvogs$contig %in% unique(jumbophage_contigs$name),]
jumbophages_pvogs <- jumbophages_pvogs[,names(jumbophages_pvogs) != "contig"]

names(jumbophages_pvogs) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                              "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                              "clu", "ov", "env", "dom", "rep", "inc")
jumbophages_pvogs <- filterBySmallestEvalue(jumbophages_pvogs)

# Join descriptions
pvogs_metadata <- read.table("../db/AllvogHMMprofiles/VOGTable_singleprot.txt", stringsAsFactors = FALSE, sep = "\t", quote = "")
jumbophages_pvogs <- left_join(jumbophages_pvogs, pvogs_metadata, by = c("target_name"="V1")) %>%
  rename(description = V7) %>%
  select(-c("V2", "V3", "V4", "V5", "V6"))

#### HMMER on Pfam
jumbophages_pfam <- read.table(file="../data/jumbophage_proteins_pfam.out", sep = "", stringsAsFactors = F, header = FALSE, skip = 3)
names(jumbophages_pfam) <- c("query_name", "query_accession", "target_accession", "target_name", "evalue", 
                             "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                             "clu", "ov", "env", "dom", "rep", "inc")
jumbophages_pfam <- filterBySmallestEvalue(jumbophages_pfam)

# Join descriptions
pfam_metadata <- data.frame(ids = unique(jumbophages_pfam$target_name))
description <- rep(NA, nrow(pfam_metadata))
for(i in 1:nrow(pfam_metadata)) {
  tmp <- system(paste0("grep -A2 ", pfam_metadata$ids[i], "$ ../db/PFAM/Pfam-A.hmm.dat"), intern = TRUE)
  description[i] <- gsub("#=GF DE   ", "", tmp[grep("DE   ", tmp)])
}
pfam_metadata$description <- description
jumbophages_pfam <- left_join(jumbophages_pfam, pfam_metadata, by = c("target_name"="ids"))

#### HMMER on TIGRFAM (with no eggnog annotation)
jumbophages_tigrfams <- read.table(file="../data/jumbophage_proteins_tigrfams.out", sep = "", stringsAsFactors = FALSE, header = FALSE)
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
jumbophages_tigrfams <- left_join(jumbophages_tigrfams, tigrfams_metadata, by = c("target_name"="ids")) 

# Combine results and select protein with highest bit score
jumbophages_proteins <- rbind(as.data.frame(jumbophages_pvogs), as.data.frame(jumbophages_pfam), as.data.frame(jumbophages_tigrfams))
jumbophages_proteins <- jumbophages_proteins %>% 
  group_by(query_name) %>%
  filter(score == max(score)) %>%
  filter(score_domain == max(score_domain)) %>%
  ungroup()

# Remove the few remaining duplicates
jumbophages_proteins <- jumbophages_proteins[!duplicated(jumbophages_proteins$query_name),]

jumbophages_proteins <- jumbophages_proteins %>% select(query_name, target_name, description)
header_names <- c("coding_region", "protein", "protein_description")
names(jumbophages_proteins) <- header_names

# Extract contig id
jumbophages_proteins$contig <- sub("_[^_]+$", "", jumbophages_proteins$coding_region)

# Modify prodigal GFF to include protein descriptions
prodigal_gff <- read.table("../data/jumbophage_contigs_prodigal.gff", stringsAsFactors = FALSE)
prodigal_gff$coding_region <- paste0(prodigal_gff$V1, "_", gsub(".*_", "", gsub(";.*", "", prodigal_gff$V9)))
prodigal_gff <- left_join(prodigal_gff, jumbophages_proteins, by = "coding_region")

# Read protein family lookup
prot_lookup <- read.table("../db/ONTOLOGIES/protein_family_lookup.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote="\"")
prot_lookup <- prot_lookup[order(prot_lookup$description_contains),] # Order so more unique lookups get selected

prodigal_gff$Name <- rep(NA, nrow(prodigal_gff))
prodigal_gff$Parent <- rep(NA, nrow(prodigal_gff))
for (i in 1:nrow(prot_lookup)) {
  tmp_inds <- grep(prot_lookup$description_contains[i], prodigal_gff$protein_description)
  prodigal_gff$Name[tmp_inds] <- prot_lookup$name[i]
  prodigal_gff$Parent[tmp_inds] <- prot_lookup$parent[i]
}
write.csv(prodigal_gff, "../data/jumbophage_prodigal.gff", row.names = FALSE)

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
phages_w_cp <- unique(prodigal_gff$V1[prodigal_gff$Name %in% "Major capsid protein"])
linear_phages <- jumbophage_contigs$name[!jumbophage_contigs$circular]
jumbophage_contigs <- jumbophage_contigs[jumbophage_contigs$name %in% unique(c(phages_w_cp, linear_phages)),]
write.table(jumbophage_contigs, file = "../data/jumbophage_contigs_filtered.txt")

# Output jumbo phage headers
write.table(jumbophage_contigs$name, file = "../data/jumbophage_contigs_for_annot.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Remove filtered jumbophages from catalogue
contig_data <- readRDS("../data/contig_data.RDS")
contig_data_filtered <- contig_data[contig_data$name %in% c(jumbophage_contigs$name, contig_data$name[contig_data$size < 200000]),]
saveRDS(contig_data_filtered, file = "../data/contig_data_filtered.RDS")

# Jumbophage headers with capsid proteins
jumbophages_proteins$protein_description <- tolower(jumbophages_proteins$protein_description)
jumbophage_mcps <- jumbophages_proteins[grepl("major capsid protein", jumbophages_proteins$protein_description),]
jumbophage_contig_meta <- inner_join(metadata, jumbophage_contigs, by = c("ID"="sample"))
jumbophage_mcps <- jumbophage_mcps[jumbophage_mcps$contig %in% jumbophage_contig_meta$name,]
write.table(jumbophage_mcps$coding_region, "../data/jumbophage_mcp_headers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

########## Get MCP annotations for IMG/VR ####
imgvr_pvogs <- read.table("../data/imgvr_uvigs_pvogs_hmm_tblout.txt", stringsAsFactors = FALSE)
imgvr_pvogs <- imgvr_pvogs[,-ncol(imgvr_pvogs)]
names(imgvr_pvogs) <- c("target_name", "target_accession", "query_name", "query_accession", "evalue", 
                        "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                        "clu", "ov", "env", "dom", "rep", "inc")
imgvr_pvogs <- filterBySmallestEvalue(imgvr_pvogs)

# Join descriptions
imgvr_pvogs <- left_join(imgvr_pvogs, pvogs_metadata, by = c("target_name"="V1")) %>%
  rename(description = V7) %>%
  select(-c("V2", "V3", "V4", "V5", "V6"))

#### HMMER on Pfam
imgvr_pfam <- read.table(file="../data/imgvr_uvigs_pfam.out", sep = "", stringsAsFactors = F, header = FALSE, skip = 3)
names(imgvr_pfam) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                       "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                       "clu", "ov", "env", "dom", "rep", "inc")
imgvr_pfam <- filterBySmallestEvalue(imgvr_pfam)

# Join descriptions
pfam_metadata <- data.frame(ids = unique(imgvr_pfam$target_accession))
description <- rep(NA, nrow(pfam_metadata))
for(i in 1:nrow(pfam_metadata)) {
  tmp <- system(paste0("grep -A2 ", pfam_metadata$ids[i], "$ ../db/PFAM/Pfam-A.hmm.dat"), intern = TRUE)
  description[i] <- gsub("#=GF DE   ", "", tmp[grep("DE   ", tmp)])
}
pfam_metadata$description <- description
imgvr_pfam <- left_join(imgvr_pfam, pfam_metadata, by = c("target_accession"="ids"))

#### HMMER on TIGRFAM (with no eggnog annotation)
imgvr_tigrfams <- read.table(file="../data/imgvr_uvigs_tigrfams.out", sep = "", stringsAsFactors = FALSE, header = FALSE)
names(imgvr_tigrfams) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                           "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                           "clu", "ov", "env", "dom", "rep", "inc")
imgvr_tigrfams <- filterBySmallestEvalue(imgvr_tigrfams)

# Join descriptions
tigrfams_metadata <- data.frame(ids = unique(imgvr_tigrfams$target_name))
description <- rep(NA, nrow(tigrfams_metadata))
for(i in 1:nrow(tigrfams_metadata)) {
  tmp <- read.table(file = paste0("../db/TIGRFAMS/", tigrfams_metadata$ids[i], ".INFO"), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  description[i] <- gsub("DE  ", "", tmp$V1[grep("DE  ", tmp$V1)])
}
tigrfams_metadata$description <- description
imgvr_tigrfams <- left_join(imgvr_tigrfams, tigrfams_metadata, by = c("target_name"="ids"))

# Combine results and select protein with highest bit score
imgvr_proteins <- rbind(as.data.frame(imgvr_pvogs), as.data.frame(imgvr_pfam), as.data.frame(imgvr_tigrfams))
imgvr_proteins <- imgvr_proteins %>% 
  group_by(query_name) %>%
  filter(score == max(score)) %>%
  filter(score_domain == max(score_domain)) %>%
  ungroup()

# Remove the few remaining duplicates
imgvr_proteins <- imgvr_proteins[!duplicated(imgvr_proteins$query_name),]

imgvr_proteins <- imgvr_proteins %>% select(query_name, target_name, description)
header_names <- c("coding_region", "protein", "protein_description")
names(imgvr_proteins) <- header_names

# Get MCPs
imgvr_proteins$protein_description <- tolower(imgvr_proteins$protein_description)
imgvr_mcps <- imgvr_proteins[grepl("major capsid protein", imgvr_proteins$protein_description),]
write.table(imgvr_mcps$coding_region, "../data/imgvr_mcp_headers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
