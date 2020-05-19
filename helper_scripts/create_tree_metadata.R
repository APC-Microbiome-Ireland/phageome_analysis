library(data.table)
library(dplyr)
library(openxlsx)
library(viridis)
library(RColorBrewer)

# Read data
jumbophage_contigs <- read.table("../data/jumbophage_contigs.txt", stringsAsFactors = FALSE)
metadata = read.csv("../data/metadata_v2.csv", stringsAsFactors = FALSE)
jumbophage_mcps <- read.table("../data/jumbophage_mcp_headers.txt", stringsAsFactors = FALSE)
jumbophage_mcps$contig <- sub("_[^_]+$", "", jumbophage_mcps$V1)

# Body site colours
cols <- plasma(length(unique(metadata$sample_type)), end = 0.8)
names(cols) <- sort(unique(metadata$sample_type))

########## Create jumbophage MCP metadata ####
# Our study
jumbophage_contig_meta <- inner_join(jumbophage_contigs, metadata, by = c("sample"="ID")) 
jumbophage_contig_meta$environment <- NA
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "China" & jumbophage_contig_meta$sample_type == "stool"] <- "Adult fecal samples, China"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "China" & jumbophage_contig_meta$sample_type == "saliva"] <- "Adult saliva, China"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "China" & jumbophage_contig_meta$sample_type == "dental"] <- "Adult dental plaque, China"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "US" & jumbophage_contig_meta$sample_type == "dental"] <- "Adult dental plaque, USA"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "US" & jumbophage_contig_meta$sample_type == "buccal mucosa"] <- "Adult buccal mucosa, USA"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "US" & jumbophage_contig_meta$sample_type == "dorsum of tongue"] <- "Adult dorsum of tongue, USA"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "US" & jumbophage_contig_meta$sample_type == "stool"] <- "Adult stool, USA"
jumbophage_contig_meta$environment[jumbophage_contig_meta$Location == "Philippines"] <- "Adult saliva, Philippines"

jumbo_ours_meta <- jumbophage_contig_meta %>% select(name, size, circular, crispr_host, environment) %>%
  mutate(circular = ifelse(circular, "circular", "linear")) %>%
  rename("host" = "crispr_host") %>%
  filter(name %in% jumbophage_mcps$contig)

jumbo_ours_meta$host[jumbo_ours_meta$host %in% c("Veillonella", "Selenomonas")] <- "Firmicutes"
jumbo_ours_meta$host[jumbo_ours_meta$host %in% "Neisseria"] <- "Proteobacteria"
jumbo_ours_meta$host[jumbo_ours_meta$host %in% "Prevotella"] <- "Bacteroidetes"

# Al-Shayeb
jumbo_alshayeb_meta <- read.xlsx("../data/alshayeb/41586_2020_2007_MOESM4_ESM.xlsx", sheet = 1)
jumbo_alshayeb_meta4 <- read.xlsx("../data/alshayeb/41586_2020_2007_MOESM4_ESM.xlsx", sheet = 4)
jumbo_alshayeb_meta <- jumbo_alshayeb_meta %>%
  filter(Genome.name %in% gsub(".gbk", "", list.files("../data/alshayeb/major_capsid_protein_phages/"))) %>%
  left_join(jumbo_alshayeb_meta4[c("Genome", "1st.phylum")], by = c("Genome.name"="Genome")) %>%
  rename("phylum" = "1st.phylum") %>%
  mutate(circular = ifelse(is.na(Circularized), "linear", "circular")) %>%
  select(Genome.name, Length, circular, phylum, Environment.of.origin) %>%
  rename("name" = "Genome.name", "size"="Length", "host"="phylum", "environment"="Environment.of.origin")

# IMG/VR
jumbo_imgvr_meta <- read.csv("../data/imgvr_uvigs_info_200414_edit.csv", stringsAsFactors = FALSE) %>%
  mutate(host = NA, circular = "linear") %>%
  select(Scaffold.ID, Seq.Length, circular, host, Environment) %>%
  rename("name"="Scaffold.ID", "size"="Seq.Length", "environment"="Environment")

# Refseq
jumbo_refseq_headers <- read.delim("../data/refseq_r99/refseq_r99_jumbophage_headers.txt", header = FALSE)
jumbo_refseq_meta <- read.table("../data/refseq_r99/refseq_r99_jumbophage_info_edit.txt", sep = "\t", header = TRUE) %>%
  select(name, length, circular, phylum, environment) %>%
  rename("host"="phylum", "size"="length") %>%
  filter(name %in% jumbo_refseq_headers$V1)

# Combine headers
system("grep '^>' ../data/major_capsid_proteins_cdhit_1 | sed 's/>//g' > ../data/major_capsid_proteins_headers.txt")
mcp_headers <- read.delim("../data/major_capsid_proteins_headers.txt", header = FALSE)
mcp_headers$contigs <- sub("_[^_]+$", "", mcp_headers$V1)
jumbo_together_meta <- rbind(jumbo_ours_meta, jumbo_alshayeb_meta, jumbo_imgvr_meta, jumbo_refseq_meta)
jumbo_together_meta <- left_join(mcp_headers, jumbo_together_meta, by = c("contigs"="name")) %>%
  select(-contigs) %>%
  rename("name"="V1") %>%
  rename("location"="environment") 

# Add cleaned metadata
jumbo_meta_lookup <- read.delim("../db/ONTOLOGIES/jumbophage_metadata_lookup.txt", stringsAsFactors = FALSE)
jumbo_mcp_meta_edit <- left_join(jumbo_together_meta, jumbo_meta_lookup, by = "name")

# Environment
environment_cols <- rep(NA, length(unique(jumbo_mcp_meta_edit$environment)))
names(environment_cols) <- unique(jumbo_mcp_meta_edit$environment)
environment_cols[is.na(names(environment_cols))] <- "#FFFFFF"
environment_cols[["Adult buccal mucosa"]] <- cols[["buccal mucosa"]]
environment_cols[["Adult dorsum of tongue"]] <- cols[["dorsum of tongue"]]
environment_cols[["Adult saliva"]] <- cols[["saliva"]]
environment_cols[["Adult dental plaque"]] <- cols[["dental"]]
environment_cols[["Adult fecal samples"]] <- cols[["stool"]]
environment_cols[is.na(environment_cols)] <- brewer.pal(sum(is.na(environment_cols)), "Set2")

environment_cols <- environment_cols[order(names(environment_cols))]
jumbo_mcp_meta_edit$environment_cols <- sapply(jumbo_mcp_meta_edit$environment, function(x) environment_cols[x])
jumbo_mcp_meta_edit$environment_cols[is.na(jumbo_mcp_meta_edit$environment_cols)] <- "#FFFFFF"

system("rm ../data/itol/environment_legend.txt")
system(paste0("echo 'LEGEND_SHAPES\t", paste0(rep("1", length(environment_cols[!is.na(names(environment_cols))])), collapse = "\t"), "' >> ../data/itol/environment_legend.txt"))
system(paste0("echo 'LEGEND_COLORS\t", paste0(environment_cols[!is.na(names(environment_cols))], collapse = "\t"), "' >> ../data/itol/environment_legend.txt"))
system(paste0("echo 'LEGEND_LABELS\t", paste0(names(environment_cols[!is.na(names(environment_cols))]), collapse = "\t"), "' >> ../data/itol/environment_legend.txt"))

# Hosts
host_cols <- brewer.pal(length(unique(jumbo_mcp_meta_edit$host)), "Set3")
names(host_cols) <- unique(jumbo_mcp_meta_edit$host)
host_cols[is.na(names(host_cols))] <- "#FFFFFF"
host_cols <- host_cols[order(names(host_cols))]
jumbo_mcp_meta_edit$host_cols <- sapply(jumbo_mcp_meta_edit$host, function(x) host_cols[x])
jumbo_mcp_meta_edit$host_cols[is.na(jumbo_mcp_meta_edit$host_cols)] <- "#FFFFFF"

system("rm ../data/itol/host_legend.txt")
system(paste0("echo 'LEGEND_SHAPES\t", paste0(rep("1", length(host_cols[!is.na(names(host_cols))])), collapse = "\t"), "' >> ../data/itol/host_legend.txt"))
system(paste0("echo 'LEGEND_COLORS\t", paste0(host_cols[!is.na(names(host_cols))], collapse = "\t"), "' >> ../data/itol/host_legend.txt"))
system(paste0("echo 'LEGEND_LABELS\t", paste0(names(host_cols[!is.na(names(host_cols))]), collapse = "\t"), "' >> ../data/itol/host_legend.txt"))

# Study
study_cols <- brewer.pal(length(unique(jumbo_mcp_meta_edit$study)), "Set1")
names(study_cols) <- unique(jumbo_mcp_meta_edit$study)
study_cols <- study_cols[order(names(study_cols))]
study_cols[["This study"]] <- "rgb(255, 220, 0)"
jumbo_mcp_meta_edit$study_cols <- sapply(jumbo_mcp_meta_edit$study, function(x) study_cols[x])

system("rm ../data/itol/study_legend.txt")
system(paste0("echo 'LEGEND_SHAPES\t", paste0(rep("1", length(study_cols)), collapse = "\t"), "' >> ../data/itol/study_legend.txt"))
system(paste0("echo 'LEGEND_COLORS\t", paste0(study_cols, collapse = "\t"), "' >> ../data/itol/study_legend.txt"))
system(paste0("echo 'LEGEND_LABELS\t", paste0(names(study_cols), collapse = "\t"), "' >> ../data/itol/study_legend.txt"))

# Phage cluster
jumbo_mcp_meta_edit$contig <- sub("_[^_]+$", "", jumbo_mcp_meta_edit$name)
jumbo_mcp_meta_edit$vcontact_cluster[jumbo_mcp_meta_edit$study == "This study"] <- sapply(jumbo_mcp_meta_edit$contig[jumbo_mcp_meta_edit$study == "This study"], function(x) as.character(jumbophage_contig_meta$vcontact_cluster[jumbophage_contig_meta$name %in% x]))
jumbo_mcp_meta_edit$vcontact_cluster[jumbo_mcp_meta_edit$vcontact_cluster %in% ""] <- NA
cluster_cols <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(jumbo_mcp_meta_edit$vcontact_cluster[!is.na(jumbo_mcp_meta_edit$vcontact_cluster)])))
#plot(1:length(cluster_cols), rep(0, length(cluster_cols)), col = cluster_cols, pch = 20, cex = 3)
names(cluster_cols) <- sort(unique(jumbo_mcp_meta_edit$vcontact_cluster[!is.na(jumbo_mcp_meta_edit$vcontact_cluster)]))
jumbo_mcp_meta_edit$cluster_cols <- NA
jumbo_mcp_meta_edit$cluster_cols[!is.na(jumbo_mcp_meta_edit$vcontact_cluster)] <- sapply(jumbo_mcp_meta_edit$vcontact_cluster[!is.na(jumbo_mcp_meta_edit$vcontact_cluster)], function(x) cluster_cols[names(cluster_cols) %in% x])
jumbo_mcp_meta_edit$vcontact_cluster[is.na(jumbo_mcp_meta_edit$vcontact_cluster)] <- jumbo_mcp_meta_edit$contig[is.na(jumbo_mcp_meta_edit$vcontact_cluster)]

jumbo_mcp_meta_edit <- jumbo_mcp_meta_edit %>% select(-"contig")
write.table(jumbo_mcp_meta_edit, "../data/jumbophage_mcp_metadata_edit_cols.txt", row.names = FALSE, quote = FALSE, sep = "\t")
