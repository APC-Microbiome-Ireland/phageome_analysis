library(dplyr)
library(reshape2)

# Read in metadata and data
metadata = read.csv("../data/metadata_v2.csv", stringsAsFactors = FALSE)
contig_ids = read.table("../data/virsorter_positive.ids")
contig_ids = unique(as.character(contig_ids$V1))
contig_data <- readRDS("../data/contig_data.RDS")
jumbophage_contigs <- read.table("../data/jumbophage_contigs.txt", stringsAsFactors = FALSE)

########## Create read data ####
#Read alignment counts
vir_counts = read.table("../data/bowtie2_read_counts.txt", header=TRUE, row.names = 1, sep = "\t", check.names = FALSE)
vir_coverage = read.table("../data/breadth_cov_collated.tbl", header=TRUE, row.names = 1, check.names = FALSE)

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
vir_counts_prop_melt_agg <- vir_counts_prop_melt_agg[!is.na(vir_counts_prop_melt_agg$vcontact_cluster),]
vir_counts_prop_melt_agg <- vir_counts_prop_melt_agg[vir_counts_prop_melt_agg$vcontact_cluster != "",]
saveRDS(vir_counts_prop_melt_agg,  file = "../data/vir_counts_prop_melt_agg.RDS")

# Aggregate counts by CRISPR host
vir_counts_prop_melt_agg2 = vir_counts_prop_melt[, sum(value), by=.(Var2, crispr_host, sample_type, Location)]
saveRDS(vir_counts_prop_melt_agg2, file = "../data/vir_counts_prop_melt_agg2.RDS")

# Re-cast counts matrix by vConTACT2 clusters
vir_counts_prop_agg = dcast.data.table(vir_counts_prop_melt_agg, vcontact_cluster ~ Var2, value.var = "V1", fun.aggregate = sum)
vir_counts_prop_agg = as.data.frame(vir_counts_prop_agg)
rownames(vir_counts_prop_agg) = vir_counts_prop_agg[,1]
vir_counts_prop_agg = vir_counts_prop_agg[,-1]
vir_counts_prop_agg = as.matrix(vir_counts_prop_agg)
vir_counts_prop_agg <- vir_counts_prop_agg[,colSums(vir_counts_prop_agg) != 0]
saveRDS(vir_counts_prop_agg,  file = "../data/vir_counts_prop_agg.RDS")

# Re-cast counts matrix by phages
vir_counts_prop_agg_phage <- dcast.data.table(vir_counts_prop_melt, Var1 ~ Var2, value.var = "value", fun.aggregate = sum)
rownames(vir_counts_prop_agg_phage) <- unlist(vir_counts_prop_agg_phage[,1])
vir_counts_prop_agg_phage <- vir_counts_prop_agg_phage[,-1]
vir_counts_prop_agg_phage <- as.matrix(vir_counts_prop_agg_phage)
saveRDS(vir_counts_prop_agg_phage, file = "../data/vir_counts_prop_agg_phage.RDS")

# Size of all contigs
scaf_len_files <- list.files("../db/SCAFFOLD_LENGTHS", pattern = "*txt", full.names = TRUE)
scaf_len <- data.frame()
for (i in 1:length(scaf_len_files)) {
  tmp <- read.delim(scaf_len_files[i], header = FALSE)
  tmp$ID <- gsub(".*_LENGTHS/", "", gsub("_scaffolds.*", "", scaf_len_files[i]))
  scaf_len <- rbind(scaf_len, tmp)
}
names(scaf_len) <- c("size", "ID")
saveRDS(scaf_len, file = "../data/scaffold_lengths.RDS")
