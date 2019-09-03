#############Libraries########################################
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(tsne)
#library(ggMarginal)
library(parallelDist)
library(cluster)
library(limma)
library(dplyr)
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(circlize)
library(ComplexHeatmap)
library(maSigPro)
library(purrr)
library(ggrepel)
library(stringr)
set.seed(1)

########## Read in metadata ####
metadata = read.csv("data/metadata_v2.csv")
rownames(metadata) = metadata$ID

########## Read in and process raw data######################
# #Contigs
# contig_ids = read.table("data/virsorter_positive.ids")
# contig_ids = unique(as.character(contig_ids$V1))
# 
# #Read alignment counts
# vir_counts = read.table("data/bowtie2_read_counts.txt", header=TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# #vir_counts = fread("data/bowtie2_read_counts.txt", header=TRUE, sep = "\t", check.names = FALSE)
# vir_coverage = read.table("data/breadth_cov_collated.tbl", header=TRUE, row.names = 1, check.names = FALSE)
# #vir_coverage = fread("data/breadth_cov_collated.tbl", header=TRUE, sep = "\t", check.names = FALSE)
# counts_total = vir_counts["Total_reads",]
# vir_counts = vir_counts[-c(nrow(vir_counts), nrow(vir_counts)-1),]
# vir_coverage = vir_coverage[rownames(vir_counts), colnames(vir_counts)]
# vir_counts[vir_coverage < 0.75] = 0 #Apply breadth of coverage filter (<75% coverage are removed)
# rm(vir_coverage)
# 
# vir_counts_unaligned = counts_total - apply(vir_counts, 2, sum)
# vir_counts = rbind(vir_counts, vir_counts_unaligned)
# vir_counts_prop = apply(vir_counts, 2, function(x) x/sum(x))
# rm(vir_counts)
# vir_counts_prop = vir_counts_prop[-nrow(vir_counts_prop),]
# 
# #Remove samples with missing total counts
# vir_counts_prop = vir_counts_prop[,-which(apply(vir_counts_prop, 2, function(x) any(is.na(x))))]
# 
# #Remove samples with missing metadata
# vir_counts_prop = vir_counts_prop[,which(colnames(vir_counts_prop) %in% rownames(metadata))]
# 
# #Remove rows with all zeroes
# vir_counts_prop = vir_counts_prop[-which(apply(vir_counts_prop, 1, function(x) all(x == 0))),]
# 
# ########## Read in and process contig annotations############
# 
# contig_data = data.frame(
#   row.names = contig_ids,
#   name = contig_ids,
#   size = as.numeric(sapply(strsplit(contig_ids, split = "_"), "[[", 5)),
#   coverage = as.numeric(sapply(strsplit(contig_ids, split = "_"), "[[", 7)),
#   circular = FALSE,
#   gc_content = NA,
#   hit_to_vir_refseq = FALSE,
#   hit_to_crasslike = FALSE,
#   demovir = NA,
#   vir_refseq = NA,
#   crasslike = NA,
#   ribo_prot_count = 0,
#   pvogs_count = 0,
#   integrase = FALSE,
#   vcontact_cluster = NA,
#   sample = NA,
#   country = NA
# )
# 
# circular = read.table("data/circular_contigs.ids", header=FALSE)
# vir_refseq = read.table("data/refseq_annot.txt", sep = "\t", header=FALSE)
# crasslike = read.table("data/crasslike_annot.txt", sep = "\t", header=FALSE)
# demovir = read.table("data/DemoVir_assignments.txt", sep = "\t", header=TRUE, row.names = 1)
# gc_content = read.table("data/virsorter_positive_gc_content.txt", header=FALSE, row.names = 1)
# ribo_prot_count = read.table("data/ribosomal_prot_counts.txt", header=FALSE)
# rownames(ribo_prot_count) = as.character(ribo_prot_count$V2)
# pvogs_count = read.table("data/pvogs_counts_per_contig.txt", header=FALSE)
# rownames(pvogs_count) = as.character(pvogs_count$V2)
# vcontact = read.table("data/new_clusters_vcontact_w_taxa.txt", header=TRUE, sep = "\t")
# rownames(vcontact) = as.character(vcontact$contig_ID)
# integrase = read.table("data/integrase_contigs.ids", header=FALSE)
# contig_data[rownames(contig_data) %in% circular$V1,"circular"] = TRUE
# contig_data[rownames(contig_data) %in% integrase$V1,"integrase"] = TRUE
# contig_data$gc_content = as.numeric(as.character(gc_content[rownames(contig_data), "V2"]))
# contig_data[rownames(contig_data) %in% vir_refseq$V1,"hit_to_vir_refseq"] = TRUE
# contig_data[rownames(contig_data) %in% crasslike$V1,"hit_to_crasslike"] = TRUE
# contig_data$demovir = demovir[rownames(contig_data),"Family"]
# contig_data$demovir = as.character(contig_data$demovir)
# contig_data[contig_data$demovir %in%
#               c("Flaviviridae", "Iridoviridae", "Poxviridae", "Retroviridae", "Mimiviridae", "Pithoviridae",
#                 "Ascoviridae", "Nudiviridae", "Baculoviridae", "Marseilleviridae", "Picornaviridae",
#                 "Phycodnaviridae", "Herpesviridae", "Alloherpesviridae"), "demovir"] = "Other"
# contig_data[contig_data$demovir %in% c(NA, "Unassigned"), "demovir"] = "Unassigned"
# contig_data[contig_data$hit_to_crasslike, "demovir"] = "crAss-like"
# contig_data$demovir = as.factor(contig_data$demovir)
# contig_data$demovir = factor(contig_data$demovir, levels =
#                                c("Bicaudaviridae", "Tectiviridae", "Inoviridae", "Microviridae",
#                                  "Siphoviridae", "Myoviridae", "Podoviridae", "crAss-like",
#                                  "Other", "Unassigned")) #Manually reordering levels. Check carefully!!!
# 
# for(i in levels(vir_refseq$V1))
# {
#   contig_data[i, "vir_refseq"] = as.character(vir_refseq[vir_refseq$V1 == i,][1,2]) #use first hit
# }
# contig_data$vir_refseq = as.factor(contig_data$vir_refseq)
# 
# for(i in levels(crasslike$V1))
# {
#   contig_data[i, "crasslike"] = as.character(crasslike[crasslike$V1 == i,][1,2]) #use first hit
# }
# contig_data$crasslike = as.factor(contig_data$crasslike)
# 
# contig_data$ribo_prot_count = ribo_prot_count[rownames(contig_data), "V1"]
# contig_data$pvogs_count = pvogs_count[rownames(contig_data), "V1"]
# contig_data$vcontact_cluster = vcontact[rownames(contig_data), "Cluster_id"]
# 
# country = read.table("data/sample_countries.txt")
# rownames(country) = country$V2
# contig_data$sample = sapply(strsplit(as.character(contig_data$name), split = "_"), "[[", 1)
# 
# #contig_data[contig_data$sample %in% unlist(merged_samples), "sample"] =
# #  as.character(sapply(contig_data[contig_data$sample %in% unlist(merged_samples), "sample"],
# #                      function(x) merged_samples[[grep(x, merged_samples)]][1]))
# 
# contig_data$country = metadata[contig_data$sample,"Location"]
# contig_data[which((is.na(contig_data$country))),"country"] =
#   country[contig_data[which((is.na(contig_data$country))),"sample"],"V1"]
# 
# #Filter contigs based on number of ribosomal protein hits and pVOGs
# contig_data[is.na(contig_data$ribo_prot_count), "ribo_prot_count"] = 0
# contig_data[is.na(contig_data$pvogs_count), "pvogs_count"] = 0
# #contig_data = contig_data[which(contig_data$ribo_prot_count <= 3),]
# #contig_data = contig_data[which(contig_data$pvogs_count/(contig_data$size/1000) >= 0.3),]
# 
# ########## Host information from CRISPR and tRNA matches#################
# #Read in and process CRISPR BLAST data
# crispr_bitscore_cutoff = 45
# trna_bitscore_cutoff = 65
# crispr_evalue_cutoff = 1E-05
# trna_evalue_cutoff = 1E-10
# 
# #crispr_hits1 =
# #  read.table(file = 'refseq89_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)
# crispr_hits2 =
#   read.table(file = 'hmp_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)
# 
# #crispr_hits = rbind(crispr_hits1, crispr_hits2)
# crispr_hits = crispr_hits2
# 
# colnames(crispr_hits) = c("qseqid", "sseqid", "bitscore", "pident", "length", "qlen", "evalue", "sstart", "send", "sname")
# crispr_hits = crispr_hits[crispr_hits$evalue < crispr_evalue_cutoff,]
# crispr_hits = crispr_hits[crispr_hits$bitscore >= crispr_bitscore_cutoff,]#Remove low score hits
# crispr_hits$sname = gsub("UNVERIFIED: ", "", crispr_hits$sname)
# crispr_hits$sname = gsub("TPA: ", "", crispr_hits$sname)
# crispr_hits$sname = gsub("Candidatus ", "", crispr_hits$sname)
# crispr_hits$sname = gsub("^complete chromosome ", "", crispr_hits$sname)
# crispr_hits$genus =
#   sapply(strsplit(as.character(crispr_hits$sname), split = " "), "[[", 2) #Check manually!
# 
# contig_data$crispr_host = NA
# for(i in unique(as.character(crispr_hits$sseqid))) {
#   temp = crispr_hits[crispr_hits$sseqid == i,c("genus", "bitscore")]
#   temp = aggregate(bitscore ~ genus, data = temp, sum)
#   contig_data[i,"crispr_host"] = temp[which.max(temp$bitscore),"genus"]
# }
# 
# ########## Taxonomy from co-clustering with RefSeq######################
# clusters_tax = list()
# for(i in unique(as.character(contig_data$vcontact_cluster))) {
#   tax_string = as.character(vcontact[which(vcontact$Cluster_id == i),"Taxonomy"])
#   if(length(tax_string[!is.na(tax_string)]) > 0) {
#     clusters_tax[[i]] = tax_string[!is.na(tax_string)]
#   }
# }
# sink("clusters_taxonomy.txt")
# print(clusters_tax)
# sink()
# 
# #Dereplicate contigs according to the set used in read alignment step
# contig_data = contig_data[rownames(vir_counts_prop),]
# 
# #Cleanup
# rm(integrase,circular,vir_refseq,crasslike,demovir,gc_content,ribo_prot_count,
#    pvogs_count,vcontact,crispr_hits2)
# #saveRDS(contig_data,  file = "data/contig_data.RDS")
# #write.table(contig_data, file = "data/contig_data.txt")

########## Contig scatterplot####################################
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

########## Virome composition barplots###########################
#Make long format table
vir_counts_prop_melt = melt(vir_counts_prop)
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, contig_data[,c("name", "demovir", "vcontact_cluster", "crispr_host", "integrase")], by = c("Var1"="name"))
vir_counts_prop_melt <- left_join(vir_counts_prop_melt, metadata[,c("ID", "sample_type", "Location")], by = c("Var2"="ID"))

#Aggregate counts by Demovir/vConTACT2
vir_counts_prop_melt <- vir_counts_prop_melt %>% filter(value != 0)
#saveRDS(vir_counts_prop_melt, file = "data/vir_counts_prop_melt.RDS")
vir_counts_prop_melt_agg <- vir_counts_prop_melt %>% group_by(Var2, demovir, vcontact_cluster, sample_type, Location) %>%
  summarise(V1 = sum(value))
#vir_counts_prop_melt_agg = vir_counts_prop_melt[, sum(value), by=.(Var2, demovir, vcontact_cluster, sample_type, Location)]
#saveRDS(vir_counts_prop_melt_agg,  file = "data/vir_counts_prop_melt_agg.RDS")

#Aggregate counts by CRISPR host
vir_counts_prop_melt_agg2 = vir_counts_prop_melt[, sum(value), by=.(Var2, crispr_host, sample_type, Location)]
vir_counts_prop_melt_agg2$crispr_host = as.factor(vir_counts_prop_melt_agg2$crispr_host)
#saveRDS(vir_counts_prop_melt_agg2, file = "data/vir_counts_prop_melt_agg2.RDS")
rm(vir_counts_prop_melt)

barplot_demovir = ggplot(vir_counts_prop_melt_agg, aes(x = Var2, y = V1, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(brewer.pal(12, "Paired")[c(1,2,4,5,7,8,9,10,12)], "lightgray")) +
  facet_wrap(~sample_type, scales = "free_x")

########## Virome beta-diversity###########################
#Re-cast counts matrix by vConTACT2 clusters
vir_counts_prop_agg = dcast.data.table(vir_counts_prop_melt_agg, vcontact_cluster ~ Var2, value.var = "V1", fun.aggregate = sum)
vir_counts_prop_agg = as.data.frame(vir_counts_prop_agg)
vir_counts_prop_agg = vir_counts_prop_agg[-which(is.na(vir_counts_prop_agg[,1])),]
rownames(vir_counts_prop_agg) = vir_counts_prop_agg[,1]
vir_counts_prop_agg = vir_counts_prop_agg[,-1]
vir_counts_prop_agg = as.matrix(vir_counts_prop_agg)

#saveRDS(vir_counts_prop_agg,  file = "data/vir_counts_prop_agg.RDS")

# Compute distance matrix and run t-SNE
vir_counts_prop_agg <- readRDS("data/vir_counts_prop_agg.RDS")
vir_counts_prop_agg <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% metadata$ID]
vir_counts_prop_agg <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% metadata$ID[metadata$Visit_Number == 1]]

# T-sne
set.seed(1)
cluster_counts_dist = vegdist(t(vir_counts_prop_agg), method = "bray")
clusters_samples_tsne = tsne(cluster_counts_dist)
# # PCA
# clusters_samples_pca <- prcomp(cluster_counts_dist, scale. = TRUE, center = TRUE)

#Annotate t-SNE coordinates
clusters_samples_tsne = as.data.frame(clusters_samples_tsne)
rownames(clusters_samples_tsne) = colnames(vir_counts_prop_agg)
clusters_samples_tsne$ID <- rownames(clusters_samples_tsne)
colnames(clusters_samples_tsne) = c("tsne1", "tsne2", "ID")
clusters_samples_tsne <- left_join(clusters_samples_tsne, metadata, by = "ID")

# #De novo clustering
# samples_clust = pam(cluster_counts_dist, 7)

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 2:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(cluster_counts_dist, k = k), cluster_counts_dist)
  avg_sil[k-1] <- mean(tmp[,3])
}
# New clusters
samples_clust <- pam(cluster_counts_dist, which.max(avg_sil)+1)

# Read contig data
contig_data <- readRDS("data/contig_data.RDS")

clusters_samples_tsne$cluster = as.factor(samples_clust$cluster[clusters_samples_tsne$ID])
clusters_samples_tsne$n_contigs = sapply(rownames(clusters_samples_tsne), function(x) length(which(contig_data$sample == x)))
clusters_samples_tsne$total_reads = as.numeric(counts_total[clusters_samples_tsne$ID])

#saveRDS(clusters_samples_tsne,  file = "data/clusters_samples_tsne.RDS")

#Make t-SNE plots
clusters_samples_tsne <- readRDS("data/clusters_samples_tsne.RDS")
tsne_plot1 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = cluster, shape = Location)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual("Country", values = c(21, 24, 22)) +
  scale_fill_manual(values = brewer.pal(8, "Set2"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("Dim 1") + ylab("Dim 2")

tsne_plot2 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = n_contigs, size = log10(total_reads), shape = Location)) +
  theme_classic() + geom_point(alpha = 0.7) +
  scale_shape_manual("Country", values = c(21, 24, 22)) +
  xlab("Dim 1") + ylab("Dim 2")

tsne_plot3 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = sample_type, shape = Location)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual("Country", values = c(21, 24, 22)) +
  scale_fill_brewer("body site", palette = "Set3",
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("Dim 1") + ylab("Dim 2")

tsne_plot4 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = Age, shape = Gender)) +
  theme_classic() + geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual("Sex", values = c(21, 22)) +
  #scale_fill_gradient2()
  scale_fill_distiller(palette = "Spectral") +
  xlab("Dim 1") + ylab("Dim 2")

########## GLM of variables###########
glm_samples <- clusters_samples_tsne[!is.na(clusters_samples_tsne$Age),]
glm(cluster ~ Location + sample_type + Gender + Age + Health, family = "binomial", data = glm_samples)
# Difficult as variables colinear

########## Differentially abundant clusters#####################
clusters_samples_tsne <- readRDS("data/clusters_samples_tsne.RDS")
vir_counts_prop_agg <- readRDS("data/vir_counts_prop_agg.RDS")
vir_counts_prop_agg <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID]
contig_data <- readRDS("data/contig_data.RDS")
clusters_samples_kw = apply(vir_counts_prop_agg, 1, function(x) kruskal.test(x, clusters_samples_tsne$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
vir_counts_prop_agg_diff = vir_counts_prop_agg[names(which(clusters_samples_kw < 0.001)),]

viral_clusters_df = data.frame(row.names = rownames(vir_counts_prop_agg),
                               demovir = sapply(rownames(vir_counts_prop_agg), function(x) names(which.max(table(contig_data[which(contig_data$vcontact_cluster == x),"demovir"]))))
)

########## Differentially abundant megaphages ####
clusters_samples_tsne <- readRDS("data/clusters_samples_tsne.RDS")
vir_counts_prop_agg <- readRDS("data/vir_counts_prop_agg.RDS")
megaphage_contigs <- read.table("data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
megaphage_contigs <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
megaphage_contigs$Health <- as.character(megaphage_contigs$Health)
megaphage_contigs$Health[megaphage_contigs$Health == "unhealthy"] <- "rheumatoid arthritis"
megaphage_contigs <- megaphage_contigs[megaphage_contigs$Visit_Number == 1,]
megaphage_contigs <- megaphage_contigs[!is.na(megaphage_contigs$circular),]
megaphage_contigs <- megaphage_contigs[megaphage_contigs$circular,]
vir_counts_prop_agg <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID]
vir_counts_prop_agg <- vir_counts_prop_agg[row.names(vir_counts_prop_agg) %in% megaphage_contigs$vcontact_cluster,]
clusters_samples_kw = apply(vir_counts_prop_agg, 1, function(x) kruskal.test(x, clusters_samples_tsne$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
vir_counts_prop_agg_diff = vir_counts_prop_agg[names(which(clusters_samples_kw < 0.001)),]

megaphage_clusters_df = data.frame(row.names = rownames(vir_counts_prop_agg),
                                   demovir = sapply(rownames(vir_counts_prop_agg), function(x) names(which.max(table(megaphage_contigs[which(megaphage_contigs$vcontact_cluster == x),"demovir"])))))

# Taxa list
vcontact = read.table("data/new_clusters_vcontact_w_taxa.txt", header=TRUE, sep = "\t")
rownames(vcontact) = as.character(vcontact$contig_ID)
clusters_tax_megaphage = list()
for(i in unique(as.character(megaphage_contigs$vcontact_cluster))) {
  tax_string = as.character(vcontact[which(vcontact$Cluster_id == i),"Taxonomy"])
  if(length(tax_string[!is.na(tax_string)]) > 0) {
    clusters_tax_megaphage[[i]] = tax_string[!is.na(tax_string)]
  }
}

# Clusters heatmap
vir_counts_prop_agg_diff_log <- log10(vir_counts_prop_agg_diff)
#vir_counts_prop_agg_diff_log <- vir_counts_prop_agg_diff_log[,colSums(vir_counts_prop_agg_diff_log) != 0]
vir_counts_prop_agg_diff_log[is.infinite(vir_counts_prop_agg_diff_log)] = -8

# Heatmap
tiff(file = "figures/heatplot_clusters_megaphage.tiff", width=2000, height=1500, res=150)
heatmap.2(vir_counts_prop_agg_diff_log,
          margins = c(5,20),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(min(vir_counts_prop_agg_diff_log), 0, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA,
          key.xlab = "log10(count)",
          RowSideColors = c(brewer.pal(12, "Paired"), "lightgrey")[c(10,4,5,8,12,9,7,2,13)][as.numeric(megaphage_clusters_df[rownames(vir_counts_prop_agg_diff),"demovir"])],
          ColSideColors = brewer.pal(8, "Set2")[as.numeric(metadata[colnames(vir_counts_prop_agg_diff), "sample_type"])],
          labRow = gsub("NULL", "", as.character(sapply(clusters_tax_megaphage[rownames(vir_counts_prop_agg_diff)], "[[", 1))),
          labCol = NA,
          cexCol = 1,
          cexRow = 0.5,
          xlab = "Samples",
          ylab = "Contig clusters",
          main = NA
)
legend("topright", legend = levels(metadata[colnames(vir_counts_prop_agg_diff), "sample_type"]),
       col = brewer.pal(8, "Set2"), bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.8, title = "Body Site")
legend("bottomright", legend = levels(megaphage_clusters_df[rownames(vir_counts_prop_agg_diff),"demovir"]),
       col = c(brewer.pal(12, "Paired"), "lightgrey")[c(10,4,5,8,12,9,7,2,13)], bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.8, title = "Predicted Taxonomy")
dev.off()

########## Host range############################
#Re-cast counts matrix by CRISPR hosts
vir_counts_prop_agg2 = dcast.data.table(vir_counts_prop_melt_agg2, crispr_host ~ Var2, value.var = "V1", fun.aggregate = sum)
vir_counts_prop_agg2 = as.data.frame(vir_counts_prop_agg2)
vir_counts_prop_agg2 = vir_counts_prop_agg2[-which(is.na(vir_counts_prop_agg2[,1])),]
rownames(vir_counts_prop_agg2) = vir_counts_prop_agg2[,1]
vir_counts_prop_agg2 = vir_counts_prop_agg2[,-1]
vir_counts_prop_agg2 = as.matrix(vir_counts_prop_agg2)

#saveRDS(vir_counts_prop_agg2,  file = "data/vir_counts_prop_agg2.RDS")

########## Saving output#######################################
pdf(file="figures/contigs_scatter.pdf", paper="A4r", width=11, height=8.5)
contig_scatter_demovir
dev.off()

pdf(file="figures/samples_barplot.pdf", paper="A4r", width=11, height=8.5)
barplot_demovir
dev.off()

pdf(file="figures/samples_tsne.pdf", paper="A4r", width=11, height=8.5)
grid.arrange(tsne_plot1, tsne_plot2, tsne_plot3, tsne_plot4, nrow = 2)
dev.off()

tiff("figures/Figure1a.tiff", width = 1500, height = 500, res = 150)
grid.arrange(tsne_plot1, tsne_plot3, nrow = 1)
dev.off()

tiff("figures/Supplementary_Figure1.tiff", width = 1500, height = 500, res = 100)
grid.arrange(tsne_plot2, tsne_plot4, nrow = 1)
dev.off()

#Clusters heatmap
vir_counts_prop_agg_diff_log = log10(vir_counts_prop_agg_diff+1)
vir_counts_prop_agg_diff_log[is.infinite(vir_counts_prop_agg_diff_log)] = -8

pdf(file = "figures/heatplot_clusters.pdf", width=8.5, height=11, pointsize = 12)
heatmap.2(vir_counts_prop_agg_diff_log,
          margins = c(10,10),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(-8, 0, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA,
          RowSideColors = c(brewer.pal(12, "Paired"), "lightgrey")[c(10,4,5,8,12,9,7,2,13)][as.numeric(viral_clusters_df[rownames(vir_counts_prop_agg_diff),"demovir"])],
          ColSideColors = brewer.pal(9, "Set1")[as.numeric(metadata[colnames(vir_counts_prop_agg_diff), "sample_type"])],
          labRow = gsub("NULL", "", as.character(sapply(clusters_tax[rownames(vir_counts_prop_agg_diff)], "[[", 1))),
          labCol = NA,
          cexCol = 0.5,
          cexRow = 0.5,
          xlab = "Samples",
          ylab = "Contig clusters",
          main = NA
)
legend(x = 0.01, y = 0.8, legend = levels(metadata[colnames(vir_counts_prop_agg_diff), "sample_type"]),
       col = brewer.pal(9, "Set1"), bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Sample clusters:")
legend(x = 0.01, y = 0.2, legend = levels(viral_clusters_df[rownames(vir_counts_prop_agg_diff),"demovir"]),
       col = c(brewer.pal(12, "Paired"), "lightgrey")[c(10,4,5,8,12,9,7,2,13)], bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.4, title = "Taxonomy by Demovir")
dev.off()

#Hosts heatmap
vir_counts_prop_agg2_log = log10(vir_counts_prop_agg2)
vir_counts_prop_agg2_log[is.infinite(vir_counts_prop_agg2_log)] = -8

pdf(file = "figures/heatplot_hosts.pdf", width=8.5, height=11, pointsize = 12)
heatmap.2(vir_counts_prop_agg2_log,
          margins = c(10,10),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(-8, 0, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA,
          ColSideColors = brewer.pal(9, "Set1")[as.numeric(metadata[colnames(vir_counts_prop_agg2_log), "sample_type"])],
          labRow = rownames(vir_counts_prop_agg2_log),
          labCol = NA,
          cexCol = 0.5,
          cexRow = 1,
          xlab = "Samples",
          ylab = "Predicted hosts",
          main = NA
)
legend(x = 0.01, y = 0.8, legend = levels(metadata[colnames(vir_counts_prop_agg2_log), "sample_type"]),
       col = brewer.pal(9, "Set1"), bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Sample clusters:")
dev.off()

save.image()

########## Venn diagram showing overview of total phages for each group ##########
unique_clusters <- unique(clusters_samples_tsne$cluster)

for (i in 1:length(unique_clusters)) {
  vir_tmp <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == unique_clusters[i]]]
  vir_group <- rownames(vir_tmp)[rowSums(vir_tmp) != 0]
}

vir_group_list <- sapply(unique_clusters, function(x) {
  vir_tmp <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == x]]
  vir_group <- rownames(vir_tmp)[rowSums(vir_tmp) != 0]
  return(vir_group)
  })
names(vir_group_list) <- c("Group 1\n(dorsum \nof tongue \n& saliva)", "Group 2\n(buccal mucosa)", "Group 3\n(dental)", "Group 4\n(stool)")
max_group_length <- max(sapply(vir_group_list, function(x) length(x)))

vir_group_list <- lapply(vir_group_list, function(x) c(x, rep(NA, max_group_length - length(x))))

vir_group_all <- do.call(cbind, vir_group_list)

tiff("figures/Figure1b.tiff", width = 600, height = 600)
suma2Venn(vir_group_all, zcolor = brewer.pal(length(unique_clusters), "Set2"), cexil = 1.5, cexsn = 1.4)
dev.off()

########## Longitudinal analysis and core phageome #########
# Create mock IDs so different time points
metadata_us <- metadata[metadata$Location == "US",]
metadata_longus <- metadata_us %>% group_by(Sample.name, sample_type) %>%
  filter(any(Visit_Number == 3))

vir_counts_prop_melt <- readRDS("data/vir_counts_prop_melt.RDS")
vir_counts_longus <- left_join(metadata_longus[,c("ID", "Sample.name", "Visit_Number")], vir_counts_prop_melt, by = c("ID"="Var2"))
vir_counts_longus <- vir_counts_longus[!is.na(vir_counts_longus$Var1),]

# countConsecutive <- function(dt) {
#   min_timepoint <- min(dt$Visit_Number)
#   max_timepoint <- max(dt$Visit_Number)
#   count <- 0
#   final_counts <- c()
#   for (i in min_timepoint:max_timepoint) {
#     if (i %in% dt$Visit_Number) {
#       count <- count + 1
#     } else {
#       final_counts <- c(final_counts, count)
#       count <- 0
#     }
#   }
#   final_counts <- c(final_counts, count)
#   final_counts <- final_counts[final_counts != 0]
#   final_counts_df <- as.data.frame(mean(final_counts))
#   names(final_counts_df) <- "avg_consecutive_tp"
#   return(final_counts_df)
# }

# Count average consecutive time points of phages
contig_no_tp <- vir_counts_longus %>% group_by(sample_type, Sample.name, Var1) %>%
  summarise(no_tp = n_distinct(Visit_Number)) 

# Count all distinct contigs
contig_all_summary <- vir_counts_longus %>% group_by(sample_type, Sample.name) %>%
  summarise(total_contig = n_distinct(Var1))

# Count contigs in at least x number of time points
contig_cum_tp <- map_df(.x = unique(contig_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(contig_no_tp$sample_type[contig_no_tp$no_tp == .x], .x),
                    Sample.name = rep(contig_no_tp$Sample.name[contig_no_tp$no_tp == .x], .x),
                    Var1 = rep(contig_no_tp$Var1[contig_no_tp$no_tp == .x]),
                    no_tp = rep(contig_no_tp$no_tp[contig_no_tp$no_tp == .x], .x))
  tmp$cum_tp <- c(1:.x)
  return(tmp)
})

contig_cum_tp_summary <- contig_cum_tp %>% group_by(sample_type, Sample.name, cum_tp) %>%
  summarise(total_contig_least = n_distinct(Var1))

contig_one <- contig_cum_tp_summary %>%
  filter(cum_tp == 1) %>%
  rename(total_contig = total_contig_least) %>%
  select(-cum_tp)

contig_cum_tp_summary <- left_join(contig_cum_tp_summary, contig_one, by = c("sample_type", "Sample.name")) %>%
  mutate(contig_frac = total_contig_least/total_contig)

# Plot
tiff("figures/Figure2a.tiff", height = 400, width = 1100, res = 100)
ggplot(contig_cum_tp_summary, aes(as.factor(cum_tp), contig_frac, fill = sample_type)) +
  geom_boxplot() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of contigs in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = brewer.pal(length(unique(metadata_us$sample_type)), "Spectral"))
dev.off()

# Proportion of reads mapped to contigs
contig_cum_reads<- vir_counts_longus %>% 
  group_by(sample_type, Sample.name, Var1) %>%
  arrange(desc(Visit_Number)) %>%
  mutate(cumsum_reads = cumsum(value)) %>%
  select(sample_type, Sample.name, Var1, cumsum_reads, Visit_Number)

contig_cum_reads_tp1 <- contig_cum_reads %>%
  filter(Visit_Number == min(Visit_Number)) %>%
  rename(tp1 = cumsum_reads) %>%
  select(sample_type, Sample.name, Var1, tp1)

contig_cum_reads_summary <- left_join(contig_cum_reads, contig_cum_reads_tp1, by = c("sample_type", "Var1", "Sample.name")) %>%
  mutate(frac_reads = cumsum_reads/tp1) %>%
  group_by(sample_type, Sample.name, Var1) %>%
  mutate(cum_tp = order(Visit_Number))

tiff("figures/Figure2b.tiff", height = 400, width = 1100, res = 90)
ggplot(contig_cum_reads_summary, aes(x = as.factor(cum_tp), y = frac_reads, fill = sample_type)) +
  geom_boxplot() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of reads mapped in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = brewer.pal(length(unique(metadata_us$sample_type)), "Spectral"))
dev.off()

# Select persistent contigs (3 at least, 2 doesnt converge)
contig_persist <- left_join(vir_counts_longus, contig_no_tp, by = c("sample_type", "Sample.name", "Var1")) %>%
  filter(no_tp >= 3)
contig_persist <- dcast(contig_persist, ID ~ Var1, sum, value.var = "value") 
rownames(contig_persist) <- contig_persist$ID
contig_persist <- contig_persist[,names(contig_persist) != "ID"]

# Run NMDS
set.seed(1)
nmds_persist <- metaMDS(contig_persist, distance = "bray", k = 2, trymax = 20)
df_nmds_persist <- as.data.frame(nmds_persist$points)
df_nmds_persist$ID <- row.names(df_nmds_persist)
df_nmds_persist$Sample.name <- as.character(sapply(df_nmds_persist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_nmds_persist$Location <- sapply(df_nmds_persist$ID, function(x) metadata$Location[metadata$ID == x])
df_nmds_persist$sample_type <- sapply(df_nmds_persist$ID, function(x) metadata$sample_type[metadata$ID == x])

# Plot NMDS
group_colours <- brewer.pal(length(unique(df_nmds_persist$sample_type)), "Spectral")
tiff("figures/Figure2c.tiff", width = 2000, height = 1000, res = 200)
ggplot(df_nmds_persist, aes(MDS1, MDS2, fill = Sample.name, colour = sample_type)) +
  geom_point() +
  geom_density2d(alpha=0.3) +
  theme_classic() +
  guides(fill = FALSE) +
  scale_colour_manual("body site", values = group_colours) +
  scale_x_continuous(limits = c(-2,2), breaks = seq(-2,2,0.5)) + 
  scale_y_continuous(limits = c(-1.5,1.2), breaks = seq(-2,2,0.5))
dev.off()

# Select nonpersistent contigs
contig_nonpersist <- left_join(vir_counts_longus, contig_no_tp, by = c("sample_type", "Sample.name", "Var1")) %>%
  filter(no_tp < 3)
contig_nonpersist <- dcast(contig_nonpersist, ID ~ Var1, sum, value.var = "value") 
rownames(contig_nonpersist) <- contig_nonpersist$ID
contig_nonpersist <- contig_nonpersist[,names(contig_nonpersist) != "ID"]

# Run NMDS
set.seed(1)
nmds_nonpersist <- metaMDS(contig_nonpersist, distance = "bray", k = 2, trymax = 20)
df_nmds_nonpersist <- as.data.frame(nmds_nonpersist$points)
df_nmds_nonpersist$ID <- row.names(df_nmds_nonpersist)
df_nmds_nonpersist$Sample.name <- as.character(sapply(df_nmds_nonpersist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_nmds_nonpersist$Location <- sapply(df_nmds_nonpersist$ID, function(x) metadata$Location[metadata$ID == x])
df_nmds_nonpersist$sample_type <- sapply(df_nmds_nonpersist$ID, function(x) metadata$sample_type[metadata$ID == x])

# Plot NMDS
tiff("figures/Figure2d.tiff", width = 2000, height = 700, res = 200)
ggplot(df_nmds_nonpersist, aes(MDS1, MDS2, fill = Sample.name, colour = sample_type)) +
  geom_point() +
  geom_density2d(alpha=0.2) +
  theme_classic() +
  guides(fill = FALSE) +
  scale_colour_manual("body site", values = group_colours) +
  scale_x_continuous(limits = c(-1.5,1), breaks = seq(-2,2,0.5)) + 
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-2,2,0.5))
dev.off()

########## Procrustes analysis #########
metaphlan <- read.csv("data/all_metaphlan.csv", stringsAsFactors = FALSE)
row.names(metaphlan) <- metaphlan$X
metaphlan <- metaphlan[,names(metaphlan) != "X"]
metaphlan <- metaphlan[,names(metaphlan) %in% colnames(vir_counts_prop_agg)]
metaphlan_filter <- metaphlan[grepl("s__", row.names(metaphlan)) & !grepl("t__", row.names(metaphlan)),]
metaphlan_dist = vegdist(t(metaphlan_filter), method = "bray")
set.seed(1)
metaphlan_tsne <- tsne(metaphlan_dist)

protest_res <- protest(clusters_samples_tsne[,c(1:2)], metaphlan_tsne[,c(1:2)], scale = TRUE)

tiff("figures/Figure1d.tiff")
plot(protest_res)
points(protest_res, display = "target", col = "red")
dev.off()

########## Correlation analysis ###########

runSpearmanCorrelation <- function(phage_data, taxa_data) {
  
  phage_data <- phage_data[rowSums(phage_data != 0) > ncol(phage_data)/2,]
  taxa_data <- taxa_data[rowSums(taxa_data != 0) > ncol(taxa_data)/2,]
  
  cor_data <- c()
  cor_test <- c()
  x <- c()
  y <- c()
  # Spearman Correlation p values
  for(i in 1:nrow(taxa_data)){
    for(j in 1:nrow(phage_data)){
      cor <- cor.test(as.numeric(taxa_data[i,]), as.numeric(phage_data[j,]), method = "spearman")
      cor_data <- c(cor_data, cor$estimate)
      cor_test <- c(cor_test, cor$p.value)
      x <- c(x, row.names(taxa_data)[i])
      y <- c(y, row.names(phage_data)[j])
    }
  }
  
  # Remove repeat tests and direct comparisons
  cor_df <- data.frame(x = x, y = y, rho = cor_data, p_val = cor_test)
  cor_df <- cor_df[paste0(x, y) != paste0(y, x),]
  
  # Benjamini-Hochberg correction
  cor_df$FDR <- p.adjust(cor_df$p_val, method = "BH")
  cor_df <- cor_df[cor_df$FDR < 0.05,]
  
  # Re-label
  cor_df$phylum <- sapply(as.character(cor_df$x), function(x) strsplit(x, "p__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$class <- sapply(as.character(cor_df$x), function(x) strsplit(x, "c__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$order <- sapply(as.character(cor_df$x), function(x) strsplit(x, "o__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$family <- sapply(as.character(cor_df$x), function(x) strsplit(x, "f__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$genus <- sapply(as.character(cor_df$x), function(x) strsplit(x, "g__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$species <- sapply(as.character(cor_df$x), function(x) strsplit(x, "s__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$strain <- sapply(as.character(cor_df$x), function(x) strsplit(x, "t__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$x <- gsub("Streptococcus_mitis_oralis_pneumoniae", "Streptococcus_mitis/oralis/pneumoniae", cor_df$x)
  
  return(cor_df)
}

correlationHeatmap <- function(cor_df, bottom_margin = 0, left_margin = 0, phyla, demovir){
  sparse_matrix <- dcast(data = cor_df, formula = y + demovir ~ species + phylum, fun.aggregate = sum, value.var = "rho")
  row.names(sparse_matrix) <- paste0(sparse_matrix$y, "_", sparse_matrix$demovir)
  sparse_matrix <- sparse_matrix[,-c(1,2)]
  
  row_cols <- c(brewer.pal(length(phyla), "Spectral"))
  names(row_cols) <- phyla
  
  phyla_samples <- rep(NA, ncol(sparse_matrix))
  phyla_cols <- rep(NA, ncol(sparse_matrix))
  for(i in 1:length(phyla)){
    phyla_samples[grep(phyla[i], names(sparse_matrix))] <- phyla[i]
    phyla_cols[grep(phyla[i], names(sparse_matrix))] <- row_cols[names(row_cols) == phyla[i]]
    names(sparse_matrix) <- gsub(paste0("_", phyla[i]), "", names(sparse_matrix))
  }
  names(phyla_cols) <- phyla_samples
  
  ha <- rowAnnotation(type = phyla_samples, col = list(type = phyla_cols),
                         annotation_legend_param = list(type = list(title = "Bacterial Phylum",
                                                                    title_gp = gpar(fontsize = 20),
                                                                    labels_gp = gpar(fontsize = 16),
                                                                    grid_height = unit(8, "mm"))))
  
  col_cols <- c(brewer.pal(length(demovir), "Paired"))
  names(col_cols) <- demovir
  
  demovir_samples <- rep(NA, nrow(sparse_matrix))
  demovir_cols <- rep(NA, nrow(sparse_matrix))
  for(i in 1:length(demovir)){
    demovir_samples[grep(demovir[i], row.names(sparse_matrix))] <- demovir[i]
    demovir_cols[grep(demovir[i], row.names(sparse_matrix))] <- col_cols[names(col_cols) == demovir[i]]
    row.names(sparse_matrix) <- gsub(paste0("_", demovir[i]), "", row.names(sparse_matrix))
  }
  names(demovir_cols) <- demovir_samples
  
  ta <- HeatmapAnnotation(type = demovir_samples, col = list(type = demovir_cols),
                          annotation_legend_param = list(type = list(title = "Viral Family",
                                                                 title_gp = gpar(fontsize = 20),
                                                                 labels_gp = gpar(fontsize = 16),
                                                                 grid_height = unit(8, "mm"))))
  
  set.seed(42)
  ht <- Heatmap(t(sparse_matrix), top_annotation = ta, left_annotation = ha, name = "rho", show_row_names = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, 
                col = colorRamp2(c(-1,0,1),  c("blue", "white", "red")), column_names_max_height = unit(bottom_margin, "mm"),
                row_title_rot = 0, row_title_gp = gpar(fontsize = 5), show_column_names = TRUE, row_names_max_width = unit(left_margin, "mm"),
                heatmap_legend_param = list(color_bar = "continuous", title_gp = gpar(fontsize = 20),
                                            labels_gp = gpar(fontsize=16), legend_height = unit(8, "cm")))
  return(ht)
}

metaphlan_species <- metaphlan[grepl("s__", row.names(metaphlan)),]
viral_clusters_df$y <- row.names(viral_clusters_df)

cor_group1 <- runSpearmanCorrelation(vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 1]],
                                     metaphlan_species[,colnames(metaphlan_species) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 1]])
cor_group2 <- runSpearmanCorrelation(vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 2]],
                                     metaphlan_species[,colnames(metaphlan_species) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 2]])
cor_group3 <- runSpearmanCorrelation(vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 3]],
                                     metaphlan_species[,colnames(metaphlan_species) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 3]])
cor_group4 <- runSpearmanCorrelation(vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 4]],
                                     metaphlan_species[,colnames(metaphlan_species) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == 4]])
cor_all <- runSpearmanCorrelation(vir_counts_prop_agg, metaphlan_species)

save(cor_group1, cor_group2, cor_group3, cor_group4, file = "data/cor_groups.RData")
load("cor_groups.RData")

phyla <- unique(c(cor_group1$phylum, cor_group2$phylum, cor_group3$phylum, cor_group4$phylum))
cor_group1 <- left_join(cor_group1, viral_clusters_df, by = "y")
cor_group2 <- left_join(cor_group2, viral_clusters_df, by = "y")
cor_group3 <- left_join(cor_group3, viral_clusters_df, by = "y")
cor_group4 <- left_join(cor_group4, viral_clusters_df, by = "y")
demovir <- as.character(unique(viral_clusters_df$demovir))

tiff("figures/Supplementary_Figure3a.tiff", width = 1500, height = 1100)
draw(correlationHeatmap(cor_group1[abs(cor_group1$rho) > 0.6,], 100, 100, phyla = phyla, demovir = demovir), column_title = "Group 1")
dev.off()

tiff("figures/Supplementary_Figure3b.tiff", width = 600, height = 300)
draw(correlationHeatmap(cor_group2[abs(cor_group2$rho) > 0.6,], 100, 100, phyla = phyla, demovir = demovir), column_title = "Group 2")
dev.off()

tiff("figures/Supplementary_Figure3c.tiff", width = 1600, height = 1200)
draw(correlationHeatmap(cor_group3[abs(cor_group3$rho) > 0.6,], 100, 100, phyla = phyla, demovir = demovir), column_title = "Group 3")
dev.off()

tiff("figures/Supplementary_Figure3d.tiff", width = 1100, height = 500)
draw(correlationHeatmap(cor_group4[abs(cor_group4$rho) > 0.6,], 100, 100, phyla = phyla, demovir = demovir), column_title = "Group 4")
dev.off()

########## Alpha-diversity and phage richness ####
## Functions ##
createPairedData <- function(df_map, pair_list){
  
  # Create groups e.g. stool vs dental
  df_map_pairs <- data.frame()
  for(i in 1:length(pair_list)){
    samples_one <- unique(df_map$Sample.name[df_map$sample_type == pair_list[[i]][1]])
    samples_two <- unique(df_map$Sample.name[df_map$sample_type == pair_list[[i]][2]])
    df_map_pair <- df_map[(df_map$Sample.name %in% Reduce(intersect,list(samples_one, samples_two))) & (df_map$sample_type %in% pair_list[[i]]),]
    
    df_map_pair$group <- paste(pair_list[[i]][1], "vs.", pair_list[[i]][2])
    if(nrow(df_map_pairs) == 0) {
      df_map_pairs <- df_map_pair
    } else {
      df_map_pairs <- rbind(df_map_pairs, df_map_pair)
    }
  }
  
  # Order by Location and change characters in group
  df_map_pairs <- df_map_pairs[order(df_map_pairs$Location),]
  df_map_pairs$group <- as.character(df_map_pairs$group)
  
  return(df_map_pairs)
}

runTtest <- function(df_paired){
  
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  
  # T-test
  p_values <- c()
  ttest_groups <- df_paired_richness[!duplicated(paste0(df_paired_richness$Location, df_paired_richness$group)),]
  for(i in 1:nrow(ttest_groups)){
    y <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & df_paired_richness$group == ttest_groups$group[i],]
    y <- y[order(y$Sample.name),]
    type1 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][1]
    type2 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][2]
    y1 <- y$richness[y$sample_type == type1]
    y2 <- y$richness[y$sample_type == type2]
    p_values <- c(p_values, wilcox.test(y1, y2, paired = TRUE, alternative = "two.sided")$p.value)
  }
  ttest_groups$pvalue <- p_values
  asterisk <- rep(NA, length(p_values))
  asterisk[p_values < 0.05] <- "*"
  asterisk[p_values < 0.01] <- "**"
  asterisk[p_values < 0.001] <- "***"
  ttest_groups$asterisk <- asterisk
  
  # Order ttest_groups
  ttest_groups <- ttest_groups[order(ttest_groups$Location, ttest_groups$group),]
  ttest_groups <- ttest_groups[,c("Location", "group", "pvalue", "asterisk")]
  return(ttest_groups)
}

plotMultipleRichnessGraphs <- function(ttest_groups, df_paired){
  
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  g <- list()
  for(i in 1:nrow(ttest_groups)){
    df_paired_richness_group <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & df_paired_richness$group == ttest_groups$group[i],]
    g[[i]] <- ggplot(df_paired_richness_group, aes(sample_type, richness)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.8) +
      theme_classic() +
      ylab("MAP Richness") +
      xlab("") +
      ylim(c(0, max(df_paired_richness_group$richness) + 20)) +
      ggtitle(paste(ttest_groups[i,]$Location, sep = "\n")) +
      geom_text(data = ttest_groups[i,], aes(label=asterisk),
                x = 1.5, y = max(df_paired_richness_group$richness)+10, size = 7,
                inherit.aes = FALSE)
  }
  return(g)
}

# Read data
vir_counts_prop_melt <- readRDS("data/vir_counts_prop_melt.RDS")
vir_counts_prop_melt <- left_join(metadata, vir_counts_prop_melt, by = c("ID"="Var2","sample_type","Location"))
vir_counts_prop_melt <- vir_counts_prop_melt[vir_counts_prop_melt$Visit_Number == 1,]

# Get phage cluster richness
richness_summary <- vir_counts_prop_melt %>% 
  group_by(ID, sample_type, Location, Sample.name) %>% 
  summarise(richness = n_distinct(vcontact_cluster))

pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
richness_paired <- createPairedData(richness_summary, pair_list)
richness_ttest <- runTtest(richness_paired)
richness_graphs <- plotMultipleRichnessGraphs(richness_ttest, richness_paired)

# Plot graph
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
tiff("figures/Figure3a.tiff", width = 2100, height = 1500, res = 220)
grid.arrange(grobs = richness_graphs, layout_matrix = lay)
dev.off()

########## Jumbo/mega circular phages ####
megaphage_contigs <- read.delim("data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
megaphage_contigs <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
megaphage_contigs$Health <- as.character(megaphage_contigs$Health)
megaphage_contigs$Health[megaphage_contigs$Health == "unhealthy"] <- "rheumatoid arthritis"
megaphage_contigs <- megaphage_contigs[megaphage_contigs$Visit_Number == 1,]
megaphage_contigs <- megaphage_contigs[!is.na(megaphage_contigs$circular),]
megaphage_contigs <- megaphage_contigs[megaphage_contigs$circular,]
megaphage_contigs$Location_Health <- paste(megaphage_contigs$Location, "-", megaphage_contigs$Health)
megaphage_contigs$numeric_samplename <- as.numeric(factor(megaphage_contigs$Sample.name))
megaphage_contigs$numeric_samplename[!megaphage_contigs$numeric_samplename %in% unique(megaphage_contigs$numeric_samplename[duplicated(megaphage_contigs$numeric_samplename)])] <- NA
megaphage_contigs$numeric_samplename <- as.numeric(factor(megaphage_contigs$numeric_samplename))
megaphage_contigs$numeric_samplename[is.na(megaphage_contigs$numeric_samplename)] <- c((max(megaphage_contigs$numeric_samplename, na.rm = TRUE)+1):(max(megaphage_contigs$numeric_samplename, na.rm = TRUE)+sum(is.na(megaphage_contigs$numeric_samplename))))

label_colours <- c("blue", "orange", "turquoise", "darkgreen")
pd <- position_dodge(0.8)
h <- 200
tiff("figures/megaphage_individual.tiff", width = 2500, height = 1000, res = 200)
ggplot(megaphage_contigs, aes(sample_type, size/1000, colour = Location_Health, group = numeric_samplename)) +
  geom_line(position=pd, colour="grey90", linetype = "dashed") +
  geom_point(position=pd, alpha = 0.6) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  geom_text(aes(1, h, label = "Megaphage - 200 kb", vjust = 1), colour = "red", size = 3) + 
  theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Country - Health", values = label_colours) 
dev.off()

tiff("figures/megaphage_cluster.tiff", width = 2500, height = 1000, res = 200)
ggplot(megaphage_contigs, aes(sample_type, size/1000, colour = Location_Health, group = vcontact_cluster)) +
  geom_line(position=pd, colour="grey90", linetype = "dashed") +
  geom_point(position = pd, alpha = 0.6) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  geom_text(aes(1, h, label = "Megaphage - 200 kb", vjust = 1), colour = "red", size = 3) + 
  theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Country - Health", values = label_colours) 
dev.off()  

########## Phages and ARGs ####
# Read data
contig_data <- readRDS("data/contig_data.RDS")
contig_data <- left_join(metadata, contig_data, by = c("ID"= "sample", "Location"="country"))
contig_data <- contig_data[contig_data$Visit_Number == 1,]
contig_data$Health <- as.character(contig_data$Health)
contig_data <- contig_data[!is.na(contig_data$circular),]
contig_data$Health[contig_data$Health == "unhealthy"] <- "rheumatoid arthritis"
contig_data$Location_Health <- paste(contig_data$Location, "-", contig_data$Health)
contig_data <- contig_data %>% group_by(Location, sample_type, Health) %>% mutate(num_samples = n_distinct(ID))

args_data <- read.csv("data/all_assemblies_card_90.csv", stringsAsFactors = FALSE)
args_data$ID <- gsub("_card.out", "", gsub(".*/", "", args_data$filename))
args_data$name <- paste0(args_data$ID, "_", args_data$qseqid)

# Combine phage and arg results where only phage does/doesnt contain ARG
phages_args_yn <- left_join(contig_data, args_data, by = c("name", "ID"), suffix = c(".phage", ".arg"))
phages_args_yn$ARG_Present <- NA
phages_args_yn$ARG_Present[is.na(phages_args_yn$ARO.Name)] <- "NO" 
phages_args_yn$ARG_Present[!is.na(phages_args_yn$ARO.Name)] <- "YES"
phages_args_yn <- phages_args_yn[!duplicated(paste(phages_args_yn$ID, phages_args_yn$name)),]

# Proportion that are circular
sum(phages_args_yn$circular)/nrow(phages_args_yn)*100

phages_args_yn$circular <- ifelse(phages_args_yn$circular, "Circular", "Linear")

tiff("figures/Figure6a.tiff", width = 500, height = 500, res = 100)
ggplot(phages_args_yn, aes(ARG_Present, size/1000)) +
  geom_boxplot() +
  theme_classic() + xlab("Contains ARG") + ylab("Contig Size (kb)") +
  scale_y_continuous(breaks = seq(0, 1200, 100)) +
  theme(panel.grid.major.y = element_line(colour = "grey90")) +
  facet_grid(~ circular)
dev.off()

# Proportion of clusters that contain ARG per Location, sample type
phages_args_summary <- phages_args_yn %>% group_by(Location, sample_type, Health, ARG_Present) %>%
  filter(ARG_Present == "YES") %>%
  summarise(prop_pol = n_distinct(ID)/unique(num_samples)*100)

# Plot proportion of population
tiff("figures/Figure6b.tiff", width = 800, height = 500, res = 100)
ggplot(phages_args_summary, aes(sample_type, prop_pol)) +
  geom_bar(stat="identity") +
  theme_classic() + xlab("Body Site") + ylab("% samples") +
  facet_grid(~Location+Health, scale="free", space="free")
dev.off()

# Combine phage and arg results where only phage contains an ARG
phages_args <- inner_join(contig_data, args_data, by = c("name", "ID"), suffix = c(".phage", ".arg"))
nrow(phages_args)/nrow(contig_data)*100
args <- phages_args %>% 
  select(Location, sample_type, Health, num_samples, ARO.Name, ID, vcontact_cluster) %>% 
  rename(arg_phage = ARO.Name)

# Add classes
args_class <- args_data %>% select(ARO.Name, Drug.Class)
args_class <- args_class[!(duplicated(args_class$ARO.Name)),]
args <- left_join(args, args_class, by = c("arg_phage" = "ARO.Name"))
args$arg_phage <- factor(args$arg_phage, levels = unique(args$arg_phage))
args$Drug.Class[sapply(args$Drug.Class, function(x) str_count(x, ";")) > 2] <- "multidrug class"

# Collapse ARGs for each phage
args <- args %>% group_by(Location, sample_type, Health, num_samples, ID, vcontact_cluster) %>% 
  summarise(arg_phage = paste(unique(arg_phage), collapse="/"), Drug.Class = paste(unique(Drug.Class), collapse="/"))

phage <- phages_args %>% ungroup() %>% select(Location, sample_type, Health, num_samples, vcontact_cluster) %>% 
  distinct() %>% left_join(contig_data, by = c("Location", "sample_type", "Health", "num_samples", "vcontact_cluster")) %>%
  select(Location, sample_type, Health, num_samples, ID, vcontact_cluster) %>%
  mutate(arg_phage = vcontact_cluster, Drug.Class = "Phage cluster (No ARG Class)") 
phages_args <- bind_rows(args, phage)

# Proportion of samples that contain phage and ARG
phages_args_pairs <- phages_args %>%
  group_by(Location, sample_type, Health, num_samples, vcontact_cluster, arg_phage, Drug.Class) %>%
  summarise(n_ID = n_distinct(ID)) %>%
  mutate(percentage = n_ID/num_samples * 100) %>%
  filter(!is.na(arg_phage) & !is.na(vcontact_cluster))

# Relabel vcontact cluster and ARGs
phages_args_labels <- phages_args_pairs %>% ungroup() %>% 
  select(Location, sample_type, Health, vcontact_cluster, arg_phage) %>%
  filter(vcontact_cluster != arg_phage) %>%
  rename(arg = arg_phage)
phages_phage_pairs <- phages_args_pairs %>% filter(Drug.Class == "Phage cluster (No ARG Class)") %>%
  left_join(phages_args_labels, by = c("Location", "sample_type", "Health", "vcontact_cluster"))
arg_pairs <- phages_args_pairs %>% filter(Drug.Class != "Phage cluster (No ARG Class)") %>%
  mutate(arg = arg_phage)
phages_args_pairs <- bind_rows(phages_phage_pairs, arg_pairs)
phages_args_pairs$vcontact_cluster_arg <- paste(phages_args_pairs$vcontact_cluster, "-", phages_args_pairs$arg)

# Relabel
phages_args_pairs$arg_phage_label <- phages_args_pairs$arg_phage
phages_args_pairs$arg_phage_label[phages_args_pairs$Drug.Class == "Phage cluster (No ARG Class)"] <- NA
phages_args_pairs$Location_sampletype <- paste(phages_args_pairs$Location, "-", phages_args_pairs$sample_type)
unique_location_sampletype <- unique(phages_args_pairs$Location_sampletype)
phages_args_pairs$arg_phage_group <- ifelse(phages_args_pairs$Drug.Class == "Phage cluster (No ARG Class)", "Phage", "ARG")

g_phage_arg <- list()
for (i in 1:length(unique_location_sampletype)) {
  g_phage_arg[[i]] <- ggplot(phages_args_pairs[phages_args_pairs$Location_sampletype == unique_location_sampletype[i],], 
         aes(vcontact_cluster_arg, percentage, fill = arg_phage_group, label = arg_phage_label)) +
    geom_bar(stat = "identity", position = "identity", width = 1) +
    geom_text(stat = "summary", fun.y = max, hjust = -0.1, size = 2) +
    coord_flip() + xlab("MAP") +
    theme_classic() + 
    scale_fill_manual("Type", values = c("coral2", "grey80")) + ylab("% samples") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 125, 10)) +
    theme(axis.text.x = element_text(angle=70, hjust = 1),
          strip.background = element_blank(), strip.text.y = element_text(angle = 180),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    facet_grid(vcontact_cluster~., scales = "free", space = "free", switch = "y") +
    ggtitle(unique_location_sampletype[i])
}

# Plot phage-arg
tiff("figures/Supplementary_Figure4.tiff", width = 3000, height = 4000, res = 200)
ggplot(phages_args_pairs, aes(vcontact_cluster_arg, percentage, fill = arg_phage_group, label = arg_phage_label)) +
  geom_bar(stat = "identity", position = "identity", width = 1) +
  geom_text(stat = "summary", fun.y = max, hjust = -0.1, size = 2) +
  coord_flip() + xlab("MAP") +
  theme_classic() + 
  scale_fill_manual("Type", values = c("coral2", "grey80")) + ylab("% samples") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 125, 10)) +
  theme(axis.text.x = element_text(angle=70, hjust = 1),
        strip.background = element_blank(), strip.text.y = element_text(angle = 180),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_grid(Location + Health + sample_type + vcontact_cluster~., scales = "free", space = "free", switch = "y")
dev.off()

########## Megaphages and ARGs ####
# Read data
metaphage_arg_blast <- read.table("data/megaphage_proteins_card.out", stringsAsFactors = FALSE)
column_headers <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
names(metaphage_arg_blast) <- column_headers

# megaphage_contigs <- read.table("data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
# megaphage_contigs <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
# megaphage_contigs <- megaphage_contigs[megaphage_contigs$Visit_Number == 1,]
# megaphage_contigs$Health <- as.character(megaphage_contigs$Health)
# megaphage_contigs <- megaphage_contigs[!is.na(megaphage_contigs$circular),]
# megaphage_contigs$Health[megaphage_contigs$Health == "unhealthy"] <- "rheumatoid arthritis"
# megaphage_contigs$Location_Health <- paste(megaphage_contigs$Location, "-", megaphage_contigs$Health)
# megaphage_contigs <- megaphage_contigs %>% group_by(Location, sample_type, Health) %>% mutate(num_samples = n_distinct(ID))
# args_data <- read.csv("data/all_assemblies_card_90.csv", stringsAsFactors = FALSE)
# args_data$ID <- gsub("_card.out", "", gsub(".*/", "", args_data$filename))
# args_data$name <- paste0(args_data$ID, "_", args_data$qseqid)
# 
# # Combine phage and arg results where only phage contains an ARG
# megaphages_args <- inner_join(megaphage_contigs, args_data, by = c("name", "ID"), suffix = c(".phage", ".arg"))
# nrow(phages_args)/nrow(contig_data)*100
# args <- phages_args %>% 
#   select(Location, sample_type, Health, num_samples, ARO.Name, ID, vcontact_cluster) %>% 
#   rename(arg_phage = ARO.Name)

########## Megaphages and CRISPR hosts ####
megaphage_contigs <- read.delim("megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
megaphage_contigs <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
megaphage_contigs$Health <- as.character(megaphage_contigs$Health)
megaphage_contigs$Health[megaphage_contigs$Health == "unhealthy"] <- "rheumatoid arthritis"
megaphage_contigs <- megaphage_contigs[megaphage_contigs$Visit_Number == 1,]
megaphage_contigs <- megaphage_contigs[!is.na(megaphage_contigs$circular),]
megaphage_contigs <- megaphage_contigs[megaphage_contigs$circular,]

megaphages_crispr_summary <- megaphage_contigs %>%
  group_by(sample_type, Location, Health, crispr_host) %>%
  summarise(no_crispr_hosts = n())

megaphages_crispr_summary$crispr_host[is.na(megaphages_crispr_summary$crispr_host)] <- "unclassified"

crispr_names <- unique(megaphages_crispr_summary$crispr_host[megaphages_crispr_summary$crispr_host != "unclassified"])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
crispr_colours = c(rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
crispr_colours <- crispr_colours[c(1:length(crispr_names))]
crispr_colours <- c(crispr_colours, "white")
names(crispr_colours) <- c(crispr_names, "unclassified")

tiff("figures/megaphages_crispr_host.tiff", width = 2200, height = 2000, res = 200)
ggplot(megaphages_crispr_summary, aes(sample_type, no_crispr_hosts, fill = crispr_host)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  facet_grid(~ Location + Health, space = "free", scale = "free") +
  theme_classic() +
  ylab("Number of megaphages") + xlab("Body Site") +
  scale_fill_manual("Predicted host", values = crispr_colours)
dev.off()

########## Eggnog analysis ####
megaphage_prots <- read.delim("data/megaphage_proteins_headers.txt", stringsAsFactors = FALSE, header = FALSE)
names(megaphage_prots) <- "query_name"
eggnog <- read.delim("data/megaphage_annot.emapper.annotations", stringsAsFactors = FALSE, header = FALSE)
names(eggnog) <- c("query_name", "seed_eggnog_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "predicted_taxonomic_group",
                   "predicted_protein_name", "gene_ontology_terms", "EC_number", "KEGG_ko", "KEGG_pathway", "KEGG_module", "KEGG_reaction",
                   "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_reaction", "tax_scope", "eggnog_ogs", "bestOG", "COG_functional_cat", "eggnog_description")
megaphage_contigs <- read.delim("data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)

megaphage_eggnog <- left_join(megaphage_prots, eggnog, by = "query_name")
megaphage_eggnog$name <- sub("_[^_]+$", "", megaphage_eggnog$query_name)
megaphage_eggnog <- left_join(megaphage_eggnog, megaphage_contigs, by = "name")
megaphage_eggnog <- left_join(megaphage_eggnog, metadata, by = c("sample"="ID"))

# Remove linear contigs
megaphage_eggnog <- megaphage_eggnog[megaphage_eggnog$circular,]

# Summarise proportion of bacterial/viral proteins
megaphage_egng_summary <- megaphage_eggnog %>% group_by(name, tax_scope, size, Location, sample_type, Health) %>%
  summarise(tax_n = n()) %>% group_by(name, size, Location, sample_type, Health) %>%
  mutate(total_n = sum(tax_n)) %>%
  mutate(tax_perc = tax_n/total_n * 100)
megaphage_func <- megaphage_egng_summary[!is.na(megaphage_egng_summary$tax_scope),]
megaphage_nonviral <- megaphage_func[megaphage_func$tax_scope != "Viruses",]

# Some megaphages have a lot of bacterial genes - as much as 43% of genes bacterial

# Look at largest phage
large_phage <- megaphage_eggnog[which.max(megaphage_eggnog$size),]

# Select non-annorated proteins
proteins_nanno <- megaphage_eggnog$query_name[is.na(megaphage_eggnog$seed_eggnog_ortholog)]
write.table(proteins_nanno, "data/megaphage_proteins_no_annot.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
