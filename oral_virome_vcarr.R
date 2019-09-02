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

set.seed(1)

########Read in and process raw data######################
#Metadata
metadata = read.csv("metadata_andrey.csv")
rownames(metadata) = metadata$ID
merged_samples = strsplit(na.omit(as.character(metadata$Merged_with)), split = '|', fixed = TRUE)

#Contigs
contig_ids = read.table("virsorter_positive.ids")
contig_ids = unique(as.character(contig_ids$V1))

#Read alignment counts
vir_counts = read.table("bowtie2_read_counts.txt", header=TRUE, row.names = 1, sep = "\t", check.names = FALSE)
#vir_counts = fread("bowtie2_read_counts.txt", header=TRUE, sep = "\t", check.names = FALSE)
vir_coverage = read.table("breadth_cov_collated.tbl", header=TRUE, row.names = 1, check.names = FALSE)
#vir_coverage = fread("breadth_cov_collated.tbl", header=TRUE, sep = "\t", check.names = FALSE)
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

#Remove samples with missing metadata
vir_counts_prop = vir_counts_prop[,which(colnames(vir_counts_prop) %in% rownames(metadata))]

#Remove rows with all zeroes
vir_counts_prop = vir_counts_prop[-which(apply(vir_counts_prop, 1, function(x) all(x == 0))),]

########Read in and process contig annotations############

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

circular = read.table("circular_contigs.ids", header=FALSE)
vir_refseq = read.table("refseq_annot.txt", sep = "\t", header=FALSE)
crasslike = read.table("crasslike_annot.txt", sep = "\t", header=FALSE)
demovir = read.table("DemoVir_assignments.txt", sep = "\t", header=TRUE, row.names = 1)
gc_content = read.table("virsorter_positive_gc_content.txt", header=FALSE, row.names = 1)
ribo_prot_count = read.table("ribosomal_prot_counts.txt", header=FALSE)
rownames(ribo_prot_count) = as.character(ribo_prot_count$V2)
pvogs_count = read.table("pvogs_counts_per_contig.txt", header=FALSE)
rownames(pvogs_count) = as.character(pvogs_count$V2)
vcontact = read.table("new_clusters_vcontact_w_taxa.txt", header=TRUE, sep = "\t")
rownames(vcontact) = as.character(vcontact$contig_ID)
integrase = read.table("integrase_contigs.ids", header=FALSE)
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

country = read.table("sample_countries.txt")
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

#######Host information from CRISPR and tRNA matches#################
#Read in and process CRISPR BLAST data
crispr_bitscore_cutoff = 45
trna_bitscore_cutoff = 65
crispr_evalue_cutoff = 1E-05
trna_evalue_cutoff = 1E-10

#crispr_hits1 =
#  read.table(file = 'refseq89_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)
crispr_hits2 =
  read.table(file = 'hmp_pilercr_blasthits_selected_contigs_taxa.txt', header = FALSE, sep = '\t', quote="", fill=FALSE)

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

contig_data$crispr_host = NA
for(i in unique(as.character(crispr_hits$sseqid))) {
  temp = crispr_hits[crispr_hits$sseqid == i,c("genus", "bitscore")]
  temp = aggregate(bitscore ~ genus, data = temp, sum)
  contig_data[i,"crispr_host"] = temp[which.max(temp$bitscore),"genus"]
}

##########Taxonomy from co-clustering with RefSeq######################
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
contig_data = contig_data[rownames(vir_counts_prop),]

#Cleanup
rm(integrase,circular,vir_refseq,crasslike,demovir,gc_content,ribo_prot_count,
   pvogs_count,vcontact,crispr_hits2)
saveRDS(contig_data,  file = "contig_data.RDS")
write.table(contig_data, file = "contig_data.txt")

#########Contig scatterplot####################################
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

#contig_scatter_demovir = ggMarginal(contig_scatter_demovir, size = 10, type = "density", fill = "grey") 

#########Virome composition barplots###########################
#Make long format table
vir_counts_prop_melt = melt(vir_counts_prop)
vir_counts_prop_melt$demovir = contig_data[as.character(vir_counts_prop_melt$Var1), "demovir"]
vir_counts_prop_melt$vcontact_cluster = contig_data[as.character(vir_counts_prop_melt$Var1), "vcontact_cluster"]
vir_counts_prop_melt$crispr_host = contig_data[as.character(vir_counts_prop_melt$Var1), "crispr_host"]
vir_counts_prop_melt$integrase = contig_data[as.character(vir_counts_prop_melt$Var1), "integrase"]
vir_counts_prop_melt$sample_type = metadata[as.character(vir_counts_prop_melt$Var2), "sample_type"]
vir_counts_prop_melt$Location = metadata[as.character(vir_counts_prop_melt$Var2), "Location"]

#Aggregate counts by Demovir/vConTACT2
vir_counts_prop_melt = as.data.table(vir_counts_prop_melt)
vir_counts_prop_melt_agg = vir_counts_prop_melt[, sum(value), by=.(Var2, demovir, vcontact_cluster, sample_type, Location)]
saveRDS(vir_counts_prop_melt_agg,  file = "vir_counts_prop_melt_agg.RDS")

#Aggregate counts by CRISPR host
vir_counts_prop_melt_agg2 = vir_counts_prop_melt[, sum(value), by=.(Var2, crispr_host, sample_type, Location)]
vir_counts_prop_melt_agg2$crispr_host = as.factor(vir_counts_prop_melt_agg2$crispr_host)
saveRDS(vir_counts_prop_melt_agg2, file = "vir_counts_prop_melt_agg2.RDS")
rm(vir_counts_prop_melt)

barplot_demovir = ggplot(vir_counts_prop_melt_agg, aes(x = Var2, y = V1, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(brewer.pal(12, "Paired")[c(1,2,4,5,7,8,9,10,12)], "lightgray")) +
  facet_wrap(~sample_type, scales = "free_x")

#bar_colors_vir =
#  colorRampPalette(c(brewer.pal(8, "Spectral"),brewer.pal(8, "RdBu"), brewer.pal(8, "Dark2")))(10000)
#barplot_vcontact = ggplot(vir_counts_prop_melt_agg, aes(x = Var2, y = V1, fill = vcontact_cluster)) + theme_classic() +
#  geom_bar(stat = "identity", show.legend = FALSE) +
#  scale_fill_manual(values = bar_colors_vir[sample(1:10000,length(levels(vir_counts_prop_melt_agg$vcontact_cluster)))]) +
#  facet_wrap(~sample_type, scales = "free_x")

#########Virome beta-diversity###########################
#Re-cast counts matrix by vConTACT2 clusters
vir_counts_prop_agg = dcast.data.table(vir_counts_prop_melt_agg, vcontact_cluster ~ Var2, value.var = "V1", fun.aggregate = sum)  
vir_counts_prop_agg = as.data.frame(vir_counts_prop_agg)
vir_counts_prop_agg = vir_counts_prop_agg[-which(is.na(vir_counts_prop_agg[,1])),]
rownames(vir_counts_prop_agg) = vir_counts_prop_agg[,1]
vir_counts_prop_agg = vir_counts_prop_agg[,-1]
vir_counts_prop_agg = as.matrix(vir_counts_prop_agg)

saveRDS(vir_counts_prop_agg,  file = "vir_counts_prop_agg.RDS")

#Compute distance matrix and run t-SNE
cluster_counts_dist = vegdist(t(vir_counts_prop_agg), method = "bray")
clusters_samples_tsne = tsne(cluster_counts_dist)

#Annotate t-SNE coordinates
clusters_samples_tsne = as.data.frame(clusters_samples_tsne)
rownames(clusters_samples_tsne) = colnames(vir_counts_prop_agg)
colnames(clusters_samples_tsne) = c("tsne1", "tsne2")
clusters_samples_tsne$country = metadata[rownames(clusters_samples_tsne),"Location"]
clusters_samples_tsne$body_site = metadata[rownames(clusters_samples_tsne),"sample_type"]
clusters_samples_tsne$subject = metadata[rownames(clusters_samples_tsne),"Sample.name"]
clusters_samples_tsne$sex = metadata[rownames(clusters_samples_tsne),"Gender"]
clusters_samples_tsne$Age = metadata[rownames(clusters_samples_tsne),"Age"]

#De novo clustering
samples_clust = pam(cluster_counts_dist, 7)

clusters_samples_tsne$cluster = as.factor(samples_clust$cluster[rownames(clusters_samples_tsne)])
clusters_samples_tsne$n_contigs = sapply(rownames(clusters_samples_tsne), function(x) length(which(contig_data$sample == x)))
clusters_samples_tsne$total_reads = as.numeric(counts_total[rownames(clusters_samples_tsne)])

saveRDS(vir_counts_prop_agg,  file = "vir_counts_prop_agg.RDS")
saveRDS(clusters_samples_tsne,  file = "clusters_samples_tsne.RDS")

#Make t-SNE plots
tsne_plot1 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = cluster, shape = country)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = brewer.pal(7, "Set1"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 3)))

tsne_plot2 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = n_contigs, size = log10(total_reads), shape = country)) +
  theme_classic() + geom_point(alpha = 0.7) +
  scale_shape_manual(values = c(21, 24, 22))

tsne_plot3 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = body_site, shape = country)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_brewer(palette = "Paired",
                    guide = guide_legend(override.aes = list(shape = 21, size = 3)))

tsne_plot4 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = Age, shape = sex)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  #scale_fill_gradient2()
  scale_fill_distiller(palette = "Spectral")

#########Differentially abundant clusters#####################

clusters_samples_kw = apply(vir_counts_prop_agg, 1, function(x) kruskal.test(x, clusters_samples_tsne$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
vir_counts_prop_agg_diff = vir_counts_prop_agg[names(which(clusters_samples_kw < 0.001)),]

viral_clusters_df = data.frame(row.names = rownames(vir_counts_prop_agg),
                               demovir = sapply(rownames(vir_counts_prop_agg), function(x) names(which.max(table(contig_data[which(contig_data$vcontact_cluster == x),"demovir"]))))
)

########################Host range############################
#Re-cast counts matrix by CRISPR hosts
vir_counts_prop_agg2 = dcast.data.table(vir_counts_prop_melt_agg2, crispr_host ~ Var2, value.var = "V1", fun.aggregate = sum)  
vir_counts_prop_agg2 = as.data.frame(vir_counts_prop_agg2)
vir_counts_prop_agg2 = vir_counts_prop_agg2[-which(is.na(vir_counts_prop_agg2[,1])),]
rownames(vir_counts_prop_agg2) = vir_counts_prop_agg2[,1]
vir_counts_prop_agg2 = vir_counts_prop_agg2[,-1]
vir_counts_prop_agg2 = as.matrix(vir_counts_prop_agg2)

saveRDS(vir_counts_prop_agg2,  file = "vir_counts_prop_agg2.RDS")

#########Saving output#######################################
pdf(file="contigs_scatter.pdf", paper="A4r", width=11, height=8.5)
contig_scatter_demovir
dev.off()

pdf(file="samples_barplot.pdf", paper="A4r", width=11, height=8.5)
barplot_demovir
dev.off()

pdf(file="samples_tsne.pdf", paper="A4r", width=11, height=8.5)
grid.arrange(tsne_plot1, tsne_plot2, tsne_plot3, tsne_plot4, nrow = 2)
dev.off()

#Clusters heatmap
vir_counts_prop_agg_diff_log = log10(vir_counts_prop_agg_diff)
vir_counts_prop_agg_diff_log[is.infinite(vir_counts_prop_agg_diff_log)] = -8

pdf(file = "heatplot_clusters.pdf", width=8.5, height=11, pointsize = 12)
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

pdf(file = "heatplot_hosts.pdf", width=8.5, height=11, pointsize = 12)
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