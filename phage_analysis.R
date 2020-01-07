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
library(devtools)
library(circlize)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#BiocManager::install("maSigPro")
library(maSigPro)
library(purrr)
library(ggrepel)
library(stringr)
library(tidyr)
library(gtable)
library(igraph)
library(ComplexHeatmap)
set.seed(1)

########## Read data, metadata and define parameters ####
# Data
vir_counts_prop_melt <- readRDS("data/vir_counts_prop_melt.RDS")
vir_counts_prop_agg <- readRDS("data/vir_counts_prop_agg.RDS")
vir_counts_prop_melt_agg <- readRDS("data/vir_counts_prop_melt_agg.RDS")
vir_counts_prop_melt_agg2 <- readRDS("data/vir_counts_prop_melt_agg2.RDS")
contig_data <- readRDS("data/contig_data.RDS")
counts_total <- readRDS("data/counts_total.RDS")
megaphage_contigs <- read.delim("data/megaphage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
metaphlan <- read.csv("data/all_metaphlan.csv", stringsAsFactors = FALSE)

# Metadata
metadata = read.csv("data/metadata_v2.csv", stringsAsFactors = FALSE)
rownames(metadata) = metadata$ID
metadata$Location_sampletype <- paste(metadata$Location, "-", metadata$sample_type)

# n summary of metadata excluding longitudinal samples
metadata_summary <- metadata %>% filter(Visit_Number == 1) %>% group_by(Location, sample_type) %>% summarise(n())

# Body site colours
cols <- brewer.pal(11, "Spectral")[c(2, 4, 9, 10, 11)]
names(cols) <- unique(metadata$sample_type)

########## Virome composition barplots ####
demovir_cols <- rev(c(brewer.pal(length(unique(vir_counts_prop_melt_agg$demovir)), "Set3")))
names(demovir_cols) <- sort(as.character(unique(vir_counts_prop_melt_agg$demovir)))
demovir_cols[names(demovir_cols) == "crAss-like"] <- "grey50"
demovir_cols[names(demovir_cols) == "Unassigned"] <- "grey80"
demovir_cols[demovir_cols == "#FFFFB3"] <- "#8DD3C7"

tiff("figures/barplot_demovir.tiff", width = 5000, height = 2000, res = 400)
ggplot(vir_counts_prop_melt_agg, aes(x = Var2, y = V1, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab("Proportion of mapped reads") + xlab("Sample") + ylim(0, max(aggregate(vir_counts_prop_melt_agg$V1, by=list(ID=vir_counts_prop_melt_agg$Var2), FUN=sum)$x))
dev.off()

########## Virome beta-diversity###########################
# Compute distance matrix and run t-SNE
vir_counts_prop_agg_meta <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% metadata$ID]
vir_counts_prop_agg_meta <- vir_counts_prop_agg_meta[,colnames(vir_counts_prop_agg_meta) %in% metadata$ID[metadata$Visit_Number == 1]]

# T-sne
set.seed(1)
cluster_counts_dist = vegdist(t(vir_counts_prop_agg_meta), method = "bray")
clusters_samples_tsne = tsne(cluster_counts_dist)

#Annotate t-SNE coordinates
clusters_samples_tsne = as.data.frame(clusters_samples_tsne)
rownames(clusters_samples_tsne) = colnames(vir_counts_prop_agg_meta)
clusters_samples_tsne$ID <- rownames(clusters_samples_tsne)
colnames(clusters_samples_tsne) = c("tsne1", "tsne2", "ID")
clusters_samples_tsne <- left_join(clusters_samples_tsne, metadata, by = "ID")

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 2:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(clusters_samples_tsne[,c("tsne1", "tsne2")], k = k), clusters_samples_tsne[,c("tsne1", "tsne2")])
  avg_sil[k-1] <- mean(tmp[,3])
}
# Group by silhouette width
samples_clust <- pam(cluster_counts_dist, which.max(avg_sil)+1)

clusters_samples_tsne$cluster = as.factor(samples_clust$cluster[clusters_samples_tsne$ID])
clusters_samples_tsne$n_contigs <- sapply(clusters_samples_tsne$ID, function(x) length(which(vir_counts_prop_melt$Var2 == x)))
clusters_samples_tsne$total_reads = as.numeric(counts_total[clusters_samples_tsne$ID])

saveRDS(clusters_samples_tsne,  file = "data/clusters_samples_tsne.RDS")

# Make t-SNE plots
# Label by cluster
tsne_plot1 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = cluster)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7, pch = 21) +
  scale_fill_manual("Group", values = brewer.pal(8, "Set2"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("Dim 1") + ylab("Dim 2")

# Label by contig size
tsne_plot2 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = n_contigs, size = log10(total_reads))) +
  theme_classic() + geom_point(alpha = 0.7, pch = 21) +
  scale_fill_continuous("No. contigs") + scale_size_continuous("Log10(No. reads)") +
  xlab("Dim 1") + ylab("Dim 2")

# Label by body sites, geographical location and health
tsne_plot3 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = sample_type, shape = Location)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7) +
  scale_shape_manual("Country", values = c(21, 24, 22)) +
  scale_fill_manual("Body Site", values = cols,
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("Dim 1") + ylab("Dim 2")

# Calculate proportion of samples
cluster_res <- clusters_samples_tsne %>% group_by(cluster, Location, sample_type) %>% summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop_cluster = n/total_n*100, Location_sampletype = paste(sample_type, "-", Location))

# Label by sex and age
tsne_plot4 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = Age, shape = Gender)) +
  theme_classic() + geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual("Sex", values = c(21, 22)) +
  #scale_fill_gradient2()
  scale_fill_distiller(palette = "Spectral") +
  xlab("Dim 1") + ylab("Dim 2")

# Save plots
tiff("figures/tsne_clusters.tiff", width = 700, height = 500, res = 150)
tsne_plot1
dev.off()

cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "yellowgreen", brewer.pal(9, "YlOrRd")[c(5,7)], brewer.pal(9, "RdPu")[c(3,5)])
tiff("figures/tsne_bodysite_location.tiff", width = 1000, height = 500, res = 150)
ggplot(cluster_res, aes(cluster, prop_cluster, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  theme_classic() + xlab("Group") + ylab("Percentage") +
  scale_fill_manual("Body Site - Geographical Location", values = cohort_cols)
dev.off()

# Set the same legend widths
g2 <- ggplotGrob(tsne_plot2)
g4 <- ggplotGrob(tsne_plot4)
g3 <- ggplotGrob(tsne_plot3)
g4$widths <- g2$widths
g3$widths <- g2$widths

tiff("figures/tsne_contigs_bodysite_age_gender.tiff", width = 2400, height = 500, res = 130)
grid.arrange(g3, g2, g4, nrow = 1)
dev.off()

########## GLM of variables###########
glm_samples <- clusters_samples_tsne[!is.na(clusters_samples_tsne$Age),]
glm(cluster ~ Location + sample_type + Gender + Age, family = "binomial", data = glm_samples)
# Difficult as variables colinear

########## Differentially abundant clusters#####################
clusters_samples_kw = apply(vir_counts_prop_agg_meta, 1, function(x) kruskal.test(x, clusters_samples_tsne$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
vir_counts_prop_agg_meta_diff = vir_counts_prop_agg_meta[rownames(vir_counts_prop_agg_meta) %in% names(which(clusters_samples_kw < 0.001)),]

viral_clusters_df = data.frame(row.names = rownames(vir_counts_prop_agg_meta),
                               demovir = sapply(rownames(vir_counts_prop_agg_meta), 
                                                function(x) names(which.max(table(contig_data[which(contig_data$vcontact_cluster == x),"demovir"]))))
)

#Clusters heatmap
vir_counts_prop_agg_meta_diff_log = log10(vir_counts_prop_agg_meta_diff)
vir_counts_prop_agg_meta_diff_log[is.infinite(vir_counts_prop_agg_meta_diff_log)] = -8

tiff(file = "figures/heatplot_clusters.tiff", width=1500, height=2000, res = 150)
heatmap.2(vir_counts_prop_agg_meta_diff_log,
          margins = c(10,10),
          trace = "none",
          scale = "none",
          density.info = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(min(vir_counts_prop_agg_meta_diff_log), 0, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA,
          key.xlab = "log10(prop. reads mapped)",
          key.ylab = NA,
          RowSideColors = unlist(sapply(viral_clusters_df[rownames(viral_clusters_df) %in% rownames(vir_counts_prop_agg_meta_diff),"demovir"], function(x) demovir_cols[names(demovir_cols) == x])),
          ColSideColors = cols[as.factor(metadata[colnames(vir_counts_prop_agg_meta_diff_log), "sample_type"])],
          labRow = NA,
          #labRow = gsub("NULL", "", as.character(sapply(clusters_tax[rownames(vir_counts_prop_agg_meta_diff)], "[[", 1))),
          labCol = NA,
          cexCol = 0.5,
          cexRow = 0.5,
          xlab = "Samples",
          ylab = "Contig clusters",
          main = NA
)
legend("topright", legend = levels(factor(metadata[colnames(vir_counts_prop_agg_meta_diff), "sample_type"])),
       col = cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Sample clusters:")

legend(x = -0.05, y = 0.95, xpd=TRUE, legend = levels(factor(viral_clusters_df[rownames(viral_clusters_df) %in% rownames(vir_counts_prop_agg_meta_diff),"demovir"], levels = names(demovir_cols))),
       col = demovir_cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.4, title = "Taxonomy by Demovir")

dev.off()

########## Venn diagram showing overview of total phages for each group ##########
unique_clusters <- unique(clusters_samples_tsne$cluster)

vir_group_list <- sapply(unique_clusters, function(x) {
  vir_tmp <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% clusters_samples_tsne$ID[clusters_samples_tsne$cluster == x]]
  vir_group <- rownames(vir_tmp)[rowSums(vir_tmp) != 0]
  return(vir_group)
  })
names(vir_group_list) <- paste("Group", unique_clusters)
max_group_length <- max(sapply(vir_group_list, function(x) length(x)))

vir_group_list <- lapply(vir_group_list, function(x) c(x, rep(NA, max_group_length - length(x))))

vir_group_all <- do.call(cbind, vir_group_list)

tiff("figures/venn_diagram.tiff", width = 600, height = 600)
suma2Venn(vir_group_all, zcolor = brewer.pal(length(unique_clusters), "Set2"), cexil = 1.5, cexsn = 1.4)
dev.off()

########## Longitudinal analysis and core phageome #########
# Create mock IDs so different time points
metadata_us <- metadata[metadata$Location == "US",]
metadata_longus <- metadata_us %>% group_by(Sample.name, sample_type) %>%
  filter(any(Visit_Number == 3))

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

# Count average time points of phages
contig_no_tp <- vir_counts_longus %>% group_by(sample_type, Sample.name, Var1) %>%
  summarise(no_tp = n_distinct(Visit_Number)) 

# Count all distinct contigs
contig_all_summary <- vir_counts_longus %>% group_by(sample_type, Sample.name) %>%
  summarise(total_contig = n_distinct(Var1))

# Count contigs in at least x number of time points
contig_count_tp <- map_df(.x = unique(contig_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(contig_no_tp$sample_type[contig_no_tp$no_tp == .x], .x),
                    Sample.name = rep(contig_no_tp$Sample.name[contig_no_tp$no_tp == .x], .x),
                    Var1 = rep(contig_no_tp$Var1[contig_no_tp$no_tp == .x]),
                    no_tp = rep(contig_no_tp$no_tp[contig_no_tp$no_tp == .x], .x))
  tmp$cum_tp <- c(1:.x)
  return(tmp)
})

contig_count_tp_summary <- contig_count_tp %>% group_by(sample_type, Sample.name, cum_tp) %>%
  summarise(total_contig_least = n_distinct(Var1))

contig_one <- contig_count_tp_summary %>%
  filter(cum_tp == 1) %>%
  rename(total_contig = total_contig_least) %>%
  select(-cum_tp)

contig_count_tp_summary <- left_join(contig_count_tp_summary, contig_one, by = c("sample_type", "Sample.name")) %>%
  mutate(contig_frac = total_contig_least/total_contig)

# Plot
tiff("figures/longitudinal_contig_count.tiff", height = 400, width = 1100, res = 100)
set.seed(1)
ggplot(contig_count_tp_summary, aes(as.factor(cum_tp), contig_frac, fill = sample_type)) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of contigs in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = cols)
dev.off()

# Proportion of reads mapped to contigs
contig_count_reads<- vir_counts_longus %>% 
  group_by(sample_type, Sample.name, Var1) %>%
  arrange(desc(Visit_Number)) %>%
  mutate(cumsum_reads = cumsum(value)) %>%
  select(sample_type, Sample.name, Var1, cumsum_reads, Visit_Number)

contig_count_reads_tp1 <- contig_count_reads %>%
  filter(Visit_Number == min(Visit_Number)) %>%
  rename(tp1 = cumsum_reads) %>%
  select(sample_type, Sample.name, Var1, tp1)

contig_count_reads_summary <- left_join(contig_count_reads, contig_count_reads_tp1, by = c("sample_type", "Var1", "Sample.name")) %>%
  mutate(frac_reads = cumsum_reads/tp1) %>%
  group_by(sample_type, Sample.name, Var1) %>%
  mutate(cum_tp = order(Visit_Number))

tiff("figures/longitudinal_read_count.tiff", height = 400, width = 1100, res = 90)
ggplot(contig_count_reads_summary, aes(x = as.factor(cum_tp), y = frac_reads, fill = sample_type)) +
  geom_boxplot() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of reads mapped in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = cols)
dev.off()

# Select persistent contigs (3 at least, 2 doesnt converge)
contig_persist <- left_join(vir_counts_longus, contig_no_tp, by = c("sample_type", "Sample.name", "Var1")) %>%
  filter(no_tp >= 3)

# Summarise n
persist_summary <- metadata %>% filter(ID %in% unique(contig_persist$ID)) %>% group_by(sample_type, Visit_Number) %>% summarise(n())

# Get phage taxonomy
tiff("figures/barplot_demovir_persistent.tiff", width = 5000, height = 2000, res = 400)
ggplot(contig_persist, aes(x = ID, y = value, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab("Proportion of mapped reads") + xlab("Sample") + ylim(0, max(aggregate(contig_persist$value, by=list(ID=contig_persist$ID), FUN=sum)$x))
dev.off()

# Cast for NMDS
contig_persist_cast <- dcast(contig_persist, ID ~ Var1, sum, value.var = "value") 
rownames(contig_persist_cast) <- contig_persist_cast$ID
contig_persist_cast <- contig_persist_cast[,names(contig_persist_cast) != "ID"]

# Run NMDS
set.seed(1)
nmds_persist <- metaMDS(contig_persist_cast, distance = "bray", k = 2, trymax = 20)
df_nmds_persist <- as.data.frame(nmds_persist$points)
df_nmds_persist$ID <- row.names(df_nmds_persist)
df_nmds_persist$Sample.name <- as.character(sapply(df_nmds_persist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_nmds_persist$Location <- sapply(df_nmds_persist$ID, function(x) metadata$Location[metadata$ID == x])
df_nmds_persist$sample_type <- sapply(df_nmds_persist$ID, function(x) metadata$sample_type[metadata$ID == x])

# Plot NMDS
tiff("figures/persistent_phages.tiff", width = 2000, height = 1000, res = 200)
ggplot(df_nmds_persist, aes(MDS1, MDS2, fill = Sample.name, colour = sample_type)) +
  geom_point() +
  geom_density2d(alpha=0.3) +
  theme_classic() +
  guides(fill = FALSE) +
  scale_colour_manual("body site", values = cols) +
  scale_x_continuous(limits = c(-2,2), breaks = seq(-2,2,0.5)) + 
  scale_y_continuous(limits = c(-1.5,1.2), breaks = seq(-2,2,0.5))
dev.off()

# Select nonpersistent contigs
contig_nonpersist <- left_join(vir_counts_longus, contig_no_tp, by = c("sample_type", "Sample.name", "Var1")) %>%
  filter(no_tp < 3)

# Summarise n
nonpersist_summary <- metadata %>% filter(ID %in% unique(contig_nonpersist$ID)) %>% group_by(sample_type, Visit_Number) %>% summarise(n())

# Get phage taxonomy
tiff("figures/barplot_demovir_transient.tiff", width = 5000, height = 2000, res = 400)
ggplot(contig_nonpersist, aes(x = ID, y = value, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab("Proportion of mapped reads") + xlab("Sample") + ylim(0, max(aggregate(contig_persist$value, by=list(ID=contig_persist$ID), FUN=sum)$x))
dev.off()

# Cast for NMDS
contig_nonpersist_cast <- dcast(contig_nonpersist, ID ~ Var1, sum, value.var = "value") 
rownames(contig_nonpersist_cast) <- contig_nonpersist_cast$ID
contig_nonpersist_cast <- contig_nonpersist_cast[,names(contig_nonpersist_cast) != "ID"]

# Run NMDS
set.seed(1)
nmds_nonpersist <- metaMDS(contig_nonpersist_cast, distance = "bray", k = 2, trymax = 20)
df_nmds_nonpersist <- as.data.frame(nmds_nonpersist$points)
df_nmds_nonpersist$ID <- row.names(df_nmds_nonpersist)
df_nmds_nonpersist$Sample.name <- as.character(sapply(df_nmds_nonpersist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_nmds_nonpersist$Location <- sapply(df_nmds_nonpersist$ID, function(x) metadata$Location[metadata$ID == x])
df_nmds_nonpersist$sample_type <- sapply(df_nmds_nonpersist$ID, function(x) metadata$sample_type[metadata$ID == x])

# Plot NMDS
tiff("figures/transient_phages.tiff", width = 2000, height = 1000, res = 200)
ggplot(df_nmds_nonpersist, aes(MDS1, MDS2, fill = Sample.name, colour = sample_type)) +
  geom_point() +
  geom_density2d(alpha=0.2) +
  theme_classic() +
  guides(fill = FALSE) +
  scale_colour_manual("body site", values = cols) +
  scale_x_continuous(limits = c(-1.5,1), breaks = seq(-2,2,0.5)) + 
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-2,2,0.5))
dev.off()

########## Host range############################
#Re-cast counts matrix by CRISPR hosts
not_genera <- c(NA, "Lachnospiraceae", "Bacteroidetes")
vir_counts_prop_agg2 = dcast.data.table(vir_counts_prop_melt_agg2[vir_counts_prop_melt_agg2$Var2 %in% metadata$ID[metadata$Visit_Number == 1] & 
                                                                                                                    !(vir_counts_prop_melt_agg2$crispr_host %in% not_genera),], 
                                        crispr_host ~ Var2, value.var = "V1", fun.aggregate = sum)
vir_counts_prop_agg2 = as.data.frame(vir_counts_prop_agg2)
vir_counts_prop_agg2 = vir_counts_prop_agg2[!is.na(vir_counts_prop_agg2$crispr_host),]
rownames(vir_counts_prop_agg2) = vir_counts_prop_agg2$crispr_host
vir_counts_prop_agg2 = vir_counts_prop_agg2[,-1]
vir_counts_prop_agg2 = as.matrix(vir_counts_prop_agg2)

#Hosts heatmap
vir_counts_prop_agg2_log = log10(vir_counts_prop_agg2)
vir_counts_prop_agg2_log[is.infinite(vir_counts_prop_agg2_log)] = -8

tiff("figures/heatplot_hosts_genus.tiff", width = 1000, height = 1500, res = 150)
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
          key.title = NA, key.xlab = "log10(prop. reads mapped)", key.ylab = NA,
          ColSideColors = cols[as.factor(metadata[colnames(vir_counts_prop_agg2_log), "sample_type"])],
          #labRow = as.character(rownames(vir_counts_prop_agg2_log)),
          labCol = NA,
          cexCol = 0.5,
          #cexRow = 1,
          xlab = "Samples",
          ylab = "Predicted host genera",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(vir_counts_prop_agg2_log), "sample_type"])),
       col = cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Body sites")
dev.off()

########## Microbial composition #########
row.names(metaphlan) <- metaphlan$X
metaphlan <- metaphlan[,names(metaphlan) != "X"]
metaphlan <- metaphlan[,names(metaphlan) %in% clusters_samples_tsne$ID]
metaphlan <- metaphlan[,names(metaphlan) %in% metadata$ID[metadata$Visit_Number == 1]]
metaphlan_filter <- metaphlan[grepl("s__", row.names(metaphlan)) & !grepl("t__", row.names(metaphlan)),]
metaphlan_filter <- metaphlan_filter[grepl("k__Bacteria", row.names(metaphlan_filter)) | 
                                       grepl("k__Archaea", row.names(metaphlan_filter)),]
metaphlan_filter <- metaphlan_filter[,colSums(metaphlan_filter) != 0]
metaphlan_meta <- data.frame(ID = colnames(metaphlan_filter)) %>%
  left_join(metadata, by = "ID")
metaphlan_summary <- metaphlan_meta %>%
  group_by(Location, Health, sample_type) %>%
  summarise(n = n_distinct(ID))

metaphlan_dist = vegdist(t(metaphlan_filter), method = "bray")
set.seed(1)
metaphlan_tsne <- tsne(metaphlan_dist)
metaphlan_tsne <- as.data.frame(metaphlan_tsne)
rownames(metaphlan_tsne) = names(metaphlan)
metaphlan_tsne$ID <- rownames(metaphlan_tsne)
names(metaphlan_tsne) = c("tsne1", "tsne2", "ID")
metaphlan_tsne <- left_join(metaphlan_tsne, metadata, by = "ID")

# Plot microbiome t-sne
tiff("figures/metaphlan_tsne.tiff", height= 500, width = 800, res = 150)
ggplot(metaphlan_tsne, aes(x = tsne1, y = tsne2, fill = sample_type)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7, pch = 21) +
  scale_fill_manual("Body Site", values = cols,
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("Dim 1") + ylab("Dim 2")
dev.off()

genera <- unique(contig_data$crispr_host[!is.na(contig_data$crispr_host)])
genera <- genera[!genera %in% c(NA, "Lachnospiraceae", "Bacteroidetes")]
metaphlan_genera <- gsub("\\|.*", "", gsub(".*\\|g__", "", row.names(metaphlan_filter)))
match_ind <- c()
genera_split <- rep(NA, length(metaphlan_genera))
for (i in 1:length(genera)) {
  tmp_match_ind <- grep(genera[i], metaphlan_genera)
  match_ind <- c(match_ind, tmp_match_ind)
  genera_split[tmp_match_ind] <- genera[i]
}
genera_split <- genera_split[!is.na(genera_split)]
metaphlan_host <- as.matrix(metaphlan_filter[match_ind,])
metaphlan_host_log <- log10(metaphlan_host)
metaphlan_host_log[is.infinite(metaphlan_host_log)] <- -6

ha = HeatmapAnnotation(type = metaphlan_meta$sample_type,
                       col = list(type = cols),
                       annotation_legend_param = list(type = list(title = "Body Site")),
                       show_annotation_name = FALSE)

# metaphlan_genera <- gsub("\\|.*", "", gsub(".*\\|g__", "", row.names(metaphlan_host)))
tiff("figures/heatplot_metaphlan.tiff", width = 2000, height = 3000, res = 300)
Heatmap(metaphlan_host_log, na_col = "white", top_annotation = ha, name = "log10(rel. abundance)", show_row_names = FALSE, cluster_rows = FALSE,
        col = colorRamp2(c(min(metaphlan_host_log, na.rm = TRUE), max(metaphlan_host_log, na.rm = TRUE)),  c("#f2f2f2", "red4")),
        split = genera_split, row_title_rot = 0, row_title_gp = gpar(fontsize = 5), show_column_names = FALSE,
        heatmap_legend_param = list(color_bar = "continuous"))
dev.off()

# Procrustes analysis
protest_res <- protest(clusters_samples_tsne[,c(1:2)], metaphlan_tsne[,c(1:2)], scale = TRUE)

tiff("figures/procrustes.tiff")
plot(protest_res)
points(protest_res, display = "target", col = "red")
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
  ttest_groups <- df_paired_richness[!duplicated(paste0(df_paired_richness$Location, 
                                                        df_paired_richness$Health,
                                                        df_paired_richness$group)),]
  for(i in 1:nrow(ttest_groups)){
    y <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & 
                              df_paired_richness$Health == ttest_groups$Health[i] &
                              df_paired_richness$group == ttest_groups$group[i],]
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
  ttest_groups <- ttest_groups[order(ttest_groups$Location, ttest_groups$Health, ttest_groups$group),]
  ttest_groups <- ttest_groups[,c("Location", "Health", "group", "pvalue", "asterisk")]
  return(ttest_groups)
}

plotRichnessGraph <- function(df_paired_richness_group, ttest_group, cols) {
  set.seed(1) # for jitter
  g <- ggplot(df_paired_richness_group, aes(sample_type, richness, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.8) +
    theme_classic() +
    ylab("Phage Cluster Richness") +
    xlab("") +
    ggtitle(paste(ttest_group$Location, ttest_group$Health, sep = " - ")) +
    geom_text(data = ttest_group, aes(label=asterisk),
              x = 1.5, y = max(df_paired_richness_group$richness)+10, size = 7,
              inherit.aes = FALSE) +
    scale_fill_manual("Body Site", values = cols[names(cols) %in% unique(df_paired_richness_group$sample_type)]) +
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  return(g)
}

plotMultipleRichnessGraphs <- function(ttest_groups, df_paired, cols){
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  g <- list()
  for(i in 1:nrow(ttest_groups)){
    df_paired_richness_group <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & 
                                                     df_paired_richness$Health == ttest_groups$Health[i] &
                                                     df_paired_richness$group == ttest_groups$group[i],]
    g[[i]] <- plotRichnessGraph(df_paired_richness_group, ttest_groups[i,], cols) +
      theme(legend.position = "none") + ylim(c(0, max(df_paired_richness$richness) + 20))
  }
  return(g)
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Join with metadata
vir_counts_prop_melt_meta <- left_join(metadata, vir_counts_prop_melt, by = c("ID"="Var2","sample_type","Location"))
vir_counts_prop_melt_meta <- vir_counts_prop_melt_meta[vir_counts_prop_melt_meta$Visit_Number == 1,]

# Sample - phage cluster matrix
vir_cluster_counts <- dcast(vir_counts_prop_melt_meta, ID ~ vcontact_cluster, length)
rownames(vir_cluster_counts) <- vir_cluster_counts[,1]
vir_cluster_counts <- vir_cluster_counts[,-1]

# Remove samples with lower than 3 or fewer phage clusters (for subsampling)
remove_ids <- rownames(vir_cluster_counts)[rowSums(vir_cluster_counts > 0) <= 3]
vir_cluster_counts <- vir_cluster_counts[!rownames(vir_cluster_counts) %in% remove_ids,]
metadata_richness <- metadata[!metadata$ID %in% remove_ids,]

pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
paired_metadata <- createPairedData(metadata_richness[metadata_richness$Visit_Number == 1,], pair_list)

# Get richness from matrix
richness_paired <- data.frame(ID = rownames(vir_cluster_counts), richness = rowSums(vir_cluster_counts > 0)) %>%
  right_join(paired_metadata, by = "ID")

richness_ttest <- runTtest(richness_paired)
richness_graphs <- plotMultipleRichnessGraphs(richness_ttest, richness_paired, cols)
richness_graphs[[length(richness_graphs)+1]] <- g_legend(plotRichnessGraph(richness_paired, richness_ttest, cols))

# Plot graph
lay <- rbind(c(1,2,3,13), c(4,5,6,13), c(7,8,9,13), c(10,11,12,13))
tiff("figures/alpha_diversity.tiff", width = 2400, height = 2000, res = 220)
grid.arrange(grobs = richness_graphs, layout_matrix = lay)
dev.off()

# # Subsample matrix and calculate richness
# richness_paired_ss <- data.frame()
# unique_groups <- unique(richness_paired$group)
# for (i in 1:length(unique_groups)) {
# 
#   group_ids <- unique(richness_paired$ID[richness_paired$group %in% unique_groups[i]])
#   vir_cluster_counts_tmp <- vir_cluster_counts[rownames(vir_cluster_counts) %in% group_ids,]
#   min_clusters <- min(rowSums(vir_cluster_counts_tmp))
#   max_clusters <- max(rowSums(vir_cluster_counts_tmp))
# 
#   for(j in 1:(max_clusters - min_clusters)) {
#     vir_cluster_counts_tmp <- t(apply(vir_cluster_counts_tmp, 1, function(x) {
#       if (sum(x) > min_clusters) {
#         ss_index <- sample(1:length(x), 1, prob = ifelse(x > 0, x/sum(x), 0))
#         x[ss_index] <- x[ss_index] - 1
#       }
#       return(x)
#     }))
#   }
# 
#   richness_paired_tmp <- data.frame(ID = rownames(vir_cluster_counts_tmp), richness = rowSums(vir_cluster_counts_tmp > 0)) %>%
#     left_join(metadata_richness, by = "ID") %>%
#     mutate(group = unique_groups[i])
# 
#   richness_paired_ss <- rbind(richness_paired_ss, richness_paired_tmp)
# }
# 
# saveRDS(richness_paired_ss, file = "data/subsampled_phage_cluster_richness.RDS")

# T-test and graphs of subsampled data
richness_paired_ss <- readRDS("data/subsampled_phage_cluster_richness.RDS")
richness_ttest_ss <- runTtest(richness_paired_ss)
richness_graphs_ss <- plotMultipleRichnessGraphs(richness_ttest_ss, richness_paired_ss, cols)
richness_graphs_ss[[length(richness_graphs_ss)+1]] <- g_legend(plotRichnessGraph(richness_paired_ss, richness_ttest_ss, cols))

# Plot graph
lay <- rbind(c(1,2,3,13), c(4,5,6,13), c(7,8,9,13), c(10,11,12,13))
tiff("figures/alpha_diversity_subsampled.tiff", width = 2400, height = 3000, res = 220)
grid.arrange(grobs = richness_graphs_ss, layout_matrix = lay)
dev.off()


########## Size of phages ####
contig_data_meta <- left_join(contig_data, metadata, by = c("sample" = "ID"))
contig_data_meta <- contig_data_meta[!is.na(contig_data_meta$sample_type),]

countSizeCat <- function(contig_data_meta, size_lower, size_upper) {
  size_summary <- contig_data_meta %>%
    group_by(sample_type) %>%
    mutate(n_samples = n_distinct(sample)) %>% ungroup() %>%
    filter(circular, size <= size_upper, size > size_lower) %>%
    group_by(sample_type, n_samples) %>%
    summarise(n = n_distinct(vcontact_cluster)) %>%
    mutate(n_per_sample = n/n_samples) %>%
    mutate(size_lower = size_lower, category = paste(as.character(size_lower/1e3), "<\n<=", as.character(size_upper/1e3)))
  return(size_summary)
}

sizes_lower <- c(0, 1e4, 5e4, 2e5, 3e5)
sizes_upper <- c(sizes_lower[-1], Inf)
size_summary <- map2_df(.x = sizes_lower, .y = sizes_upper, .f = function(.x, .y) countSizeCat(contig_data_meta, .x, .y))
size_summary$category <- gsub("\n<= Inf", "", size_summary$category)

tiff("figures/circular_phage_size.tiff", width = 800, height = 500, res = 150)
ggplot(size_summary, aes(reorder(category, size_lower), n_per_sample, fill = sample_type)) +
  geom_bar(stat = "identity") + xlab("Size (kb)") + ylab("No. unqiue circular phage clusters \n/ no. samples") + 
  theme_classic() +
  scale_fill_manual("Body Site", values = cols)
dev.off()

########## Jumbo/mega circular phages ####
megaphage_contigs_meta <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
megaphage_contigs_meta$Health <- as.character(megaphage_contigs_meta$Health)
megaphage_contigs_meta$Health[megaphage_contigs_meta$Health == "unhealthy"] <- "rheumatoid arthritis"
megaphage_contigs_meta <- megaphage_contigs_meta[megaphage_contigs_meta$Visit_Number == 1,]
megaphage_contigs_meta <- megaphage_contigs_meta[!is.na(megaphage_contigs_meta$circular),]
megaphage_contigs_meta <- megaphage_contigs_meta[megaphage_contigs_meta$circular,]
megaphage_contigs_meta$Location_Health <- paste(megaphage_contigs_meta$Location, "-", megaphage_contigs_meta$Health)
megaphage_contigs_meta$numeric_samplename <- as.numeric(factor(megaphage_contigs_meta$Sample.name))
megaphage_contigs_meta$numeric_samplename[!megaphage_contigs_meta$numeric_samplename %in% unique(megaphage_contigs_meta$numeric_samplename[duplicated(megaphage_contigs_meta$numeric_samplename)])] <- NA
megaphage_contigs_meta$numeric_samplename <- as.numeric(factor(megaphage_contigs_meta$numeric_samplename))
megaphage_contigs_meta$numeric_samplename[is.na(megaphage_contigs_meta$numeric_samplename)] <- c((max(megaphage_contigs_meta$numeric_samplename, na.rm = TRUE)+1):(max(megaphage_contigs_meta$numeric_samplename, na.rm = TRUE)+sum(is.na(megaphage_contigs_meta$numeric_samplename))))

# Percentage of cohorts containing a jumbophage
megaphage_summary <- megaphage_contigs_meta %>%
  group_by(Location, sample_type, Health) %>%
  summarise(n_samples_megaphage = n_distinct(ID))
metadata_summary <- metadata %>%
  group_by(Location, sample_type, Health) %>%
  summarise(n_samples = n_distinct(ID))
megaphage_summary <- full_join(megaphage_summary, metadata_summary, by = c("Location", "sample_type","Health")) %>%
  mutate(perc_samples = 100*n_samples_megaphage/n_samples, Location_Health_sampletype = paste(sample_type, Location, Health, sep = "\n"))

tiff("figures/jumbophage_percentage_samples.tiff", width = 1000, height = 500, res = 100)
ggplot(megaphage_summary, aes(Location_Health_sampletype, perc_samples)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Cohort") + ylab("% samples") + ylim(c(0,100))
dev.off()

label_colours <- c("blue", "orange", "turquoise", "darkgreen")
pd <- position_dodge(0.8)
h <- 200
tiff("figures/jumbophage_individual.tiff", width = 2500, height = 1000, res = 200)
ggplot(megaphage_contigs_meta, aes(sample_type, size/1000, colour = Location_Health, group = numeric_samplename)) +
  geom_line(position=pd, colour="grey90", linetype = "dashed") +
  geom_point(position=pd, alpha = 0.6) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Country - Health", values = label_colours) 
dev.off()

tiff("figures/jumbophage_cluster.tiff", width = 2500, height = 1000, res = 200)
ggplot(megaphage_contigs_meta, aes(sample_type, size/1000, colour = Location_Health, group = vcontact_cluster)) +
  geom_line(position=pd, colour="grey90", linetype = "dashed") +
  geom_point(position = pd, alpha = 0.6) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Country - Health", values = label_colours) 
dev.off() 

tiff("figures/circular_jumbophages.tiff", width = 2500, height = 1000, res = 250)
ggplot(megaphage_contigs_meta, aes(sample_type, size/1000, colour = Location_Health)) +
  geom_jitter(alpha = 0.6, size = 2) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Country - Health", values = label_colours) 
dev.off() 

########## Megaphages and CRISPR hosts ####
megaphage_contigs_meta <- left_join(metadata, megaphage_contigs, by = c("ID"= "sample", "Location"="country"))
megaphage_contigs_meta$Health <- as.character(megaphage_contigs_meta$Health)
megaphage_contigs_meta$Health[megaphage_contigs_meta$Health == "unhealthy"] <- "rheumatoid arthritis"
megaphage_contigs_meta <- megaphage_contigs_meta[megaphage_contigs_meta$Visit_Number == 1,]
megaphage_contigs_meta <- megaphage_contigs_meta[!is.na(megaphage_contigs_meta$circular),]
megaphage_contigs_meta <- megaphage_contigs_meta[megaphage_contigs_meta$circular,]

megaphages_crispr_summary <- megaphage_contigs_meta %>%
  group_by(sample_type, Location, Health, genus) %>%
  summarise(no_crispr_hosts = n())

megaphages_crispr_summary$genus[is.na(megaphages_crispr_summary$genus)] <- "unclassified"

crispr_names <- unique(megaphages_crispr_summary$genus[megaphages_crispr_summary$genus != "unclassified"])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
crispr_colours = c(rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
crispr_colours <- crispr_colours[c(1:length(crispr_names))]
crispr_colours <- c(crispr_colours, "white")
names(crispr_colours) <- c(crispr_names, "unclassified")

tiff("figures/megaphages_crispr_host.tiff", width = 2200, height = 2000, res = 200)
ggplot(megaphages_crispr_summary, aes(sample_type, no_crispr_hosts, fill = genus)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  facet_grid(~ Location + Health, space = "free", scale = "free") +
  theme_classic() +
  ylab("Number of megaphages") + xlab("Body Site") +
  scale_fill_manual("Predicted host", values = crispr_colours)
dev.off()


########## Functional Annotation ####

#### DIAMOND on NCBI NR
megaphages_ncbi <- read.table(file="data/megaphage_proteins_diamond.out", sep = "\t", stringsAsFactors = FALSE)
names(megaphages_ncbi) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore") 

# Filter by best e-value
megaphages_ncbi <- megaphages_ncbi %>% group_by(qseqid) %>%
  filter(evalue == min(evalue))

# Read uniprotkb
uniprot_metadata <- read.table(file="db/UNIPROTKB/uniprotkb_refseq.tab", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Remove duplicates
uniprot_metadata_char <- uniprot_metadata[uniprot_metadata$Protein.names != "Uncharacterized protein",]
uniprot_metadata_char <- uniprot_metadata_char[!duplicated(uniprot_metadata_char$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X),]
uniprot_metadata_unchar <- uniprot_metadata[uniprot_metadata$Protein.names == "Uncharacterized protein",]
uniprot_metadata_unchar <- uniprot_metadata_unchar[!uniprot_metadata_unchar$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X %in% uniprot_metadata_char$yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X,]
uniprot_metadata <- rbind(uniprot_metadata_char, uniprot_metadata_unchar)

megaphages_ncbi <- left_join(megaphages_ncbi, uniprot_metadata, by = c("sseqid" = "yourlist.M201909045C475328CEF75220C360D524E9D456CE3FF4D7X"))
megaphages_ncbi <- megaphages_ncbi[!is.na(megaphages_ncbi$Entry),]
megaphages_ncbi <- megaphages_ncbi %>% select(qseqid, sseqid, Protein.names)

# Remove duplicates
megaphages_ncbi <- megaphages_ncbi[!duplicated(megaphages_ncbi$qseqid),]

#### HMMER on pVOGs
filterBySmallestEvalue <- function(alignment_results) {
  alignment_results_filtered <- alignment_results %>% 
    group_by(query_name) %>%
    filter(evalue == min(evalue)) %>%
    ungroup()
  return(alignment_results_filtered)
}

megaphages_pvogs <- read.table("data/megaphage_pvogs_hmm_tblout.txt", stringsAsFactors = FALSE)

names(megaphages_pvogs) <- c("target_name", "target_accession", "query_name", "query_accession", "evalue", 
                            "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                            "clu", "ov", "env", "dom", "rep", "inc", "description")
megaphages_pvogs <- filterBySmallestEvalue(megaphages_pvogs)

# Join descriptions
pvogs_metadata <- read.table("db/AllvogHMMprofiles/VOGTable_singleprot.txt", stringsAsFactors = FALSE, sep = "\t", quote = "")
megaphages_pvogs <- left_join(megaphages_pvogs, pvogs_metadata, by = c("target_name"="V1")) %>%
  select(query_name, target_name, description = V7)

#### HMMER on Pfam
megaphages_pfam <- read.table(file="data/megaphage_proteins_pfam.out", sep = "", stringsAsFactors = F, header = FALSE, skip = 3)
names(megaphages_pfam) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                            "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                            "clu", "ov", "env", "dom", "rep", "inc")
megaphages_pfam <- filterBySmallestEvalue(megaphages_pfam)

# Join descriptions
pfam_metadata <- data.frame(ids = unique(megaphages_pfam$target_accession))
description <- rep(NA, nrow(pfam_metadata))
for(i in 1:nrow(pfam_metadata)) {
  tmp <- system(paste0("grep -A2 ", pfam_metadata$ids[i], "$ db/PFAM/Pfam-A.hmm.dat"), intern = TRUE)
  description[i] <- gsub("#=GF DE   ", "", tmp[grep("DE   ", tmp)])
}
pfam_metadata$description <- description
megaphages_pfam <- left_join(megaphages_pfam, pfam_metadata, by = c("target_accession"="ids")) %>%
  select(query_name, target_accession, description)

#### HMMER on TIGRFAM (with no eggnog annotation)
megaphages_tigrfams <- read.table(file="data/megaphage_proteins_tigrfams.out", sep = "", stringsAsFactors = FALSE, header = FALSE)
names(megaphages_tigrfams) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                                "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                                "clu", "ov", "env", "dom", "rep", "inc")
megaphages_tigrfams <- filterBySmallestEvalue(megaphages_tigrfams)

# Join descriptions
tigrfams_metadata <- data.frame(ids = unique(megaphages_tigrfams$target_name))
description <- rep(NA, nrow(tigrfams_metadata))
for(i in 1:nrow(tigrfams_metadata)) {
  tmp <- read.table(file = paste0("db/TIGRFAMS/", tigrfams_metadata$ids[i], ".INFO"), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  description[i] <- gsub("DE  ", "", tmp$V1[grep("DE  ", tmp$V1)])
}
tigrfams_metadata$description <- description
megaphages_tigrfams <- left_join(megaphages_tigrfams, tigrfams_metadata, by = c("target_name"="ids")) %>%
  select(query_name, target_name, description)

#### HHBLIT on largest phages
hhblit_result <- data.frame()
for (i in 1:length(largest_phages)) {
  hhblit_result_tmp <- read.table(paste0("data/proteins/", largest_phages[i], ".1E-5.tab"), stringsAsFactors = FALSE)
  hhblit_result <- rbind(hhblit_result, hhblit_result_tmp)
}
names(hhblit_result) <- c("query", "target", "match", "tlen", "mismatch", "gap_open", "qstart", "qend", "tstart", "tend", "eval", "score")
hhblit_result_filter <- hhblit_result %>% group_by(query) %>%
  filter(eval == min(eval))

# uni_lookup <- fread("db/UNICLUST30/uniclust_uniprot_mapping.tsv", stringsAsFactors = FALSE)
# uniprot_id <- rep(NA, nrow(hhblit_result_filter))
# for (i in 1:length(uniprot_id)) {
#   uniprot_id[i] <- uni_lookup$V2[hhblit_result_filter$target[i] + 1]
# }
# write.table(uniprot_id, "data/largest_phages_uniprot_ids.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

uniprot_id <- read.delim("data/largest_phages_uniprot_ids.txt", stringsAsFactors = FALSE, header = FALSE)
uniprot_lookup <- read.delim("db/UNIPROTKB/uniprot_uniclust_largest_phages_lookup.tab", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
hhblit_result_filter$Entry <- uniprot_id$V1
megaphages_hhblit <- left_join(hhblit_result_filter, uniprot_lookup, by = "Entry") %>%
  select(query, Entry, Protein.names)

#### Combine results (in order ncbi, pvogs, pfam, tigrfam)
header_names <- c("coding_region", "protein", "protein_description")
names(megaphages_ncbi) <- header_names
names(megaphages_pvogs) <- header_names
names(megaphages_pfam) <- header_names
names(megaphages_tigrfams) <- header_names
names(megaphages_hhblit) <- header_names
megaphages_proteins <- rbind(as.data.frame(megaphages_ncbi[megaphages_ncbi$protein_description != "Uncharacterized protein",]),
                             as.data.frame(megaphages_pvogs), as.data.frame(megaphages_pfam), as.data.frame(megaphages_tigrfams),
                             as.data.frame(megaphages_hhblit))

# Check duplications sum(duplicated(megaphages_proteins$contig))
megaphages_proteins <- megaphages_proteins[!duplicated(megaphages_proteins$coding_region),]

# Extract contig id
megaphages_proteins$contig <- sub("_[^_]+$", "", megaphages_proteins$coding_region)

# Modify prodigal GFF to include protein descriptions
prodigal_gff <- read.table("data/megaphage_contigs_prodigal.gff", stringsAsFactors = FALSE)
prodigal_gff$coding_region <- paste0(prodigal_gff$V1, "_", gsub(".*_", "", gsub(";.*", "", prodigal_gff$V9)))
prodigal_gff <- left_join(prodigal_gff, megaphages_proteins, by = "coding_region")

# Read protein family lookup
prot_lookup <- read.table("db/ONTOLOGIES/protein_family_lookup.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
prot_lookup <- prot_lookup[order(prot_lookup$description_contains),] # Order so more unique lookups get selected

prodigal_gff$Name <- rep(NA, nrow(prodigal_gff))
prodigal_gff$Parent <- rep(NA, nrow(prodigal_gff))
for (i in 1:nrow(prot_lookup)) {
  tmp_inds <- grep(prot_lookup$description_contains[i], prodigal_gff$protein_description)
  prodigal_gff$Name[tmp_inds] <- prot_lookup$name[i]
  prodigal_gff$Parent[tmp_inds] <- prot_lookup$parent[i]
}

# Check all jumbo phages have structural proteins
struct_prots <- prodigal_gff %>% filter(Parent == "Structural proteins")

struct_prot_count <- prodigal_gff %>% group_by(V1, Parent) %>%
  summarise(n = n()) %>%
  filter(Parent == "Structural proteins") # All jumbophages have structural proteins making them putative phages

struct_prot_list <- prodigal_gff %>%
  filter(Parent == "Structural proteins")

# Add Name and Parent to attributes
prodigal_gff$V9 <- ifelse(is.na(prodigal_gff$Name), 
                          paste0(prodigal_gff$V9, "Name=hypothetical protein;"),
                          paste0(prodigal_gff$V9, "Name=", gsub(" \\[.*", "", prodigal_gff$Name), ";"))
prodigal_gff$V9 <- ifelse(is.na(prodigal_gff$Parent), 
                          paste0(prodigal_gff$V9, "Parent=None;"),
                          paste0(prodigal_gff$V9, "Parent=", prodigal_gff$Parent, ";"))
prodigal_gff <- prodigal_gff[,c(1:9)]
names(prodigal_gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Add tRNA
megaphages_trna <- read.table("data/megaphage_trna_clean.tab", stringsAsFactors = FALSE)
megaphages_trna$strand <- ifelse(grepl("c", megaphages_trna$V3), "-", "+")
megaphages_trna <- cbind(megaphages_trna, map_df(.x = gsub("c", "", gsub("\\]", "", gsub("\\[", "", megaphages_trna$V3))), 
                                                 function(.x) {
                                                   split <- strsplit(.x, ",")
                                                   return(data.frame(start = as.numeric(split[[1]][1]), 
                                                                     end = as.numeric(split[[1]][2])))
                                                   }))
megaphages_trna$attributes <- paste0("Name=", megaphages_trna$V2, ";")
trna_gff <- megaphages_trna %>% mutate(source = "ARAGORN_v1.2.36", type = "tRNA", phase = 0, score = format(as.numeric(V4), nsmall=1)) %>%
  select(seqid = V6, source, type, start, end, score, strand, phase, attributes)

# Combine prodigal CDS and trna
comb_gff <- rbind(prodigal_gff, trna_gff)

# Write GFF file for largest phages
largest_phages <- megaphage_contigs_meta$name[order(megaphage_contigs_meta$size, decreasing = TRUE)][c(1:3)]
for (i in 1:length(largest_phages)) {
  largest_phage_gff <- comb_gff[comb_gff$seqid == largest_phages[i],]
  write.table(largest_phage_gff, paste0("data/", largest_phages[i], ".gff"), 
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")                          
}

######## EXTRA #########
########## Phages and ARGs 
# Read data
contig_data_meta <- left_join(metadata, contig_data, by = c("ID"= "sample", "Location"="country"))
contig_data_meta <- contig_data_meta[contig_data_meta$Visit_Number == 1,]
contig_data_meta$Health <- as.character(contig_data_meta$Health)
contig_data_meta <- contig_data_meta[!is.na(contig_data_meta$circular),]
contig_data_meta$Health[contig_data_meta$Health == "unhealthy"] <- "rheumatoid arthritis"
contig_data_meta$Location_Health <- paste(contig_data_meta$Location, "-", contig_data_meta$Health)
contig_data_meta <- contig_data_meta %>% group_by(Location, sample_type, Health) %>% mutate(num_samples = n_distinct(ID))

args_data <- read.csv("data/all_assemblies_card_90.csv", stringsAsFactors = FALSE)
args_data$ID <- gsub("_card.out", "", gsub(".*/", "", args_data$filename))
args_data$name <- paste0(args_data$ID, "_", args_data$qseqid)

# Combine phage and arg results where only phage does/doesnt contain ARG
phages_args_yn <- left_join(contig_data_meta, args_data, by = c("name", "ID"), suffix = c(".phage", ".arg"))
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
phages_args <- inner_join(contig_data_meta, args_data, by = c("name", "ID"), suffix = c(".phage", ".arg"))
nrow(phages_args)/nrow(contig_data_meta)*100
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
  distinct() %>% left_join(contig_data_meta, by = c("Location", "sample_type", "Health", "num_samples", "vcontact_cluster")) %>%
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


########## Auxilliary Metabolic Genes 
#### EGGNOG
megaphage_prots <- read.delim("data/megaphage_proteins_headers.txt", stringsAsFactors = FALSE, header = FALSE)
names(megaphage_prots) <- "query_name"
eggnog <- read.delim("data/megaphage_annot.emapper.annotations", stringsAsFactors = FALSE, header = FALSE)
names(eggnog) <- c("query_name", "seed_eggnog_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "predicted_taxonomic_group",
                   "predicted_protein_name", "gene_ontology_terms", "EC_number", "KEGG_ko", "KEGG_pathway", "KEGG_module", "KEGG_reaction",
                   "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_reaction", "tax_scope", "eggnog_ogs", "bestOG", "COG_functional_cat", "eggnog_description")

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

# Summarise functional groups for each cohort
megaphage_eggnog_annot <- megaphage_eggnog[!is.na(megaphage_eggnog$seed_eggnog_ortholog),]
megaphage_eggnog_annot <- megaphage_eggnog_annot[megaphage_eggnog_annot$COG_functional_cat != "",]
megaphage_eggnog_cogs <- megaphage_eggnog_annot %>%
  mutate(COG_functional_cat = strsplit(as.character(COG_functional_cat), "")) %>%
  unnest(COG_functional_cat)
functional_cats <- c("Translation, ribosomal structure and biogenesis", "RNA processing and modification", 
                     "Transcription", "Replication, recombination and repair", "Chromatin structure and dynamics",
                     "Cell cycle control, cell division, chromosome partitioning", "Nuclear structure", "Defense mechanisms",
                     "Signal transduction mechanisms", "Cell wall/membrane/envelope biogenesis", "Cell motility",
                     "Cytoskeleton", "Extracellular structures", "Intracellular trafficking, secretion, and vesicular transport",
                     "Posttranslational modification, protein turnover, chaperones", "Energy production and conversion", 
                     "Carbohydrate transport and metabolism", "Amino acid transport and metabolism", "Nucleotide transport and metabolism", 
                     "Coenzyme transport and metabolism", "Lipid transport and metabolism", "Inorganic ion transport and metabolism", 
                     "Secondary metabolites biosynthesis, transport and catabolism", "General function prediction only", "Function unknown")
names(functional_cats) <- c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S")
megaphage_eggnog_cogs$functional_cats <- sapply(megaphage_eggnog_cogs$COG_functional_cat, function(x) functional_cats[names(functional_cats) == x])
metabolic_summary <- megaphage_eggnog_cogs %>% group_by(Location, Health, sample_type, functional_cats) %>%
  summarise(n_cog = n()) %>%
  group_by(Location, Health, sample_type) %>%
  mutate(sum_n_cog = sum(n_cog)) %>%
  mutate(per_cog = n_cog/sum_n_cog * 100)

# Plot functional groups for each body site and cohort
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
functional_colours <- c(rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
names(functional_colours) <- functional_cats
functional_colours <- functional_colours[!is.na(names(functional_colours))]
functional_colours[names(functional_colours) == "Function unknown"] <- "white"

tiff("figures/functional_categories_bodysite.tiff", width = 2500, height = 1000, res = 150)
ggplot(metabolic_summary, aes(sample_type, per_cog, fill = functional_cats)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  facet_grid(~ Location + Health, space = "free", scale = "free") +
  scale_fill_manual("Functional Categories", values = functional_colours) +
  xlab("Body Site") + ylab("Percentage") + theme_classic()
dev.off()

# Summarise functional groups for each cluster
metabolic_cluster_summary <- megaphage_eggnog_cogs %>% group_by(genus, functional_cats) %>%
  summarise(n_cog = n()) %>%
  group_by(genus) %>%
  mutate(sum_n_cog = sum(n_cog)) %>%
  mutate(per_cog = n_cog/sum_n_cog * 100)

# Plot functional groups for each cluster
metabolic_cluster_summary$genus <- factor(metabolic_cluster_summary$genus, levels=c("Neisseria", "Veillonella", "Centipeda", "Lachnoanaerobaculum", "Streptococcus", "Rothia", "Atopobium", "Selenomonas", NA))
tiff("figures/functional_categories_host.tiff", width = 2600, height = 1000, res = 150)
ggplot(metabolic_cluster_summary, aes(genus, per_cog, fill = functional_cats)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  scale_fill_manual("Functional Categories", values = functional_colours) +
  xlab("Predicted host") + ylab("Percentage") + theme_classic()
dev.off()

# Group into metabolic profiles
megaphage_cogs <- acast(megaphage_eggnog_cogs, name ~ seed_eggnog_ortholog)
cogs_dist <- dist(megaphage_cogs, "binary")
set.seed(1)
mds <- wcmdscale(cogs_dist, w=rep(1,nrow(megaphage_cogs)))
mds <- as.data.frame(mds[,1:2])
mds$name <- row.names(mds)
mds <- left_join(mds, megaphage_contigs, by = "name") %>%
  left_join(metadata, by = c("sample"="ID"))

# Plot contig size and crispr host
tiff("figures/pcoa_eggnog_crisprhost_size.tiff", width = 1000, height = 750, res = 150)
ggplot(mds, aes(V1, V2, size = size, fill = genus)) +
  theme_classic() + geom_point(pch = 21, alpha = 0.7) +
  scale_fill_brewer("Predicted host", palette = "Set2",
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("PCo 1") + ylab("PCo 2")
dev.off()

# Plot body site and Country - Health
tiff("figures/pcoa_eggnog_bodysite_cohort.tiff", width = 1000, height = 750, res = 150)
ggplot(mds, aes(V1, V2, fill = sample_type, shape = Location_Health)) +
  geom_point(size = 3, alpha = 0.7) + theme_classic() +
  scale_shape_manual("Country - Health", values = c(21, 24, 22, 23)) +
  scale_fill_brewer("Body Site", palette = "Set3",
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("PCo 1") + ylab("PCo 2") 
dev.off()

# Plot age and gender
tiff("figures/pcoa_eggnog_age_gender.tiff", width = 1000, height = 750, res = 150)
ggplot(mds, aes(V1, V2, fill = Age, shape = Gender)) +
  theme_classic() + geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual("Sex", values = c(21, 22)) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("PCo 1") + ylab("PCo 2")
dev.off()
