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
library(viridis)
library(ggalt)
library(VennDiagram)
set.seed(1)

########## Read data, metadata and define parameters ####
# Data
vir_counts_prop_melt <- readRDS("data/vir_counts_prop_melt.RDS")
vir_counts_prop_agg <- readRDS("data/vir_counts_prop_agg.RDS")
vir_counts_prop_agg_phage <- readRDS("data/vir_counts_prop_agg_phage.RDS")
vir_counts_prop_melt_agg <- readRDS("data/vir_counts_prop_melt_agg.RDS")
vir_counts_prop_melt_agg2 <- readRDS("data/vir_counts_prop_melt_agg2.RDS")
contig_data <- readRDS("data/contig_data.RDS")
counts_total <- readRDS("data/counts_total.RDS")
jumbophage_contigs <- read.delim("data/jumbophage_contigs.txt", sep = " ", stringsAsFactors = FALSE)
metaphlan <- read.csv("data/all_metaphlan.csv", stringsAsFactors = FALSE)

# Metadata
metadata = read.csv("data/metadata_v2.csv", stringsAsFactors = FALSE)
rownames(metadata) = metadata$ID
metadata <- metadata[metadata$ID %in% unique(vir_counts_prop_melt$Var2),]
metadata$Location_sampletype <- paste(metadata$Location, "-", metadata$sample_type)

# Create Supplementary Table of metadata
metadata$n_contigs <- sapply(metadata$ID, function(x) length(which(vir_counts_prop_melt$Var2 == x)))
metadata$total_reads <- as.numeric(counts_total[metadata$ID])
write.csv(metadata, "data/Supplementary_Table1.csv", row.names = FALSE)

# n summary of metadata excluding longitudinal samples
metadata_summary <- metadata %>% filter(Visit_Number == 1) %>% group_by(Location, sample_type) %>% summarise(n())

# Body site colours
cols <- plasma(length(unique(metadata$sample_type)), end = 0.8)
names(cols) <- sort(unique(metadata$sample_type))

# Cohort colours
metadata$Location_sampletype <- paste(metadata$sample_type, "-", metadata$Location)
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "gold", brewer.pal(9, "YlOrRd")[c(5,7)], brewer.pal(9, "RdPu")[c(3,5)])
names(cohort_cols) <- sort(unique(metadata$Location_sampletype)) 

########## Virome composition barplots ####
demovir_cols <- rev(c(brewer.pal(length(unique(vir_counts_prop_melt_agg$demovir)), "Set3")))
names(demovir_cols) <- sort(as.character(unique(vir_counts_prop_melt_agg$demovir)))
demovir_cols[names(demovir_cols) == "crAss-like"] <- "grey50"
demovir_cols[names(demovir_cols) == "Unassigned"] <- "grey80"
demovir_cols[demovir_cols == "#FFFFB3"] <- "#8DD3C7"

vir_counts_prop_melt_agg_exus <- vir_counts_prop_melt_agg %>%
  filter(Var2 %in% metadata$ID[metadata$Visit_Number == 1])

tiff("figures/barplot_demovir.tiff", width = 5000, height = 2000, res = 400)
ggplot(vir_counts_prop_melt_agg_exus, aes(x = Var2, y = V1, fill = as.factor(demovir))) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols[sort(names(demovir_cols))]) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of reads") + xlab("Sample") + ylim(0, max(aggregate(vir_counts_prop_melt_agg_exus$V1, by=list(ID=vir_counts_prop_melt_agg_exus$Var2), FUN=sum)$x))
dev.off()

vir_counts_prop_melt_agg_rel <- vir_counts_prop_melt_agg_exus %>%
  group_by(Var2) %>%
  mutate(V1_prop = V1/sum(V1))
  
tiff("figures/barplot_demovir_rel.tiff", width = 5000, height = 2000, res = 400)
ggplot(vir_counts_prop_melt_agg_rel, aes(x = Var2, y = V1_prop, fill = as.factor(demovir))) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols[sort(names(demovir_cols))]) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of reads") + xlab("Sample")
dev.off()

plotFamilyHeatmap <- function(family_matrix, title, cols = cohort_cols) {
  heatmap.2(family_matrix,
            trace = "none",
            scale = "none",
            density.info = "none",
            hclustfun = function(x) {hclust(x, method = "ward.D2")},
            dendrogram = "both",
            col =  c("white", brewer.pal(9, "PuRd")[3:9]),
            breaks = seq(min(family_matrix), 0, length.out = 9),
            symbreaks = FALSE,
            keysize = 1,
            lhei = c(1,8),
            key.title = NA,
            key.xlab = "log10(prop. reads mapped)",
            key.ylab = NA,
            ColSideColors = cols[metadata[colnames(family_matrix), "Location_sampletype"]],
            labRow = NA,
            labCol = NA,
            cexCol = 0.5,
            cexRow = 0.5,
            xlab = "",
            ylab = "Phage Clusters",
            main = title)
}

unique_demovir <- sort(unique(vir_counts_prop_melt_agg$demovir))
unique_demovir <- as.character(unique_demovir[!unique_demovir %in% c("Other", "Unassigned")])
num_clusters <- rep(NA, length(unique_demovir))
list_family_mat <- list()

for (i in 1:length(unique_demovir)) {
  family_prop <- vir_counts_prop_melt_agg %>%
    filter(Var2 %in% metadata$ID[metadata$Visit_Number == 1]) %>%
    filter(demovir == unique_demovir[i]) %>%
    dcast(vcontact_cluster ~ Var2, sum)
  
  row.names(family_prop) <- family_prop$vcontact_cluster
  family_prop <- as.matrix(family_prop[,-1])
  num_clusters[i] <- dim(family_prop)[1]
  family_prop_log = log10(as.matrix(family_prop))
  family_prop_log[is.infinite(family_prop_log)] = -8
  list_family_mat[[i]] <- family_prop_log
}
names(list_family_mat) <- unique_demovir
names(num_clusters) <- unique_demovir

tiff("figures/Inoviridae.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["Inoviridae"]], "Inoviridae")
dev.off()

tiff("figures/Microviridae.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["Microviridae"]], "Microviridae")
dev.off()

tiff("figures/Siphoviridae.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["Siphoviridae"]], "Siphoviridae")
dev.off()

tiff("figures/Myoviridae.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["Myoviridae"]], "Myoviridae")
dev.off()

tiff("figures/Podoviridae.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["Podoviridae"]], "Podoviridae")
dev.off()

tiff("figures/crAss-like.tiff", width = 2000, height = 1400, res = 150)
plotFamilyHeatmap(list_family_mat[["crAss-like"]], "crAss-like")
dev.off()

# Bicaudaviridae
metadata[colnames(list_family_mat[["Bicaudaviridae"]]), "Location_sampletype"]

tiff("figures/cohort_legend.tiff")
plot.new()
legend("bottomleft", legend = levels(factor(metadata[, "Location_sampletype"])),
       col = cohort_cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 1, title = "Cohort:")
dev.off()

########## Virome beta-diversity###########################
# Calculate number of samples with number of unique phages and vcontact clusters
samples_with_n_phage <- vir_counts_prop_melt %>% 
  group_by(Var1) %>%
  summarise(n = n_distinct(Var2)) %>%
  group_by(n) %>%
  summarise(num_samples = n()) %>%
  mutate(perc_samples = num_samples/sum(num_samples)*100)

samples_with_n_clusters <- vir_counts_prop_melt %>% 
  filter(vcontact_cluster != "" | !is.na(vcontact_cluster)) %>%
  group_by(vcontact_cluster) %>%
  summarise(n = n_distinct(Var2)) %>%
  group_by(n) %>%
  summarise(num_samples = n()) %>%
  mutate(perc_samples = num_samples/sum(num_samples)*100)

tiff("figures/unique_phages.tiff", width= 500, height = 500, res = 100)
ggplot(samples_with_n_phage, aes(n, perc_samples)) +
  geom_bar(stat="identity") + theme_classic() +
  xlab("No. unique phages") + ylab("% samples") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  ylim(c(0, c(max(samples_with_n_phage$perc_samples, samples_with_n_clusters$perc_samples))))
dev.off()

tiff("figures/unique_clusters.tiff", width= 500, height = 500, res = 100)
ggplot(samples_with_n_clusters, aes(n, perc_samples)) +
  geom_bar(stat="identity") + theme_classic() +
  xlab("No. unique phage clusters") + ylab("% samples") + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  ylim(c(0, c(max(samples_with_n_phage$perc_samples, samples_with_n_clusters$perc_samples))))
dev.off()

# Compute distance matrix and run t-SNE on phage clusters
vir_counts_prop_agg_meta <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% metadata$ID[metadata$Visit_Number == 1]]
cluster_counts_dist = vegdist(t(vir_counts_prop_agg_meta), method = "bray")
clusters_samples_tsne = tsne(cluster_counts_dist)

#Annotate t-SNE coordinates
clusters_samples_tsne = as.data.frame(clusters_samples_tsne)
rownames(clusters_samples_tsne) = colnames(vir_counts_prop_agg_meta)
clusters_samples_tsne$ID <- rownames(clusters_samples_tsne)
colnames(clusters_samples_tsne) = c("tsne1", "tsne2", "ID")
clusters_samples_tsne <- left_join(clusters_samples_tsne, metadata, by = "ID")
rownames(clusters_samples_tsne) <- clusters_samples_tsne$ID

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 2:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(clusters_samples_tsne[,c("tsne1", "tsne2")], k = k), clusters_samples_tsne[,c("tsne1", "tsne2")])
  avg_sil[k-1] <- mean(tmp[,3])
}
# Group by silhouette width
samples_clust <- pam(clusters_samples_tsne[,c("tsne1", "tsne2")], which.max(avg_sil)+1)

clusters_samples_tsne$cluster = as.factor(samples_clust$cluster[clusters_samples_tsne$ID])
clusters_samples_tsne$n_contigs <- sapply(clusters_samples_tsne$ID, function(x) length(which(vir_counts_prop_melt$Var2 == x)))
clusters_samples_tsne$total_reads = as.numeric(counts_total[clusters_samples_tsne$ID])

saveRDS(clusters_samples_tsne,  file = "data/clusters_samples_tsne.RDS")

# Make t-SNE plots
group_cols <- viridis(length(unique(clusters_samples_tsne$cluster)))

# Calculate proportion of samples
cluster_res <- clusters_samples_tsne %>% group_by(cluster, Location, sample_type) %>% summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(total_n = sum(n)) %>%
  mutate(prop_cluster = n/total_n*100, Location_sampletype = paste(sample_type, "-", Location))

# Calculate proportion of body sites
cluster_site <- cluster_res %>% 
  group_by(cluster, sample_type) %>%
  summarise(prop_cluster = signif(sum(prop_cluster), 3)) %>%
  mutate(summary = paste(sample_type, ": ", prop_cluster, "%")) %>%
  group_by(cluster) %>%
  summarise(summary = paste(summary, collapse = "; ")) %>%
  mutate(summary = paste0("Group ", cluster, " (", summary, ")"))

# Label by cluster
tsne_plot1 <- ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, colour = sample_type, shape = cluster)) +
  theme_classic() + 
  geom_text(aes(label=cluster)) +
  scale_colour_manual("Body Site", values = cols, 
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

# Label by sex and age
tsne_plot4 = ggplot(clusters_samples_tsne, aes(x = tsne1, y = tsne2, fill = Age, shape = Gender)) +
  theme_classic() + geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual("Sex", values = c(21, 22)) +
  #scale_fill_gradient2()
  scale_fill_distiller(palette = "Spectral") +
  xlab("Dim 1") + ylab("Dim 2")

# Save plots
tiff("figures/tsne_clusters.tiff", width = 800, height = 500, res = 150)
tsne_plot1
dev.off()

tiff("figures/tsne_bodysite_location.tiff", width = 1000, height = 500, res = 150)
ggplot(cluster_res, aes(cluster, prop_cluster, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  theme_classic() + xlab("Group") + ylab("Percentage") +
  scale_fill_manual("Body Site - Geographical Location", values = cohort_cols)
dev.off()

tiff("figures/tsne_bodysite.tiff", width = 800, height = 500, res = 150)
tsne_plot3
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
          xlab = "",
          ylab = "Phage Clusters",
          main = NA
)
legend(x = 0.8, y = 0.35, legend = levels(factor(metadata[colnames(vir_counts_prop_agg_meta_diff), "sample_type"])),
       col = cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.8, title = "Sample clusters:")

legend(x = 0.8, y = 0.2, xpd=TRUE, legend = levels(factor(viral_clusters_df[rownames(viral_clusters_df) %in% rownames(vir_counts_prop_agg_meta_diff),"demovir"], levels = names(demovir_cols))),
       col = demovir_cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.8, title = "Viral Family")

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
suma2Venn(vir_group_all, zcolor = group_cols, cexil = 1.5, cexsn = 1.4)
dev.off()

########## Longitudinal analysis and core phageome #########
# Create mock IDs so different time points
metadata_us <- metadata[metadata$Location == "US",]
metadata_longus <- metadata_us %>% group_by(Sample.name, sample_type) %>%
  filter(any(Visit_Number == 3))

# Remove demovir predictions
vir_counts_prop_melt_agg3 = vir_counts_prop_melt_agg[, sum(V1), by=.(Var2, vcontact_cluster, sample_type, Location)]

vir_counts_longus <- left_join(metadata_longus[,c("ID", "Sample.name", "Visit_Number")], vir_counts_prop_melt_agg3, by = c("ID"="Var2"))

# Count average time points of phages
cluster_no_tp <- vir_counts_longus %>% group_by(sample_type, Sample.name, vcontact_cluster) %>%
  summarise(no_tp = n_distinct(Visit_Number)) 

# Count all distinct clusters
cluster_all_summary <- vir_counts_longus %>% group_by(sample_type, Sample.name) %>%
  summarise(total_cluster = n_distinct(vcontact_cluster))

# Count clusters in at least x number of time points
cluster_count_tp <- map_df(.x = unique(cluster_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(cluster_no_tp$sample_type[cluster_no_tp$no_tp == .x], .x),
                    Sample.name = rep(cluster_no_tp$Sample.name[cluster_no_tp$no_tp == .x], .x),
                    vcontact_cluster = rep(cluster_no_tp$vcontact_cluster[cluster_no_tp$no_tp == .x]),
                    no_tp = rep(cluster_no_tp$no_tp[cluster_no_tp$no_tp == .x], .x))
  tmp$tp <- c(1:.x)
  return(tmp)
})

cluster_count_tp_summary <- cluster_count_tp %>% group_by(sample_type, Sample.name, tp) %>%
  summarise(total_cluster_least = n_distinct(vcontact_cluster))

cluster_one <- cluster_count_tp_summary %>%
  filter(tp == 1) %>%
  rename(total_cluster = total_cluster_least) %>%
  select(-tp)

cluster_count_tp_summary <- left_join(cluster_count_tp_summary, cluster_one, by = c("sample_type", "Sample.name")) %>%
  mutate(cluster_frac = total_cluster_least/total_cluster)

# Plot
tiff("figures/longitudinal_cluster_count.tiff", height = 400, width = 1100, res = 100)
set.seed(1)
ggplot(cluster_count_tp_summary, aes(as.factor(tp), cluster_frac, fill = as.factor(sample_type))) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of phage clusters in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = cols)
dev.off()

# Proportion of reads mapped to clusters
cluster_count_reads<- vir_counts_longus %>% 
  group_by(sample_type, Sample.name, vcontact_cluster) %>%
  arrange(desc(Visit_Number)) %>%
  mutate(cumsum_reads = cumsum(V1)) %>%
  select(sample_type, Sample.name, vcontact_cluster, cumsum_reads, Visit_Number)

cluster_count_reads_tp1 <- cluster_count_reads %>%
  filter(Visit_Number == min(Visit_Number)) %>%
  rename(tp1 = cumsum_reads) %>%
  select(sample_type, Sample.name, vcontact_cluster, tp1)

cluster_count_reads_summary <- left_join(cluster_count_reads, cluster_count_reads_tp1, by = c("sample_type", "vcontact_cluster", "Sample.name")) %>%
  mutate(frac_reads = cumsum_reads/tp1) %>%
  group_by(sample_type, Sample.name, vcontact_cluster) %>%
  mutate(tp = order(Visit_Number))

tiff("figures/longitudinal_read_count.tiff", height = 400, width = 1100, res = 90)
ggplot(cluster_count_reads_summary, aes(x = as.factor(tp), y = frac_reads, fill = sample_type)) +
  geom_boxplot() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. of time points") + ylab("Proportion of reads mapped in at least x time points") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual("body site", values = cols)
dev.off()

vir_counts_longus_demovir <- left_join(metadata_longus[,c("ID", "Sample.name", "Visit_Number")], vir_counts_prop_melt_agg, by = c("ID"="Var2"))

# Phage clusters that are persistent/transient
cluster_all <- left_join(vir_counts_longus_demovir, cluster_no_tp, by = c("sample_type", "Sample.name", "vcontact_cluster"))

# Select persistent clusters (3 at least, 2 doesnt converge)
cluster_persist <- cluster_all %>%
  filter(no_tp >= 3)

# Percentage phage clusters that are found to be persistent
length(unique(cluster_persist$vcontact_cluster))/length(unique(cluster_all$vcontact_cluster)) * 100
length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "buccal mucosa"]))/length(unique(cluster_all$vcontact_cluster[cluster_all$sample_type == "buccal mucosa"])) * 100
length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dorsum of tongue"]))/length(unique(cluster_all$vcontact_cluster[cluster_all$sample_type == "dorsum of tongue"])) * 100
length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dental"]))/length(unique(cluster_all$vcontact_cluster[cluster_all$sample_type == "dental"])) * 100
length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "stool"]))/length(unique(cluster_all$vcontact_cluster[cluster_all$sample_type == "stool"])) * 100

# Relative abundance of mapped reads to persistent phages
cluster_persist_rel <- cluster_persist %>%
  group_by(ID) %>%
  mutate(V1_prop = V1/sum(V1))

tiff("figures/barplot_demovir_persistent_rel.tiff", width = 5000, height = 2000, res = 400)
ggplot(cluster_persist_rel, aes(x = ID, y = V1_prop, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of mapped reads") + xlab("Sample")
dev.off()

# Select nonpersistent clusters
cluster_nonpersist <- cluster_all %>%
  filter(no_tp < 3)

# Percentage of persistent phage clusters that can also be transient
sum(unique(cluster_persist$vcontact_cluster) %in% unique(cluster_nonpersist$vcontact_cluster))/length(unique(cluster_persist$vcontact_cluster))*100
sum(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "buccal mucosa"]) %in% unique(cluster_nonpersist$vcontact_cluster[cluster_nonpersist$sample_type == "buccal mucosa"]))/length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "buccal mucosa"]))*100
sum(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dental"]) %in% unique(cluster_nonpersist$vcontact_cluster[cluster_nonpersist$sample_type == "dental"]))/length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dental"]))*100
sum(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dorsum of tongue"]) %in% unique(cluster_nonpersist$vcontact_cluster[cluster_nonpersist$sample_type == "dorsum of tongue"]))/length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "dorsum of tongue"]))*100
sum(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "stool"]) %in% unique(cluster_nonpersist$vcontact_cluster[cluster_nonpersist$sample_type == "stool"]))/length(unique(cluster_persist$vcontact_cluster[cluster_persist$sample_type == "stool"]))*100

# Relative abundance of mapped reads to transient phages
cluster_nonpersist_rel <- cluster_nonpersist %>%
  group_by(ID) %>%
  mutate(V1_prop = V1/sum(V1)) 

tiff("figures/barplot_demovir_transient_rel.tiff", width = 5000, height = 2000, res = 400)
ggplot(cluster_nonpersist_rel, aes(x = ID, y = V1_prop, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of mapped reads") + xlab("Sample")
dev.off()

# Proportion of phage clusters that are persistent/transient/both
cluster_all$status <- "both"
unique_sample_types <- sort(unique(cluster_all$sample_type))
for (i in 1:length(unique_sample_types)) {
  cluster_all$status[cluster_all$vcontact_cluster %in% unique(cluster_persist$vcontact_cluster) & !(cluster_all$vcontact_cluster %in% unique(cluster_nonpersist$vcontact_cluster)) & cluster_all$sample_type == unique_sample_types[i]] <- "persistent"
  cluster_all$status[!(cluster_all$vcontact_cluster %in% unique(cluster_persist$vcontact_cluster)) & cluster_all$vcontact_cluster %in% unique(cluster_nonpersist$vcontact_cluster) & cluster_all$sample_type == unique_sample_types[i]] <- "transient"
}

cluster_all_perc <- cluster_all %>%
  group_by(sample_type) %>%
  mutate(n_individuals = n_distinct(Sample.name)) %>%
  group_by(sample_type, vcontact_cluster, n_individuals, status) %>%
  summarise(n = n_distinct(Sample.name)) %>%
  mutate(perc = n/n_individuals*100) %>% ungroup()

perc_categories <- seq(0,100,10)
cluster_all_perc$perc_cat <- NA

for (i in 1:(length(perc_categories)-1)) {
  cluster_all_perc$perc_cat[cluster_all_perc$perc > perc_categories[i] & cluster_all_perc$perc <= perc_categories[i+1]] <- paste(perc_categories[i], "<\n<=", perc_categories[i+1])
}

status_all_perc <- cluster_all_perc %>%
  group_by(sample_type) %>%
  mutate(n_clusters = n_distinct(vcontact_cluster)) %>%
  group_by(sample_type, perc_cat, n_clusters, status) %>%
  summarise(n = n()) %>%
  mutate(perc = n/n_clusters*100)

tiff("figures/both_cluster_perc.tiff", width = 6000, height = 2000, res = 300)
ggplot(status_all_perc, aes(perc_cat, perc, fill = factor(status))) +
  geom_bar(stat="identity") +
  facet_grid(~sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("% Individuals") + ylab("% Phage Clusters") + ylim(c(0,100)) +
  scale_fill_manual("", values = brewer.pal(length(unique(status_all_perc$status)), "Set2"), 
                    labels = c("Both (persistent and transient in at least one individual)", "Persistent", "Transient"))
dev.off()

cluster_all_perc_demovir <- left_join(cluster_all_perc, vir_counts_longus_demovir[!duplicated(paste(vir_counts_longus_demovir$sample_type, vir_counts_longus_demovir$vcontact_cluster)),], by = c("sample_type", "vcontact_cluster"))

# Transient percentage level
transient_perc <- cluster_all_perc %>%
  group_by(sample_type, n_individuals) %>%
  filter(status == "transient") %>%
  summarise(max_perc = max(perc)) %>%
  mutate(n = n_individuals * max_perc * 0.01)

# Persistent percentage level
persistent_perc <- cluster_all_perc %>%
  group_by(sample_type, n_individuals) %>%
  filter(status == "persistent") %>%
  summarise(max_perc = max(perc)) %>%
  mutate(n = n_individuals * max_perc * 0.01)

# Core cutoff
core_cutoff <- transient_perc$max_perc
names(core_cutoff) <- unique_sample_types

# Get core of each body site
cluster_core <- data.frame()
for (i in 1:length(unique_sample_types)) {
  cluster_core <- rbind(cluster_core, cluster_all_perc[cluster_all_perc$sample_type == unique_sample_types[i] & cluster_all_perc$perc > core_cutoff[names(core_cutoff) == unique_sample_types[i]],])
}

# Number of phage clusters in each body site
clusters_body_site <- cluster_all_perc %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(vcontact_cluster))

# number of phage clusters in core
cluster_core %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(vcontact_cluster)) %>%
  mutate(n_clusters = clusters_body_site$n) %>%
  mutate(n/n_clusters*100)

# Venn diagram of shared core phages
cluster_group_list <- sapply(unique(cluster_core$sample_type), function(x) {
  cluster_group <- unique(cluster_core$vcontact_cluster[cluster_core$sample_type == x])
  return(cluster_group)
})

names(cluster_group_list) <- unique(cluster_core$sample_type)
cluster_group_list <- lapply(cluster_group_list, function(x) c(x, rep(NA, max(sapply(cluster_group_list, function(y) length(y))) - length(x))))
total_core_clusters <- length(unique(cluster_core$vcontact_cluster))
cluster_overlap <- lapply(calculate.overlap(cluster_group_list), function(x) signif(100*length(x)/total_core_clusters, 2))
tiff("figures/venn_core_clusters.tiff")
draw.quad.venn(area.vector = c(cluster_overlap$a1, cluster_overlap$a2, cluster_overlap$a3, cluster_overlap$a4, 
                               cluster_overlap$a5, cluster_overlap$a6, cluster_overlap$a7, cluster_overlap$a8, 
                               cluster_overlap$a9, cluster_overlap$a10, cluster_overlap$a11, cluster_overlap$a12, 
                               cluster_overlap$a13, cluster_overlap$a14, cluster_overlap$a15), 
               direct.area = T, category = unique(metadata_us$sample_type), 
               fill = cols[unique(metadata_us$sample_type)],
               alpha = rep(0.3, 4), cex = 1.2, margin = 0.1)
dev.off()


# Get noncore of each body site
cluster_noncore <- data.frame()
for (i in 1:length(unique_sample_types)) {
  cluster_noncore <- rbind(cluster_noncore, cluster_all_perc[cluster_all_perc$sample_type == unique_sample_types[i] & cluster_all_perc$perc <= core_cutoff[names(core_cutoff) == unique_sample_types[i]],])
}

# Venn diagram of shared non-core phages
cluster_noncore_group_list <- sapply(unique(cluster_noncore$sample_type), function(x) {
  cluster_noncore_group <- unique(cluster_noncore$vcontact_cluster[cluster_noncore$sample_type == x])
  return(cluster_noncore_group)
})

names(cluster_noncore_group_list) <- unique(cluster_noncore$sample_type)
cluster_noncore_group_list <- lapply(cluster_noncore_group_list, function(x) c(x, rep(NA, max(sapply(cluster_noncore_group_list, function(y) length(y))) - length(x))))
total_noncore_clusters <- length(unique(cluster_noncore$vcontact_cluster))
cluster_overlap <- lapply(calculate.overlap(cluster_noncore_group_list), function(x) signif(100*length(x)/total_noncore_clusters, 2))
tiff("figures/venn_noncore_clusters.tiff")
draw.quad.venn(area.vector = c(cluster_overlap$a1, cluster_overlap$a2, cluster_overlap$a3, cluster_overlap$a4, 
                               cluster_overlap$a5, cluster_overlap$a6, cluster_overlap$a7, cluster_overlap$a8, 
                               cluster_overlap$a9, cluster_overlap$a10, cluster_overlap$a11, cluster_overlap$a12, 
                               cluster_overlap$a13, cluster_overlap$a14, cluster_overlap$a15), 
               direct.area = T, category = unique(metadata_us$sample_type), 
               fill = cols[unique(metadata_us$sample_type)],
               alpha = rep(0.3, 4), cex = 1.2, margin = 0.1)
dev.off()

# Taxonomy breakdown of core phages
cluster_core_demovir <- left_join(cluster_core, vir_counts_longus_demovir, by = c("sample_type", "vcontact_cluster")) %>%
  group_by(ID) %>%
  mutate(V1_prop = V1/sum(V1))

cluster_core_n <- left_join(cluster_core_demovir, cluster_count_tp_summary, by = c("sample_type", "Sample.name")) %>% 
  group_by(sample_type, tp) %>%
  summarise(n_distinct(Sample.name))

tiff("figures/barplot_demovir_core_rel.tiff", width = 5000, height = 2000, res = 400)
ggplot(cluster_core_demovir, aes(x = ID, y = V1_prop, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of mapped reads") + xlab("Sample")
dev.off()

# Taxonomy breakdown of noncore phages
cluster_noncore_demovir <- left_join(cluster_noncore, vir_counts_longus_demovir, by = c("sample_type", "vcontact_cluster")) %>%
  group_by(ID) %>%
  mutate(V1_prop = V1/sum(V1))

cluster_noncore_n <- left_join(cluster_noncore_demovir, cluster_count_tp_summary, by = c("sample_type", "Sample.name")) %>% 
  group_by(sample_type, tp) %>%
  summarise(n_distinct(Sample.name))

tiff("figures/barplot_demovir_noncore_rel.tiff", width = 5000, height = 2000, res = 400)
ggplot(cluster_noncore_demovir, aes(x = ID, y = V1_prop, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of mapped reads") + xlab("Sample")
dev.off()

########## Host range############################
# Percentage of contigs that have CRISPR hosts
sum(!is.na(vir_counts_prop_melt$crispr_host))/nrow(vir_counts_prop_melt) * 100

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
          labRow = as.expression(lapply(rownames(vir_counts_prop_agg2_log), function(a) bquote(italic(.(a))))),
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

# Phage host for each phage family
vir_fam_host <- left_join(vir_counts_prop_melt_agg[vir_counts_prop_melt_agg$Var2 %in% metadata$ID[metadata$Visit_Number == 1],-c("V1")], vir_counts_prop_melt_agg2[,-c("V1")], by = c("Var2", "sample_type", "Location"))
vir_fam_host <- vir_fam_host[!vir_fam_host$demovir %in% c("Other", "Unassigned"),]
vir_fam_host_summary <- vir_fam_host %>% 
  mutate(crispr_host = replace(crispr_host, is.na(crispr_host), "Unassigned")) %>%
  group_by(sample_type, demovir) %>%
  mutate(n = length(crispr_host)) %>%
  group_by(sample_type, demovir, crispr_host, n) %>%
  summarise(num_host = n()) %>%
  mutate(prop = num_host/n) %>%
  ungroup() 

vir_fam_host_top <- vir_fam_host_summary %>%
  group_by(crispr_host) %>%
  summarise(sum_prop = sum(prop)) %>%
  arrange(desc(sum_prop))

vir_fam_host_summary <- vir_fam_host_summary %>%
  mutate(crispr_host_alt = replace(crispr_host, !crispr_host %in% vir_fam_host_top$crispr_host[c(1:12)], "Other"))

crispr_cols <- brewer.pal(12, "Set3")
crispr_cols <- c(crispr_cols, "grey")
names(crispr_cols) <- c(unique(vir_fam_host_summary$crispr_host_alt[vir_fam_host_summary$crispr_host_alt != "Unassigned"]), "Unassigned")

tiff("figures/vir_family_host.tiff", width = 2000)
ggplot(vir_fam_host_summary, aes(demovir, prop, fill = factor(crispr_host_alt))) +
  geom_bar(stat = "identity") +
  facet_grid(~sample_type, scales = "free", space = "free") +
  theme_classic() + xlab("Phage Family") + ylab("Proportion of Phage Clusters") +
  scale_fill_manual("Predicted Phage Host", values = crispr_cols[names(crispr_cols) %in% unique(vir_fam_host_summary$crispr_host_alt)])
dev.off()

########## Microbial composition #########
metaphlan <- read.csv("data/all_metaphlan.csv", stringsAsFactors = FALSE)
row.names(metaphlan) <- metaphlan$X
metaphlan <- metaphlan[,names(metaphlan) != "X"]
metaphlan <- metaphlan[,names(metaphlan) %in% clusters_samples_tsne$ID]
metaphlan <- metaphlan[,names(metaphlan) %in% metadata$ID[metadata$Visit_Number == 1]]
metaphlan_filter <- metaphlan[grepl("g__", row.names(metaphlan)) & !grepl("s__", row.names(metaphlan)),]
metaphlan_filter <- metaphlan_filter[grepl("k__Bacteria", row.names(metaphlan_filter)) | 
                                       grepl("k__Archaea", row.names(metaphlan_filter)),]
metaphlan_filter <- metaphlan_filter[!grepl("unclassified", row.names(metaphlan_filter)),]
metaphlan_filter <- metaphlan_filter[!grepl("noname", row.names(metaphlan_filter)),]
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

# Procrustes analysis
protest_res <- protest(clusters_samples_tsne[,c(1:2)], metaphlan_tsne[,c(1:2)], scale = TRUE)

tiff("figures/procrustes.tiff")
plot(protest_res)
points(protest_res, display = "target", col = "red")
dev.off()

# Heatmap of metaphlan
genera <- unique(contig_data$crispr_host[!is.na(contig_data$crispr_host)])
genera <- genera[!genera %in% c(NA, "Lachnospiraceae", "Bacteroidetes")]
metaphlan_genera <- gsub("\\|.*", "", gsub(".*\\|g__", "", row.names(metaphlan_filter)))
match_ind <- c()
for (i in 1:length(genera)) {
  tmp_match_ind <- grep(paste0(genera[i], "$"), metaphlan_genera)
  match_ind <- c(match_ind, tmp_match_ind)
}
metaphlan_host <- as.matrix(metaphlan_filter[match_ind,])
metaphlan_host_log <- log10(metaphlan_host)
metaphlan_host_log[is.infinite(metaphlan_host_log)] <- -6
rownames(metaphlan_host_log) <- gsub("\\|.*", "", gsub(".*\\|g__", "", rownames(metaphlan_host_log)))

tiff("figures/heatplot_metaphlan.tiff", width = 2000, height = 3000, res = 300)
heatmap.2(metaphlan_host_log,
          margins = c(10,10),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(-6, 2, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA, key.xlab = "log10(rel. abundance)", key.ylab = NA,
          ColSideColors = cols[as.factor(metadata[colnames(metaphlan_host_log), "sample_type"])],
          labRow = as.expression(lapply(rownames(metaphlan_host_log), function(a) bquote(italic(.(a))))),
          labCol = NA,
          cexCol = 0.5,
          #cexRow = 1,
          xlab = "Samples",
          ylab = "Predicted host genera",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(metaphlan_host_log), "sample_type"])),
       col = cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Body sites")
dev.off()

# Plot species of genera that are abundant in all sites but only found in certain sites as phage hosts
metaphlan_abund <- metaphlan[grepl("s__", row.names(metaphlan)) & !grepl("t__", row.names(metaphlan)),]
metaphlan_abund <- metaphlan_abund[grepl("k__Bacteria", row.names(metaphlan_abund)) | 
                                       grepl("k__Archaea", row.names(metaphlan_abund)),]
metaphlan_abund <- metaphlan_abund[!grepl("unclassified", row.names(metaphlan_abund)),]
metaphlan_abund <- metaphlan_abund[!grepl("noname", row.names(metaphlan_abund)),]
metaphlan_abund_genera <- gsub("\\|.*", "", gsub(".*\\|g__", "", row.names(metaphlan_abund)))
abund_genera <- c("Eubacterium", "Haemophilus", "Prevotella", "Streptococcus", "Veillonella")
abund_ind <- c()
for (i in 1:length(abund_genera)) {
  abund_ind <- unique(c(abund_ind, grep(paste0(abund_genera[i], "$"), metaphlan_abund_genera)))
}
metaphlan_abund <- as.matrix(metaphlan_abund[abund_ind,])
metaphlan_abund_log <- log10(metaphlan_abund)
metaphlan_abund_log[is.infinite(metaphlan_abund_log)] <- -6
rownames(metaphlan_abund_log) <- gsub("_", " ", gsub(".*\\|s__", "", rownames(metaphlan_abund_log)))
rownames(metaphlan_abund_log) <- gsub("Candidatus Prevotella conceptionensis", "Prevotella conceptionensis", rownames(metaphlan_abund_log))
metaphlan_abund_log <- metaphlan_abund_log[order(rownames(metaphlan_abund_log)),]

tiff("figures/heatplot_metaphlan_species.tiff", width = 2000, height = 4000, res = 300)
heatmap.2(metaphlan_abund_log,
          margins = c(10,10),
          trace = "none",
          scale = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          # reorderfun = function(d, w=order(d)) reorder(d, w),
          Rowv = FALSE,
          dendrogram = "column",
          col =  c("white", brewer.pal(9, "PuRd")[3:9]),
          breaks = seq(-6, 2, length.out = 9),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,8),
          key.title = NA, key.xlab = "log10(rel. abundance)", key.ylab = NA,
          ColSideColors = cols[as.factor(metadata[colnames(metaphlan_abund_log), "sample_type"])],
          labRow=as.expression(lapply(rownames(metaphlan_abund_log), function(a) bquote(italic(.(a))))),
          labCol = NA,
          cexCol = 0.5,
          #cexRow = 1,
          xlab = "Samples",
          ylab = "Predicted host species",
          main = NA
)
legend(x = 0.85, y = 1.05, xpd=TRUE, legend = levels(as.factor(metadata[colnames(metaphlan_abund_log), "sample_type"])),
       col = cols, bg = "white", box.col = "black",
       lty = 1, lwd = 5, cex = 0.5, title = "Body sites")
dev.off()

# metaphlan_abund <- metaphlan_abund[!grepl("unclassified", row.names(metaphlan_abund)),]
# metaphlan_abund <- metaphlan_abund[!grepl("noname", row.names(metaphlan_abund)),]
metaphlan_abund <- metaphlan_abund[,colSums(metaphlan_abund) != 0]

########## Alpha-diversity and phage richness ####
## Functions ##
createPairedData <- function(df_meta, pair_list){
  
  # Create groups e.g. stool vs dental
  df_meta_pairs <- data.frame()
  for(i in 1:length(pair_list)){
    
    samples_one <- unique(df_meta$Sample.name[df_meta$sample_type == pair_list[[i]][1]])
    samples_two <- unique(df_meta$Sample.name[df_meta$sample_type == pair_list[[i]][2]])
    df_meta_pair <- df_meta[(df_meta$Sample.name %in% Reduce(intersect,list(samples_one, samples_two))) & (df_meta$sample_type %in% pair_list[[i]]),]
    
    df_meta_pair$group <- paste(pair_list[[i]][1], "vs.", pair_list[[i]][2])
    if(nrow(df_meta_pairs) == 0) {
      df_meta_pairs <- df_meta_pair
    } else {
      df_meta_pairs <- rbind(df_meta_pairs, df_meta_pair)
    }
  }
  
  # Order by Location and change characters in group
  df_meta_pairs <- df_meta_pairs[order(df_meta_pairs$Location),]
  df_meta_pairs$group <- as.character(df_meta_pairs$group)
  
  return(df_meta_pairs)
}

runTtest <- function(df_paired){
  
  df_paired_richness <- df_paired[!duplicated(paste0(df_paired$ID, df_paired$group)),]
  
  # T-test
  p_values <- c()
  ttest_groups <- df_paired_richness[!duplicated(paste0(df_paired_richness$Location,
                                                        df_paired_richness$group)),]
  for(i in 1:nrow(ttest_groups)){
    y <- df_paired_richness[df_paired_richness$Location == ttest_groups$Location[i] & 
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
  ttest_groups <- ttest_groups[order(ttest_groups$Location, ttest_groups$group),]
  ttest_groups <- ttest_groups[,c("Location", "group", "pvalue", "asterisk")]
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
    ggtitle(ttest_group$Location) +
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
                                                     df_paired_richness$group == ttest_groups$group[i],]
    g[[i]] <- plotRichnessGraph(df_paired_richness_group, ttest_groups[i,], cols) +
      theme(legend.position = "none") + ylim(c(0, max(df_paired_richness$richness) + 50))
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
vir_counts_prop_melt_meta <- left_join(vir_counts_prop_melt, metadata, by = c("Var2"="ID","sample_type","Location")) %>% rename(ID = Var2)
vir_counts_prop_melt_meta <- vir_counts_prop_melt_meta[vir_counts_prop_melt_meta$Visit_Number == 1,]

# Sample - phage cluster matrix
vir_cluster_counts <- dcast(vir_counts_prop_melt_meta, ID ~ vcontact_cluster, length)
rownames(vir_cluster_counts) <- vir_cluster_counts[,1]
vir_cluster_counts <- vir_cluster_counts[,-1]

# # Remove samples with lower than quartile for
# remove_ids <- rownames(vir_cluster_counts)[rowSums(vir_cluster_counts > 0) <= 3]
# vir_cluster_counts <- vir_cluster_counts[!rownames(vir_cluster_counts) %in% remove_ids,]
# metadata_richness <- metadata[!metadata$ID %in% unique(c(remove_ids, metadata$ID[!metadata$ID %in% vir_counts_prop_melt_meta$ID])),]

pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
paired_metadata <- createPairedData(metadata[metadata$ID %in% vir_counts_prop_melt_meta$ID,], pair_list)
paired_metadata_summary <- paired_metadata %>% group_by(Location, group, sample_type) %>% summarise(n())

# Get richness from matrix
richness_paired <- data.frame(ID = rownames(vir_cluster_counts), richness = rowSums(vir_cluster_counts > 0), no_phages = rowSums(vir_cluster_counts)) %>%
  right_join(paired_metadata, by = "ID")

# Alpha-diversity vs. number of contigs
richness <- data.frame(ID = rownames(vir_cluster_counts), richness = rowSums(vir_cluster_counts > 0), no_phages = rowSums(vir_cluster_counts)) %>%
  right_join(metadata[metadata$ID %in% unique(paired_metadata$ID),], by = "ID")
richness$n_contigs <- sapply(richness$ID, function(x) length(which(vir_counts_prop_melt$Var2 == x)))

linear_mod <- lm(richness ~ n_contigs, richness)
summary(linear_mod)
tiff("figures/richness_ncontigs.tiff")
plot(richness ~ n_contigs, richness, pch = 16, xlab = "No. contigs", ylab = "Phage Cluster Richness")
abline(linear_mod, col = "red")
dev.off()

# For each group, remove samples with less than or equal to 100 phage contigs
unique_groups <- unique(richness_paired$group)
for (i in 1:length(unique_groups)) {
  remove_samples <- richness_paired$Sample.name[(richness_paired$no_phages <= 100 & richness_paired$group %in% unique_groups[i])]
  richness_paired <- richness_paired[!(richness_paired$Sample.name %in% remove_samples & richness_paired$group %in% unique_groups[i]),]
}

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
lay <- rbind(c(4,5,6,1,2,3), c(7,8,9,10,10,10))
tiff("figures/alpha_diversity_subsampled.tiff", width = 4000, height = 2000, res = 220)
grid.arrange(grobs = richness_graphs_ss, layout_matrix = lay)
dev.off()


########## Size of phages ####
contig_data_meta <- left_join(vir_counts_prop_melt, metadata, by = c("Var2"="ID","sample_type","Location")) %>% rename(ID = Var2) %>%
  left_join(contig_data[,c("name", "circular", "size")], by = c("Var1"="name"))

countSizeCat <- function(contig_data_meta, size_lower, size_upper) {
  size_summary <- contig_data_meta %>%
    group_by(sample_type) %>%
    mutate(n_samples = n_distinct(ID)) %>% ungroup() %>%
    filter(circular, size < size_upper, size >= size_lower) %>%
    group_by(sample_type, n_samples) %>%
    summarise(n = n_distinct(vcontact_cluster)) %>%
    mutate(n_per_sample = n/n_samples) %>%
    mutate(size_lower = size_lower, category = paste(as.character(size_lower/1e3), "<=\n<", as.character(size_upper/1e3)))
  return(size_summary)
}

sizes_lower <- c(3e3, 1e4, 5e4, 2e5)
sizes_upper <- c(sizes_lower[-1], Inf)
phage_size_summary <- map2_df(.x = sizes_lower, .y = sizes_upper, .f = function(.x, .y) countSizeCat(contig_data_meta, .x, .y))
phage_size_summary$category <- gsub("\n< Inf", "", phage_size_summary$category)

tiff("figures/circular_phage_size.tiff", width = 800, height = 500, res = 150)
ggplot(phage_size_summary, aes(reorder(category, size_lower), n_per_sample, fill = sample_type)) +
  geom_bar(stat = "identity") + xlab("Size (kb)") + ylab("No. unqiue circular phage clusters \n/ no. samples") + 
  theme_classic() +
  scale_fill_manual("Body Site", values = cols)
dev.off()

# Size of all contigs
scaffold_lengths <- readRDS("data/scaffold_lengths.RDS")
scaffold_lengths_summary <- left_join(scaffold_lengths, metadata, by = "ID") %>%
  group_by(sample_type) %>%
  summarise("3 <\n<= 10" = sum(size < sizes_upper[1] & size >= sizes_lower[1])/length(size)*100,
            "10 <\n<= 50" = sum(size < sizes_upper[2] & size >= sizes_lower[2])/length(size)*100,
            "50 <\n<= 200" = sum(size < sizes_upper[3] & size >= sizes_lower[3])/length(size)*100,
            "200 <" = sum(size < sizes_upper[4] & size >= sizes_lower[4])/length(size)*100) %>%
  ungroup() %>%
  melt(variable.name = "perc_category", value.names = "perc", id.vars = "sample_type")

tiff("figures/contigs_size.tiff", width = 3000, height = 1000, res = 200)
ggplot(scaffold_lengths_summary, aes(perc_category, value, fill = sample_type)) +
  geom_bar(stat="identity") +
  facet_grid(~sample_type, scale = "free", space = "free") +
  scale_fill_manual("Body Site", values = cols, labels = names(cols)) +
  theme_classic() + xlab("Size (kb)") + ylab("% All Contigs") + ylim(c(0,100))
dev.off()

########## Jumbo/mega circular phages ####
jumbophage_contigs_meta <- left_join(vir_counts_prop_melt, metadata, by = c("Var2"="ID","sample_type","Location")) %>% rename(ID = Var2) %>%
  left_join(jumbophage_contigs[,c("name", "circular", "size")], by = c("Var1"="name"))
jumbophage_contigs_meta <- jumbophage_contigs_meta[!is.na(jumbophage_contigs_meta$circular),]
jumbophage_contigs_meta <- jumbophage_contigs_meta[jumbophage_contigs_meta$circular,]
jumbophage_contigs_meta$numeric_samplename <- as.numeric(factor(jumbophage_contigs_meta$Sample.name))
jumbophage_contigs_meta$numeric_samplename[!jumbophage_contigs_meta$numeric_samplename %in% unique(jumbophage_contigs_meta$numeric_samplename[duplicated(jumbophage_contigs_meta$numeric_samplename)])] <- NA
jumbophage_contigs_meta$numeric_samplename <- as.numeric(factor(jumbophage_contigs_meta$numeric_samplename))
jumbophage_contigs_meta$numeric_samplename[is.na(jumbophage_contigs_meta$numeric_samplename)] <- c((max(jumbophage_contigs_meta$numeric_samplename, na.rm = TRUE)+1):(max(jumbophage_contigs_meta$numeric_samplename, na.rm = TRUE)+sum(is.na(jumbophage_contigs_meta$numeric_samplename))))
jumbophage_contigs_meta$vcontact_cluster[jumbophage_contigs_meta$vcontact_cluster %in% ""] <- NA

# Percentage of cohorts containing a jumbophage
jumbophage_summary <- jumbophage_contigs_meta %>%
  filter(Visit_Number == 1) %>%
  group_by(Location, sample_type) %>%
  summarise(n_samples_jumbophage = n_distinct(ID))
metadata_summary <- metadata %>%
  filter(Visit_Number == 1) %>%
  group_by(Location, sample_type) %>%
  summarise(n_samples = n_distinct(ID))
jumbophage_summary <- full_join(jumbophage_summary, metadata_summary, by = c("Location", "sample_type")) %>%
  mutate(perc_samples = 100*n_samples_jumbophage/n_samples, Location_sampletype = paste(sample_type, Location, sep = "\n"))

tiff("figures/jumbophage_percentage_samples.tiff", width = 1000, height = 500, res = 100)
ggplot(jumbophage_summary, aes(Location_sampletype, perc_samples)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Cohort") + ylab("% samples") + ylim(c(0,100))
dev.off()

# Transient/Persistent jumbophages
jumbophage_contigs_meta$vcontact_cluster <- as.character(jumbophage_contigs_meta$vcontact_cluster)
jumbophage_contigs_meta$vcontact_cluster_Var1 <- jumbophage_contigs_meta$vcontact_cluster
jumbophage_contigs_meta$vcontact_cluster_Var1[is.na(jumbophage_contigs_meta$vcontact_cluster)] <- jumbophage_contigs_meta$Var1[is.na(jumbophage_contigs_meta$vcontact_cluster)]

jumbophages_long <- jumbophage_contigs_meta %>% 
  group_by(vcontact_cluster_Var1, Location_sampletype, Sample.name) %>%
  summarise(Timepoints = n_distinct(Visit_Number)) %>%
  ungroup() %>% group_by(vcontact_cluster_Var1, Location_sampletype) %>%
  select(-Sample.name) %>%
  arrange(Timepoints) %>%
  summarise_each(funs(paste(Timepoints, collapse = ", "))) 
  
jumbophage_persistent <- jumbophage_contigs_meta %>% 
  group_by(vcontact_cluster_Var1, Location_sampletype, Sample.name) %>%
  summarise(Timepoints = n_distinct(Visit_Number)) %>%
  ungroup() %>% group_by(vcontact_cluster_Var1, Location_sampletype) %>%
  select(-Sample.name) %>%
  filter(Timepoints > 2)

jumbophage_us <- jumbophage_contigs_meta %>%
  filter(Location == "US") %>%
  summarise(n = n_distinct(Sample.name))

# Remove phage duplicates
jumbophage_contigs_meta_nodup <- jumbophage_contigs_meta[!duplicated(paste(jumbophage_contigs_meta$Var1, jumbophage_contigs_meta$sample_type)),]

# label_colours <- c("turquoise", "darkgreen")
# pd <- position_dodge(0.8)
# h <- 200
# tiff("figures/jumbophage_individual.tiff", width = 2500, height = 1000, res = 200)
# ggplot(jumbophage_contigs_meta, aes(sample_type, size/1000, colour = Location, group = numeric_samplename)) +
#   geom_line(position=pd, colour="grey90", linetype = "dashed") +
#   geom_point(position=pd, alpha = 0.6) +
#   geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
#   theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
#   scale_y_continuous(breaks = seq(200, 800, 50)) +
#   theme(panel.grid.major.y = element_line(colour = "grey90"), 
#         axis.text.x = element_text(size = 12), 
#         axis.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 16)) +
#   scale_colour_manual("Country", values = label_colours) 
# dev.off()
# 
# tiff("figures/jumbophage_cluster.tiff", width = 2500, height = 1000, res = 200)
# ggplot(jumbophage_contigs_meta, aes(sample_type, size/1000, colour = Location, group = vcontact_cluster)) +
#   geom_line(position=pd, colour="grey90", linetype = "dashed") +
#   geom_point(position = pd, alpha = 0.6) +
#   geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
#   theme_classic() + xlab("Body Site") + ylab("Circular contig size (kb)") +
#   scale_y_continuous(breaks = seq(200, 800, 50)) +
#   theme(panel.grid.major.y = element_line(colour = "grey90"), 
#         axis.text.x = element_text(size = 12), 
#         axis.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 16)) +
#   scale_colour_manual("Country", values = label_colours) 
# dev.off()

jumbophage_contigs_meta_nodup$sampletype_Location <- paste(jumbophage_contigs_meta_nodup$sample_type, "-", jumbophage_contigs_meta_nodup$Location)
h = 200
tiff("figures/circular_jumbophages.tiff", width = 2000, height = 1000, res = 250)
ggplot(jumbophage_contigs_meta_nodup, aes(x = "Contigs", y = size/1000, colour = sampletype_Location)) +
  geom_jitter(alpha = 0.6, size = 2) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  theme_classic() + ylab("Circular contig size (kb)") + xlab("") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Body Site - Geographical Location", values = cohort_cols)
dev.off()

# Number of samples
jumbophage_samples <- jumbophage_contigs_meta %>%
  select(size, demovir, crispr_host, Location_sampletype, vcontact_cluster_Var1, Sample.name) %>%
  group_by(Location_sampletype, vcontact_cluster_Var1) %>%
  summarise(num_samples = n_distinct(Sample.name)) %>% ungroup() 
  
# Jumbophage table
jumbophage_table <- jumbophage_contigs_meta %>%
  select(size, demovir, crispr_host, Location_sampletype, vcontact_cluster_Var1) %>%
  group_by(demovir, crispr_host, Location_sampletype, vcontact_cluster_Var1) %>%
  summarise_each(funs(paste(unique(size), collapse = "\n"))) %>%
  left_join(jumbophage_samples, by = c("Location_sampletype", "vcontact_cluster_Var1")) %>%
  left_join(jumbophages_long, by = c("Location_sampletype", "vcontact_cluster_Var1")) %>%
  mutate(num_Location_sampletype = paste0(num_samples, " ", Location_sampletype, " (", Timepoints, ")")) %>%
  ungroup() %>%
  select(-num_samples, -Location_sampletype) %>%
  group_by(vcontact_cluster_Var1) %>%
  summarise_each(funs(paste(unique(.), collapse = "\n")))  %>%
  select(vcontact_cluster_Var1, size, demovir, crispr_host, num_Location_sampletype) 

jumbophage_table$vcontact_cluster_Var1[!grepl("VC", jumbophage_table$vcontact_cluster_Var1)] <- "No Phage Cluster"
names(jumbophage_table) <- c("Phage Cluster", "Size (nt)", "Phage Family", "Predicted Host Genus", "No. sites (Timepoints)")

write.csv(jumbophage_table, "data/circular_jumbophage_summary_table.csv")

######### Function of genes ####

# Summary of functional categories
jumbophage_gff <- read.csv("data/jumbophage_prodigal.gff", stringsAsFactors = FALSE)
jumbophage_gff_meta <- inner_join(jumbophage_gff, jumbophage_contigs_meta, by = c("V1" = "Var1"))
jumbophage_gff_meta$vcontact_cluster <- as.character(jumbophage_gff_meta$vcontact_cluster)
jumbophage_gff_meta$vcontact_cluster[is.na(jumbophage_gff_meta$vcontact_cluster)] <- "No Phage Cluster"
jumbophage_gff_meta$Parent[is.na(jumbophage_gff_meta$Parent)] <- "Hypothetical protein"
jumbophage_gff_meta$size <- as.character(jumbophage_gff_meta$size)

cat_summary <- jumbophage_gff_meta %>% group_by(size, vcontact_cluster, Parent) %>%
  summarise(n_cog = n()) %>%
  group_by(size, vcontact_cluster) %>%
  mutate(sum_n_cog = sum(n_cog)) %>%
  mutate(per_cog = n_cog/sum_n_cog * 100)

unique_cat <- unique(cat_summary$Parent)
cat_summary$Parent <- factor(cat_summary$Parent, levels = c("Hypothetical protein", sort(unique_cat[unique_cat != "Hypothetical protein"])))
cat_cols <- c("white", "#e0a8e0", "#d6a738", "#4d4cd4", "#6ed4c0", "#f2514c", "#b64ed9")
names(cat_cols) <- c("Hypothetical protein", sort(unique_cat[unique_cat != "Hypothetical protein"]))

# Summary of functional genes
jumbophages_func_summary <- jumbophage_gff_meta %>%
  filter(Parent != "Hypothetical protein") %>%
  group_by(size, Parent, Name) %>%
  summarise(No_genes = n()) 
names(jumbophages_func_summary) <- c("Size (nt)", "Category", "Protein", "No. Genes")
write.csv(jumbophages_func_summary, "data/circular_jumbophage_functional_summary.csv")

# Functional genes of persistent phage clusters
unique_persistent_function <- unique(jumbophage_gff_meta$protein_description[jumbophage_gff_meta$vcontact_cluster_Var1 %in% jumbophage_persistent$vcontact_cluster_Var1])
unique_transient_function <- unique(jumbophage_gff_meta$protein_description[!jumbophage_gff_meta$vcontact_cluster_Var1 %in% jumbophage_persistent$vcontact_cluster_Var1])
exclusive_pers_func <- unique_persistent_function[!unique_persistent_function %in% unique_transient_function]
excl_pers_func_summary <- jumbophage_gff_meta %>% 
  filter(vcontact_cluster_Var1 %in% jumbophage_persistent$vcontact_cluster_Var1) %>%
  filter(protein_description %in% exclusive_pers_func) %>%
  group_by(protein_description) %>%
  summarise(n = n_distinct(vcontact_cluster_Var1))

tiff("figures/functional_categories_cluster.tiff", width = 6500, height = 1000, res = 150)
ggplot(cat_summary, aes(size, per_cog, fill = Parent)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  facet_grid(~ vcontact_cluster, space = "free", scale = "free") +
  scale_fill_manual("Functional Categories", values = cat_cols, labels = names(cat_cols)) +
  xlab("Size (kb)") + ylab("Percentage") + theme_classic()
dev.off()

# GFF for Gview
jumbophage_gff$V9 <- ifelse(is.na(jumbophage_gff$Name),
                          paste0(jumbophage_gff$V9, "Name=hypothetical protein;"),
                          paste0(jumbophage_gff$V9, "Name=", gsub(" \\[.*", "", jumbophage_gff$Name), ";"))
jumbophage_gff$V9 <- ifelse(is.na(jumbophage_gff$Parent),
                          paste0(jumbophage_gff$V9, "Parent=None;"),
                          paste0(jumbophage_gff$V9, "Parent=", jumbophage_gff$Parent, ";"))
jumbophage_gff <- jumbophage_gff[,c(1:9)]
names(jumbophage_gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Add tRNA
jumbophages_trna <- read.table("data/megaphage_trna_clean.tab", stringsAsFactors = FALSE)
jumbophages_trna$strand <- ifelse(grepl("c", jumbophages_trna$V3), "-", "+")
jumbophages_trna <- cbind(jumbophages_trna, map_df(.x = gsub("c", "", gsub("\\]", "", gsub("\\[", "", jumbophages_trna$V3))),
                                                 function(.x) {
                                                   split <- strsplit(.x, ",")
                                                   return(data.frame(start = as.numeric(split[[1]][1]),
                                                                     end = as.numeric(split[[1]][2])))
                                                 }))
jumbophages_trna$attributes <- paste0("Name=", jumbophages_trna$V2, ";")
trna_gff <- jumbophages_trna %>% mutate(source = "ARAGORN_v1.2.36", type = "tRNA", phase = 0, score = format(as.numeric(V4), nsmall=1)) %>%
  select(seqid = V6, source, type, start, end, score, strand, phase, attributes)

# Combine prodigal CDS and trna
comb_gff <- rbind(jumbophage_gff, trna_gff)

# Write GFF file for largest two circular phages
circular_jumbophages <- unique(jumbophage_contigs_meta$Var1[order(jumbophage_contigs_meta$size, decreasing = TRUE)])[c(1,2)]
for (i in 1:length(circular_jumbophages)) {
  write.table(comb_gff[comb_gff$seqid == circular_jumbophages[i],], paste0("data/", circular_jumbophages[i], ".gff"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}

########## Phages and ARGs 
# Read data
args_data <- read.csv("data/all_assemblies_card_90.csv", stringsAsFactors = FALSE)
args_data$ID <- gsub("_card.out", "", gsub(".*/", "", args_data$filename))
args_data$name <- paste0(args_data$ID, "_", args_data$qseqid)

# Combine phage and arg results where only phage does/doesnt contain ARG
jumbophages_args <- left_join(jumbophage_contigs_meta, args_data, by = c("Var1"="name", "ID"), suffix = c(".phage", ".arg"))
jumbophages_args <- jumbophages_args[!duplicated(jumbophages_args$Var1),]

# Percentage jumbophages containing ARG
length(phages_args_yn$Var1[!is.na(phages_args_yn$ARO.Name)])*100/length(phages_args_yn$Var1[is.na(phages_args_yn$ARO.Name)])

######## EXTRA #########
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
jumbophage_prots <- read.delim("data/megaphage_proteins_headers.txt", stringsAsFactors = FALSE, header = FALSE)
names(jumbophage_prots) <- "query_name"
eggnog <- read.delim("data/megaphage_annot.emapper.annotations", stringsAsFactors = FALSE, header = FALSE)
names(eggnog) <- c("query_name", "seed_eggnog_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "predicted_taxonomic_group",
                   "predicted_protein_name", "gene_ontology_terms", "EC_number", "KEGG_ko", "KEGG_pathway", "KEGG_module", "KEGG_reaction",
                   "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_reaction", "tax_scope", "eggnog_ogs", "bestOG", "COG_functional_cat", "eggnog_description")

jumbophage_eggnog <- left_join(jumbophage_prots, eggnog, by = "query_name")
jumbophage_eggnog$name <- sub("_[^_]+$", "", jumbophage_eggnog$query_name)
jumbophage_eggnog_meta <- left_join(jumbophage_contigs_meta, jumbophage_eggnog, by = c("Var1"="name"))

# # Summarise proportion of bacterial/viral proteins
# jumbophage_egng_summary <- jumbophage_eggnog_meta %>% group_by(Var1, tax_scope, size, Location, sample_type) %>%
#   summarise(tax_n = n()) %>% group_by(Var1, size, Location, sample_type) %>%
#   mutate(total_n = sum(tax_n)) %>%
#   mutate(tax_perc = tax_n/total_n * 100)
# jumbophage_func <- jumbophage_egng_summary[!is.na(jumbophage_egng_summary$tax_scope),]
# jumbophage_nonviral <- jumbophage_func[jumbophage_func$tax_scope != "Viruses",]

# Summarise functional groups for each cohort
jumbophage_eggnog_meta <- jumbophage_eggnog_meta[!is.na(jumbophage_eggnog_meta$seed_eggnog_ortholog),]
jumbophage_eggnog_meta$COG_functional_cat[jumbophage_eggnog_meta$COG_functional_cat == ""] <- "S"
jumbophage_eggnog_cogs <- jumbophage_eggnog_meta %>%
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
jumbophage_eggnog_cogs$functional_cats <- sapply(jumbophage_eggnog_cogs$COG_functional_cat, function(x) functional_cats[names(functional_cats) == x])
jumbophage_eggnog_cogs$vcontact_cluster <- as.character(jumbophage_eggnog_cogs$vcontact_cluster)
jumbophage_eggnog_cogs$vcontact_cluster[is.na(jumbophage_eggnog_cogs$vcontact_cluster)] <- "No Phage Cluster"
jumbophage_eggnog_cogs$size <- as.character(jumbophage_eggnog_cogs$size)

metabolic_summary <- jumbophage_eggnog_cogs %>% group_by(size, vcontact_cluster, functional_cats) %>%
  summarise(n_cog = n()) %>%
  group_by(size, vcontact_cluster) %>%
  mutate(sum_n_cog = sum(n_cog)) %>%
  mutate(per_cog = n_cog/sum_n_cog * 100)

# Plot functional groups for each body site and cohort
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
functional_colours <- c(rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
names(functional_colours) <- functional_cats
functional_colours <- functional_colours[!is.na(names(functional_colours))]
functional_colours[names(functional_colours) == "Function unknown"] <- "white"

tiff("figures/functional_categories_bodysite.tiff", width = 2600, height = 1000, res = 150)
ggplot(metabolic_summary, aes(size, per_cog, fill = functional_cats)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  facet_grid(~ vcontact_cluster, space = "free", scale = "free") +
  scale_fill_manual("Functional Categories", values = functional_colours) +
  xlab("Size (kb)") + ylab("Percentage") + theme_classic()
dev.off()

# Summarise functional groups for each cluster
metabolic_cluster_summary <- jumbophage_eggnog_cogs %>% group_by(genus, functional_cats) %>%
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
jumbophage_cogs <- acast(jumbophage_eggnog_cogs, name ~ seed_eggnog_ortholog)
cogs_dist <- dist(jumbophage_cogs, "binary")
set.seed(1)
mds <- wcmdscale(cogs_dist, w=rep(1,nrow(jumbophage_cogs)))
mds <- as.data.frame(mds[,1:2])
mds$name <- row.names(mds)
mds <- left_join(mds, jumbophage_contigs, by = "name") %>%
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
