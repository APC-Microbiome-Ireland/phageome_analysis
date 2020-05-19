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
library(randomForest)
library(rfPermute)
set.seed(1)

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
  g <- ggplot(df_paired_richness_group, aes(sample_type, richness)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.8) +
    theme_classic() +
    ylab("Phage Cluster Richness") +
    xlab("") +
    ggtitle(ttest_group$Location) +
    geom_text(data = ttest_group, aes(label=asterisk),
              x = 1.5, y = max(df_paired_richness_group$richness)+10, size = 7,
              inherit.aes = FALSE) +
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
scaffold_lengths <- readRDS("data/scaffold_lengths.RDS")

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
metadata_summary <- metadata %>% filter(Visit_Number == 1) %>% group_by(Location, sample_type) %>% summarise(n = n())

# Number of unique phage contigs
n_phage_contigs <- length(unique(vir_counts_prop_melt$Var1))

# Number of phage clusters
n_phage_clusters <- length(unique(vir_counts_prop_melt$vcontact_cluster[!(is.na(vir_counts_prop_melt$vcontact_cluster) | vir_counts_prop_melt$vcontact_cluster %in% "")]))

# Number of phage contigs in or not in clusters 
n_phage_contigs_in_clusters <- length(unique(vir_counts_prop_melt$Var1[!(is.na(vir_counts_prop_melt$vcontact_cluster) | vir_counts_prop_melt$vcontact_cluster %in% "")]))
n_phage_contigs_not_in_clusters <- length(unique(vir_counts_prop_melt$Var1[is.na(vir_counts_prop_melt$vcontact_cluster) | vir_counts_prop_melt$vcontact_cluster %in% ""]))


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
demovir_cols[demovir_cols == "#D9D9D9"] <- "#f7b4f1"

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

# Compute distance matrix and run NMDS on phage clusters
vir_counts_prop_agg_meta <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% metadata$ID[metadata$Visit_Number == 1]]

cluster_samples_nmds <- metaMDS(t(vir_counts_prop_agg_meta), distance = "bray", k = 2, trymax = 20)
df_cluster_samples_nmds <- as.data.frame(cluster_samples_nmds$points)

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 3:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], k = k), df_cluster_samples_nmds[,c("MDS1", "MDS2")])
  avg_sil[k-1] <- mean(tmp[,3])
}

# Group by silhouette width
samples_clust <- pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], which.max(avg_sil)+1)
df_cluster_samples_nmds$cluster = as.factor(samples_clust$cluster[row.names(df_cluster_samples_nmds)])
df_cluster_samples_nmds$ID <- row.names(df_cluster_samples_nmds)
df_cluster_samples_nmds$Sample.name <- as.character(sapply(df_cluster_samples_nmds$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_cluster_samples_nmds$Location <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Location[metadata$ID == x])
df_cluster_samples_nmds$sample_type <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$sample_type[metadata$ID == x])
df_cluster_samples_nmds$Location_sampletype <- paste(df_cluster_samples_nmds$sample_type, "-", df_cluster_samples_nmds$Location)
df_cluster_samples_nmds$Age <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Age[metadata$ID == x])
df_cluster_samples_nmds$Sex <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Sex[metadata$ID == x])

# Plot NMDS by sample type and cluster
tiff("figures/nmds_clusters.tiff", width = 1600, height = 1000, res = 200)
ggplot(df_cluster_samples_nmds, aes(MDS1, MDS2, colour = Location_sampletype, shape = cluster)) +
  theme_classic() +
  geom_text(aes(label=cluster)) +
  scale_colour_manual("Body Site", values = cohort_cols, guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("NMDS 1") + ylab("NMDS 2")
dev.off()

# Calculate proportion of samples
cluster_res <- df_cluster_samples_nmds %>% 
  group_by(cluster, Location, sample_type) %>% 
  summarise(n = n()) %>%
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

# Label by sex and age
tiff("figures/nmds_age_sex.tiff", width = 1600, height = 1000, res = 200)
ggplot(df_cluster_samples_nmds, aes(x = MDS1, y = MDS2, fill = Age, shape = Sex)) +
  theme_classic() + geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("NMDS 1") + ylab("NMDS 2")
dev.off()

# Breakdown of groups
tiff("figures/nmds_group_breakdown.tiff", width = 1600, height = 1000, res = 200)
ggplot(cluster_res, aes(cluster, prop_cluster, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  theme_classic() + xlab("Group") + ylab("% Samples") +
  scale_fill_manual("Body Site - Geographical Location", values = cohort_cols)
dev.off()

########## Differentially abundant clusters#####################
clusters_samples_kw = apply(vir_counts_prop_agg_meta, 1, function(x) kruskal.test(x, df_cluster_samples_nmds$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
vir_counts_prop_agg_meta_diff = vir_counts_prop_agg_meta[rownames(vir_counts_prop_agg_meta) %in% names(which(clusters_samples_kw < 0.001)),]

viral_clusters_df = data.frame(row.names = rownames(vir_counts_prop_agg_meta),
                               demovir = sapply(rownames(vir_counts_prop_agg_meta), 
                                                function(x) names(which.max(table(contig_data[which(contig_data$vcontact_cluster == x),"demovir"])))))

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
unique_clusters <- unique(df_cluster_samples_nmds$cluster)

vir_group_list <- sapply(unique_clusters, function(x) {
  vir_tmp <- vir_counts_prop_agg[,colnames(vir_counts_prop_agg) %in% df_cluster_samples_nmds$ID[df_cluster_samples_nmds$cluster == x]]
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
metadata_longus <- metadata[metadata$Location == "US",] %>% 
  group_by(Sample.name, sample_type) %>%
  filter(any(Visit_Number == 3))

metadata_longus_summary <- metadata_longus %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

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

# Select persistent contigs (3 at least, 2 doesnt converge)
vir_counts_longus_demovir <- left_join(metadata_longus[,c("ID", "Sample.name", "Visit_Number")], vir_counts_prop_melt_agg, by = c("ID"="Var2"))

persist_graphs <- list()
transient_graphs <- list()
n_samples <- rep(NA, length(unique(vir_counts_longus_demovir$sample_type)))
n_persist_clusters <- rep(NA, length(unique(vir_counts_longus_demovir$sample_type)))
n_transient_clusters <- rep(NA, length(unique(vir_counts_longus_demovir$sample_type)))
permanova_persist <- list()
permanova_transient <- list()
cluster_persist_all <- data.frame()
cluster_transient_all <- data.frame()
cluster_sharing_stats <- data.frame()

for (i in 1:length(unique(vir_counts_longus_demovir$sample_type))) {
  cluster_persist <- left_join(vir_counts_longus_demovir, cluster_no_tp, by = c("sample_type", "Sample.name", "vcontact_cluster")) %>%
    filter(no_tp >= 3 & sample_type == unique(vir_counts_longus_demovir$sample_type)[i]) %>%
    group_by(ID) %>%
    filter(n_distinct(vcontact_cluster) > 1)
  
  # Select transient clusters
  cluster_transient <- left_join(vir_counts_longus_demovir, cluster_no_tp, by = c("sample_type", "Sample.name", "vcontact_cluster")) %>%
    filter(no_tp < 3 & sample_type == unique(vir_counts_longus_demovir$sample_type)[i]) %>%
    group_by(ID) %>%
    filter(n_distinct(vcontact_cluster) > 1)
  
  # Same individuals in transient and persistent
  cluster_persist <- cluster_persist %>%
    filter(Sample.name %in% intersect(unique(cluster_transient$Sample.name), unique(cluster_persist$Sample.name)))
  
  # Wilcoxon rank test of shared clusters
  sample_names <- unique(cluster_persist$Sample.name)
  cluster_persist_sharing <- rep(NA, length(sample_names))
  cluster_transient_sharing <- rep(NA, length(sample_names))
  for (j in 1:length(sample_names)) {
    cluster_persist_sample <- unique(cluster_persist$vcontact_cluster[cluster_persist$Sample.name == sample_names[j]])
    cluster_persist_others <- unique(cluster_persist$vcontact_cluster[cluster_persist$Sample.name != sample_names[j]])
    cluster_persist_sharing[j] <- sum(cluster_persist_sample %in% cluster_persist_others)/length(cluster_persist_sample)*100 
    
    cluster_transient_sample <- unique(cluster_transient$vcontact_cluster[cluster_transient$Sample.name == sample_names[j]])
    cluster_transient_others <- unique(cluster_transient$vcontact_cluster[cluster_transient$Sample.name != sample_names[j]])
    cluster_transient_sharing[j] <- sum(cluster_transient_sample %in% cluster_transient_others)/length(cluster_transient_sample)*100 
  }
  cluster_sharing_stats <- rbind(cluster_sharing_stats, data.frame(persist_med = median(cluster_persist_sharing), transient_med = median(cluster_transient_sharing),
             persist_iqr = IQR(cluster_persist_sharing), transient_iqr = IQR(cluster_transient_sharing),
             pval = wilcox.test(cluster_persist_sharing, cluster_transient_sharing)$p.value, sample_type = unique(vir_counts_longus_demovir$sample_type)[i]))
  
  cluster_persist_all <- rbind(cluster_persist_all, as.data.frame(cluster_persist))
  
  cluster_transient <- cluster_transient %>%
    filter(Sample.name %in% intersect(unique(cluster_transient$Sample.name), unique(cluster_persist$Sample.name)))
  
  cluster_transient_all <- rbind(cluster_transient_all, as.data.frame(cluster_transient))
  
  # Cast for persisent NMDS 
  cluster_persist_cast <- dcast(cluster_persist, ID ~ vcontact_cluster, sum, value.var = "V1") 
  rownames(cluster_persist_cast) <- cluster_persist_cast$ID
  cluster_persist_cast <- cluster_persist_cast[,names(cluster_persist_cast) != "ID"]
  
  # Run NMDS
  set.seed(1)
  nmds_persist <- metaMDS(cluster_persist_cast, distance = "bray", k = 2, trymax = 20)
  df_nmds_persist <- as.data.frame(nmds_persist$points)
  df_nmds_persist$ID <- row.names(df_nmds_persist)
  df_nmds_persist$Sample.name <- as.character(sapply(df_nmds_persist$ID, function(x) metadata$Sample.name[metadata$ID == x]))
  df_nmds_persist$Location <- sapply(df_nmds_persist$ID, function(x) metadata$Location[metadata$ID == x])
  
  persist_sample <- data.frame(Sample.name = sapply(rownames(cluster_persist_cast), function(x) unique(cluster_persist$Sample.name[cluster_persist$ID == x]))) 
  permanova_persist[[i]] <- adonis(cluster_persist_cast ~ Sample.name, data = persist_sample)
  
  # Summarise number of samples and clusters
  n_samples[i] <- length(unique(df_nmds_persist$Sample.name))
  n_persist_clusters[i] <- length(unique(cluster_persist$vcontact_cluster))
  
  # Plot NMDS
  persist_graphs[[i]] <- ggplot(df_nmds_persist, aes(MDS1, MDS2)) +
    geom_point(alpha = 0.5, colour = cols[unique(vir_counts_longus_demovir$sample_type)[i]]) +
    geom_encircle(aes(fill = Sample.name), alpha=0.3, expand = 0) +
    theme_classic() +
    guides(fill = FALSE) +
    scale_fill_manual(values = rep(cols[unique(vir_counts_longus_demovir$sample_type)[i]], length(unique(df_nmds_persist$Sample.name)))) +
    ggtitle(unique(vir_counts_longus_demovir$sample_type)[i])
  
  # Cast for transient NMDS
  cluster_transient_cast <- dcast(cluster_transient, ID ~ vcontact_cluster, sum, value.var = "V1") 
  rownames(cluster_transient_cast) <- cluster_transient_cast$ID
  cluster_transient_cast <- cluster_transient_cast[,names(cluster_transient_cast) != "ID"]
  
  # Run NMDS
  set.seed(1)
  nmds_transient <- metaMDS(cluster_transient_cast, distance = "bray", k = 2, trymax = 20)
  df_nmds_transient <- as.data.frame(nmds_transient$points)
  df_nmds_transient$ID <- row.names(df_nmds_transient)
  df_nmds_transient$Sample.name <- as.character(sapply(df_nmds_transient$ID, function(x) metadata$Sample.name[metadata$ID == x]))
  df_nmds_transient$Location <- sapply(df_nmds_transient$ID, function(x) metadata$Location[metadata$ID == x])
  
  # Plot NMDS
  transient_graphs[[i]] <- ggplot(df_nmds_transient, aes(MDS1, MDS2)) +
    geom_point(alpha = 0.5, colour = cols[unique(vir_counts_longus_demovir$sample_type)[i]]) +
    geom_encircle(aes(fill = Sample.name), alpha=0.3, expand = 0) +
    theme_classic() +
    guides(fill = FALSE) +
    scale_fill_manual(values = rep(cols[unique(vir_counts_longus_demovir$sample_type)[i]], length(unique(df_nmds_transient$Sample.name)))) +
    ggtitle(unique(vir_counts_longus_demovir$sample_type)[i])
  
  # Summarise number of clusters
  n_transient_clusters[i] <- length(unique(cluster_transient$vcontact_cluster))
  
  transient_sample <- data.frame(Sample.name = sapply(rownames(cluster_transient_cast), function(x) unique(cluster_transient$Sample.name[cluster_transient$ID == x]))) 
  permanova_transient[[i]] <- adonis(cluster_transient_cast ~ Sample.name, data = transient_sample)
}

# Number of individuals
n_individuals <- cluster_persist_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

# Number of persistent clusters
n_persist_clusters <- cluster_persist_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(vcontact_cluster))

# Number of transient clusters
n_transient_clusters <- cluster_transient_all %>% 
  group_by(sample_type) %>%
  summarise(n = n_distinct(vcontact_cluster))

tiff("figures/persistent_phage_clusters.tiff", width = 2000, height = 1000, res = 200)
grid.arrange(grobs = persist_graphs, layout_matrix = rbind(c(1,3), c(2,4)))
dev.off()

tiff("figures/transient_phage_clusters.tiff", width = 2000, height = 1000, res = 200)
grid.arrange(grobs = transient_graphs, layout_matrix = rbind(c(1,3), c(2,4)))
dev.off()

# Percentage of individuals where phage cluster is persistent, transient or doesn't exist
long_n <- metadata_longus %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

cluster_all <- rbind(cluster_persist_all, cluster_transient_all)

cluster_summary <- cluster_no_tp %>%
  filter(paste0(Sample.name, sample_type) %in% unique(paste0(cluster_all$Sample.name, cluster_all$sample_type))) %>%
  group_by(sample_type, vcontact_cluster) %>%
  summarise(persistent = n_distinct(Sample.name[no_tp >= 2]), transient = n_distinct(Sample.name[no_tp < 2])) %>%
  left_join(long_n, by = "sample_type") %>%
  mutate(perc_pers = persistent/n*100, perc_trans = transient/n*100) %>%
  mutate(perc_na = (n - (persistent + transient))/n*100)

cluster_summary_melt <- rbind(data.frame(sample_type = cluster_summary$sample_type, vcontact_cluster = cluster_summary$vcontact_cluster, 
                                   perc = cluster_summary$perc_pers, group = "Persistent"),
                              data.frame(sample_type = cluster_summary$sample_type, vcontact_cluster = cluster_summary$vcontact_cluster, 
                                         perc = cluster_summary$perc_trans, group = "Transient"),
                              data.frame(sample_type = cluster_summary$sample_type, vcontact_cluster = cluster_summary$vcontact_cluster, 
                                         perc = cluster_summary$perc_na, group = "Absent"))
cluster_summary_melt$group <- factor(cluster_summary_melt$group, levels = c("Absent", "Transient", "Persistent"))

# Plot persistency of phage clusters
cluster_summary_graphs <- list()
cluster_ranks <- data.frame()
for (i in 1:length(unique(cluster_summary$sample_type))) {
  cluster_sample_type <- cluster_summary_melt[cluster_summary_melt$sample_type == unique(cluster_summary$sample_type)[i],]
  cluster_sample_type$vcontact_cluster <- factor(cluster_sample_type$vcontact_cluster,
                                                 levels = cluster_sample_type$vcontact_cluster[cluster_sample_type$group == "Absent"][order(cluster_sample_type$perc[cluster_sample_type$group == "Absent"],
                                                                                                    cluster_sample_type$perc[cluster_sample_type$group == "Transient"])])
  
  cluster_summary_graphs[[i]] <- ggplot(cluster_sample_type, aes(vcontact_cluster, perc, fill = group)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_manual("Status", values = c("white", "lightblue", "blue"), labels = c("Absent", "Transient", "Persistent")) +
    ylab("% Individuals") + xlab("Phage clusters") +
    theme(axis.text.x = element_blank()) + 
    ggtitle(unique(cluster_summary$sample_type)[i])

}
tiff("figures/phage_pers_trans_perc.tiff", width = 3000, height = 1000, res = 200)
grid.arrange(grobs = cluster_summary_graphs, layout_matrix = rbind(c(1,2),c(3,4)))
dev.off()

# Taxonomy breakdown of persist phage clusters
cluster_persist_demovir <- cluster_persist_all %>%
  group_by(Sample.name, sample_type, demovir) %>%
  summarise(V1_sample = sum(V1)) %>%
  ungroup() %>% group_by(sample_type, Sample.name) %>%
  mutate(V1_prop = V1_sample/sum(V1_sample))

cluster_persist_n <- cluster_persist_demovir %>%
  group_by(sample_type, demovir) %>%
  summarise(n = n_distinct(Sample.name))

tiff("figures/barplot_demovir_persist_rel.tiff", width = 3000, height = 2000, res = 400)
ggplot(cluster_persist_demovir, aes(x = Sample.name, y = V1_prop, fill = demovir)) + theme_classic() +
  geom_bar(stat = "identity") +
  scale_fill_manual("Viral Family",  values = demovir_cols) +
  facet_wrap(~sample_type, scale = "free", shrink = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(face="italic")) +
  ylab("Proportion of mapped reads") + xlab("Sample")
dev.off()

# Taxonomy breakdown of transient phage clusters
cluster_transient_demovir <- cluster_transient_all %>%
  group_by(Sample.name, sample_type, demovir) %>%
  summarise(V1_sample = sum(V1)) %>%
  ungroup() %>% group_by(sample_type, Sample.name) %>%
  mutate(V1_prop = V1_sample/sum(V1_sample))

cluster_transient_n <- cluster_transient_demovir %>%
  group_by(sample_type, demovir) %>%
  summarise(n = n_distinct(Sample.name))

tiff("figures/barplot_demovir_transient_rel.tiff", width = 3000, height = 2000, res = 400)
ggplot(cluster_transient_demovir, aes(x = Sample.name, y = V1_prop, fill = demovir)) + theme_classic() +
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
metaphlan <- metaphlan[,names(metaphlan) %in% df_cluster_samples_nmds$ID]
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

# NMDS of metaphlan data
metaphlan_nmds <- metaMDS(t(metaphlan), distance = "bray", k = 2, trymax = 20)
df_metaphlan_nmds <- as.data.frame(metaphlan_nmds$points)
df_metaphlan_nmds$ID <- names(metaphlan)
names(df_metaphlan_nmds) = c("MDS1", "MDS2", "ID")
df_metaphlan_nmds <- left_join(df_metaphlan_nmds, metadata, by = "ID")

# Plot microbiome NMDS
tiff("figures/metaphlan_nmds.tiff", height= 1000, width = 1600, res = 200)
ggplot(df_metaphlan_nmds, aes(x = MDS1, y = MDS2, fill = sample_type)) +
  theme_classic() + geom_point(size = 1.5, alpha = 0.7, pch = 21) +
  scale_fill_manual("Body Site", values = cols,
                    guide = guide_legend(override.aes = list(shape = 21, size = 3))) +
  xlab("NMDS 1") + ylab("NMDS 2")
dev.off()

# Procrustes analysis and plot
protest_res <- protest(df_cluster_samples_nmds[,c(1:2)], df_metaphlan_nmds[,c(1:2)], scale = TRUE)
tiff("figures/procrustes.tiff", height = 1000, width = 1200, res = 150)
plot(protest_res, cex = 0.5, ar.col = "grey", )
points(protest_res, display = "target", col = "red", cex = 0.5)
dev.off()

# Percentage of crispr host predictions
crispr_host_percentage <- length(unique(vir_counts_prop_melt$Var1[!is.na(vir_counts_prop_melt$crispr_host)]))/length(unique(vir_counts_prop_melt$Var1))*100

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

########## Alpha-diversity and phage richness ####
# Join with metadata
vir_counts_prop_melt_meta <- left_join(vir_counts_prop_melt, metadata, by = c("Var2"="ID","sample_type","Location")) %>% rename(ID = Var2)
vir_counts_prop_melt_meta <- vir_counts_prop_melt_meta[vir_counts_prop_melt_meta$Visit_Number == 1,]

# Sample - phage cluster matrix
vir_cluster_counts <- dcast(vir_counts_prop_melt_meta[!is.na(vir_counts_prop_melt_meta$vcontact_cluster),], ID ~ vcontact_cluster, length)
rownames(vir_cluster_counts) <- vir_cluster_counts[,1]
vir_cluster_counts <- vir_cluster_counts[,-1]

# Remove samples with lower than quartile for
remove_ids <- rownames(vir_cluster_counts)[rowSums(vir_cluster_counts > 0) <= 3]
vir_cluster_counts <- vir_cluster_counts[!rownames(vir_cluster_counts) %in% remove_ids,]
metadata_richness <- metadata[!metadata$ID %in% unique(c(remove_ids, metadata$ID[!metadata$ID %in% vir_counts_prop_melt_meta$ID])),]

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
richness$n_phage_contigs <- sapply(richness$ID, function(x) length(which(vir_counts_prop_melt$Var2 == x)))

linear_mod <- lm(richness ~ n_phage_contigs, richness)
summary(linear_mod)
richness$predict <- predict(linear_mod, richness)
tiff("figures/richness_n_phage_contigs.tiff", width = 1500, height = 750, res = 150)
ggplot(richness) +
  geom_point(aes(n_phage_contigs, richness, color = Location_sampletype)) +
  xlab("No. phage contigs") + ylab("Phage Cluster Richness") +
  theme_classic() +
  scale_color_manual("Body Site - Geographical Location", values = cohort_cols) +
  geom_abline(color = "red", slope = linear_mod$coefficients[2], intercept = linear_mod$coefficients[1])
dev.off()

# For each group, remove samples with less than or equal to 100 phage contigs
unique_groups <- unique(richness_paired$group)
for (i in 1:length(unique_groups)) {
  remove_samples <- richness_paired$Sample.name[(richness_paired$no_phages <= 100 & richness_paired$group %in% unique_groups[i])]
  richness_paired <- richness_paired[!(richness_paired$Sample.name %in% remove_samples & richness_paired$group %in% unique_groups[i]),]
}

# Subsample matrix and calculate richness
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
#richness_graphs_ss[[length(richness_graphs_ss)+1]] <- g_legend(plotRichnessGraph(richness_paired_ss, richness_ttest_ss, cols))

# Plot graph
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
tiff("figures/alpha_diversity_subsampled.tiff", width = 2000, height = 2500, res = 250)
grid.arrange(grobs = richness_graphs_ss, layout_matrix = lay)
dev.off()


########## Size of phages ####
contig_data_meta <- left_join(vir_counts_prop_melt, metadata, by = c("Var2"="ID","sample_type","Location")) %>% rename(ID = Var2) %>%
  left_join(contig_data[,c("name", "circular", "size")], by = c("Var1"="name"))

countSizeCat <- function(contig_data_meta, size_lower, size_upper) {
  size_summary <- contig_data_meta %>%
    group_by(sample_type) %>%
    mutate(n_samples = n_distinct(ID)) %>% ungroup() %>%
    filter(size < size_upper, size >= size_lower) %>%
    group_by(sample_type, n_samples) %>%
    summarise(n = n_distinct(Var1)) %>%
    mutate(n_per_sample = n/n_samples) %>%
    mutate(size_lower = size_lower, category = paste(as.character(size_lower/1e3), "<=\n<", as.character(size_upper/1e3)))
  return(size_summary)
}

sizes_lower <- c(3e3, 1e4, 5e4, 2e5)
sizes_upper <- c(sizes_lower[-1], Inf)
phage_size_summary <- map2_df(.x = sizes_lower, .y = sizes_upper, .f = function(.x, .y) countSizeCat(contig_data_meta, .x, .y))
phage_size_summary$category <- gsub("\n< Inf", "", phage_size_summary$category)

tiff("figures/phage_size.tiff", width = 800, height = 500, res = 150)
ggplot(phage_size_summary, aes(reorder(category, size_lower), n_per_sample, fill = sample_type)) +
  geom_bar(stat = "identity") + xlab("Size (kb)") + ylab("No. unqiue phages") + 
  theme_classic() +
  scale_fill_manual("Body Site", values = cols)
dev.off()

# Size of all contigs
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

########## Jumbo circular phages ####
rownames(jumbophage_contigs) <- jumbophage_contigs$name
jumbophage_contigs$vcontact_cluster <- as.character(jumbophage_contigs$vcontact_cluster)
jumbophage_contigs$vcontact_cluster[jumbophage_contigs$vcontact_cluster %in% ""] <- NA
jumbophage_contigs$vcontact_cluster_Var1 <- jumbophage_contigs$vcontact_cluster
jumbophage_contigs$vcontact_cluster_Var1[is.na(jumbophage_contigs$vcontact_cluster_Var1)] <- jumbophage_contigs$name[is.na(jumbophage_contigs$vcontact_cluster_Var1)]

jumbophage_contigs_meta <- inner_join(jumbophage_contigs[,c("name", "circular", "size", "vcontact_cluster_Var1", "sample", "vcontact_cluster")], metadata, by = c("sample"="ID")) %>%
  left_join(vir_counts_prop_melt[,names(vir_counts_prop_melt) != "vcontact_cluster"], by = c("name"="Var1", "sample"="Var2","sample_type","Location")) %>% rename(ID = sample)
jumbophage_contigs_meta_circ <- jumbophage_contigs_meta[jumbophage_contigs_meta$circular,]

# Number of jumbo phages belonging to clusters
n_clusters_with_jumbophages <- length(unique(jumbophage_contigs_meta$vcontact_cluster[!is.na(jumbophage_contigs_meta$vcontact_cluster)]))
n_jumbophages_in_clusters <- length(jumbophage_contigs_meta$name[!is.na(jumbophage_contigs_meta$vcontact_cluster)])

# Percentage of cohorts containing a jumbophage
jumbophage_summary <- jumbophage_contigs_meta %>%
  filter(Visit_Number == 1) %>%
  group_by(Location, sample_type, circular) %>%
  summarise(n_samples_jumbophage = n_distinct(ID))

metadata_summary <- metadata %>%
  filter(Visit_Number == 1) %>%
  group_by(Location, sample_type) %>%
  summarise(n_samples = n_distinct(ID))

jumbophage_summary <- full_join(jumbophage_summary, metadata_summary, by = c("Location", "sample_type")) %>%
  mutate(perc_samples = 100*n_samples_jumbophage/n_samples, Location_sampletype = paste(sample_type, Location, sep = "\n"))

tiff("figures/jumbophage_percentage_samples.tiff", width = 1000, height = 500, res = 100)
ggplot(jumbophage_summary, aes(Location_sampletype, perc_samples, fill = circular)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Cohort") + ylab("% samples") + ylim(c(0,100))
dev.off()

# Persistent jumbophages
jumbophage_persistent <- jumbophage_contigs_meta %>% 
  filter(vcontact_cluster_Var1 %in% unique(as.character(cluster_persist_all$vcontact_cluster)) & Location == "US")

# Number of locations for each jumbophage
jumbophage_num_loc <- jumbophage_contigs_meta %>%
  group_by(vcontact_cluster_Var1) %>%
  summarise(num_locations = n_distinct(Location))

# Persistent and prevalent jumbophage clusters
jumbophage_prev <- jumbophage_num_loc %>% 
  filter(vcontact_cluster_Var1 %in% unique(as.character(cluster_persist_all$vcontact_cluster))) %>%
  filter(num_locations > 1)

jumbophage_prev_meta <- left_join(jumbophage_prev, jumbophage_contigs_meta, by = "vcontact_cluster_Var1")

# Remove phage duplicates
jumbophage_contigs_meta_nodup <- jumbophage_contigs_meta[!duplicated(paste(jumbophage_contigs_meta$name, jumbophage_contigs_meta$sample_type)),]

jumbophage_contigs_meta_nodup$sampletype_Location <- paste(jumbophage_contigs_meta_nodup$sample_type, "-", jumbophage_contigs_meta_nodup$Location)
h = 200
tiff("figures/jumbophages_size.tiff", width = 3000, height = 1000, res = 270)
set.seed(1)
ggplot(jumbophage_contigs_meta_nodup, aes(x = "Contigs", y = size/1000, colour = sample_type, shape = circular, size = circular)) +
  geom_jitter(alpha = 0.6) +
  geom_hline(aes(yintercept = h), colour = "red", linetype = "dotted") +
  theme_classic() + ylab("Contig size (kbp)") + xlab("") +
  scale_y_continuous(breaks = seq(200, 800, 50)) +
  scale_shape_manual(values = c(15, 19)) +
  scale_size_manual(values = c(1, 3)) +
  facet_grid(~Location, scale = "free", space = "free") +
  theme(panel.grid.major.y = element_line(colour = "grey90"), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)) +
  scale_colour_manual("Body Site", values = cols)
dev.off()


jumbophage_contigs_table <- jumbophage_contigs_meta %>%
  mutate(persistent = ifelse(vcontact_cluster %in% cluster_persist_all$vcontact_cluster & Location == "US", TRUE, FALSE)) %>%
  left_join(jumbophage_num_loc, by = "vcontact_cluster_Var1") %>%
  select(name, size, circular, vcontact_cluster, crispr_host, demovir, persistent, num_locations, Sample.name, Visit_Number, Location, sample_type) %>%
  rename(phage_family = demovir, location = Location, sample_name = Sample.name, visit_number = Visit_Number)

write.csv(jumbophage_contigs_table, "data/jumbophage_supplementary_table.csv", row.names = FALSE)

#### Jumbophage cluster network

# Extract sub-graph
vcontact_sig = fread("data/cc_sig1.0_mcl1.5.ntw")

# Only contigs connected to a jumbophage
vcontact_sig_jumbo <- vcontact_sig[V1 %in% jumbophage_contigs$name & V2 %in% jumbophage_contigs$name]

# Create graph
vcontact_igraph <- igraph::graph.data.frame(vcontact_sig_jumbo, directed = FALSE, vertices = jumbophage_contigs[unique(vcontact_sig_jumbo$V1),])

#E(vcontact_igraph)$weight = E(vcontact_igraph)$V3
#E(vcontact_igraph)$width = E(vcontact_igraph)$V3/10
vcontact_igraph <- simplify(vcontact_igraph)
#graph_layout = layout_with_kk(vcontact_igraph)
graph_layout <- layout_with_fr(vcontact_igraph)

# Taxonomy from co-clustering with RefSeq
vcontact = read.csv("data/genome_by_genome_overview.csv")
vcontact_refseq <- vcontact[!grepl("NODE", vcontact$Genome),]

clusters_tax = list()
for(i in unique(as.character(jumbophage_contigs$vcontact_cluster))) {
  tax_string = as.character(vcontact_refseq[which(vcontact_refseq$VC.Subcluster == i),"Genome"])
  if(length(tax_string[!is.na(tax_string)]) > 0) {
    clusters_tax[[i]] = tax_string[!is.na(tax_string)]
  }
}

labels_clusters = paste(as.character(sapply(clusters_tax[V(vcontact_igraph)$VC.Cluster], "[[", 1)), V(vcontact_igraph)$crispr_host)
labels_clusters = gsub("NULL", "", labels_clusters)
labels_clusters = gsub("NA", "", labels_clusters)
# labels_clusters[labels_clusters != " "] = paste(V(vcontact_igraph)$vcontact_cluster[labels_clusters != " "], labels_clusters[labels_clusters != " "])
# labels_clusters <- gsub("NA  ", "", labels_clusters)
labels_clusters <- gsub(" ", "", labels_clusters)

nclusters <- length(levels(as.factor(V(vcontact_igraph)$vcontact_cluster)))
rand_colours <- colorRampPalette(c(brewer.pal(8, "Spectral"),brewer.pal(8, "RdBu"), brewer.pal(8, "Dark2")))(nclusters)
set.seed(1)
vertex_cluster <- rand_colours[sample(1:nclusters, nclusters)]
vertex_cluster <- vertex_cluster[as.numeric(as.factor(V(vcontact_igraph)$vcontact_cluster))]
vertex_cluster[is.na(vertex_cluster)] <- "black"

rand_colours2 <- colorRampPalette(c(brewer.pal(8, "Spectral"),brewer.pal(8, "RdBu"), brewer.pal(8, "Dark2")))(length(unique(labels_clusters)))
set.seed(1)
vertex_host <- rand_colours2[sample(1:length(unique(labels_clusters)),length(unique(labels_clusters)))]
vertex_host <- vertex_host[as.numeric(as.factor(V(vcontact_igraph)$crispr_host))]
vertex_host[is.na(vertex_host)] <- "white"
names(vertex_host) <- as.factor(V(vcontact_igraph)$crispr_host)
names(vertex_host)[is.na(names(vertex_host))] <- "Unassigned"

tiff("figures/jumbophage_network.tiff", height = 2000, width = 2000, res = 200)
plot(vcontact_igraph,
     rescale = TRUE,
     layout = graph_layout,
     vertex.size = V(vcontact_igraph)$size/150000,
     vertex.frame.color = vertex_cluster,
     vertex.color = vertex_host,
     vertex.shape = c("square","circle")[as.numeric(as.factor(V(vcontact_igraph)$circular))],
     vertex.label = NA)
legend('topleft',legend=as.expression(lapply(unique(sort(names(vertex_host))), function(a) bquote(italic(.(a))))),col='black',pch=21, 
       cex = 0.8, pt.bg=unique(vertex_host[order(names(vertex_host))]), title = "Predicted Host Genus")
legend('topright',legend=c("circular", "linear"), col='black',pch=c(1,0), cex = 1)
dev.off()

######### Function of genes ####

# Summary of functional categories
jumbophage_gff <- read.csv("data/jumbophage_prodigal.gff", stringsAsFactors = FALSE)
jumbophage_gff_meta <- inner_join(jumbophage_gff, jumbophage_contigs_meta, by = c("V1" = "name"))

jumbophage_gff_summary <- jumbophage_gff_meta %>%
  filter(V1 %in% jumbophage_contigs_meta$name) %>%
  mutate(Name = ifelse(is.na(Name), "Hypothetical protein", Name)) %>%
  group_by(Name, Parent) %>%
  summarise(n = n_distinct(V1)) %>%
  ungroup() %>%
  mutate(Parent = ifelse(n >= 50, Parent, NA)) %>%
  mutate(Parent = ifelse(is.na(Parent), "Unknown", Parent)) %>%
  arrange(desc(n)) %>%
  rename("category" = "Parent")

jumbophage_gff_summary$Name <- factor(jumbophage_gff_summary$Name, levels = jumbophage_gff_summary$Name[order(jumbophage_gff_summary$n)])

fun_colours <- c("#E0A8E0", "#A6AAD1", "#D6A738", "#4D4CD4", "#6ED4C0", "#40A339", "#e88012","#B64ED9", "#d1d1d1")
names(fun_colours) <- sort(unique(jumbophage_gff_summary$category))

tiff("figures/jumbophage_functions.tiff", width = 2500, height = 2500, res = 200)
ggplot(jumbophage_gff_summary[jumbophage_gff_summary$n >= 50,], aes(Name, n, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  xlab("Protein function") + 
  scale_fill_manual("Functional Category", values = fun_colours, labels = names(fun_colours)) +
  scale_y_continuous("No. jumbo phages", breaks = seq(0, 500, 50))
dev.off()

jumbophage_fun_table <- jumbophage_gff_meta %>%
  select(V1, coding_region, protein, protein_description) %>%
  rename("name"="V1", "protein_function"="protein_description")

write.csv(jumbophage_fun_table, "data/jumbophage_function_supplementary_table.csv", row.names = FALSE)

# Group into metabolic profiles
jumbophage_fun_mat <- table(jumbophage_gff_meta$vcontact_cluster_Var1, jumbophage_gff_meta$protein)
jumbophage_fun_mat[jumbophage_fun_mat > 1] <- 1
jumbophage_fun_df <- as.data.frame.matrix(jumbophage_fun_mat)
jumbophage_cluster_crispr_host <- jumbophage_contigs_meta %>%
  group_by(vcontact_cluster_Var1, crispr_host) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
jumbophage_cluster_crispr_host <- jumbophage_cluster_crispr_host[!duplicated(jumbophage_cluster_crispr_host$vcontact_cluster_Var1),]

jumbophage_fun_df <- jumbophage_fun_df[,colSums(jumbophage_fun_df) > 0]
jumbophage_fun_df$crispr_host <- as.character(sapply(rownames(jumbophage_fun_df), function(x) 
  jumbophage_cluster_crispr_host$crispr_host[jumbophage_cluster_crispr_host$vcontact_cluster_Var1 == x]))
jumbophage_fun_df <- jumbophage_fun_df[!is.na(jumbophage_fun_df$crispr_host),]

# RF permutation test
jumbophage_fun_df$crispr_host <- as.factor(jumbophage_fun_df$crispr_host)
crisprhostRF <- rfPermute(crispr_host~., data = jumbophage_fun_df, ntree = 100, num.cores = 2, nrep = 50)
accuracyRF <- as.data.frame(confusionMatrix(crisprhostRF)[,colnames(confusionMatrix(crisprhostRF)) %in% c("pct.correct", "LCI_0.95", "UCI_0.95")])
accuracyRF$host <- rownames(accuracyRF)
accuracyRF$host <- factor(accuracyRF$host, level = accuracyRF$host)

# Importance
crisprhostRF_imp <- as.data.frame(rp.importance(crisprhostRF, scale = TRUE))
crisprhostRF_sig_imp <- crisprhostRF_imp[crisprhostRF_imp$MeanDecreaseAccuracy.pval < 0.05,]

crisprhostRF_sig_imp_gen <- crisprhostRF_sig_imp[,names(crisprhostRF_sig_imp) %in% unique(jumbophage_fun_df$crispr_host)]
crisprhostRF_sig_imp_pval <- crisprhostRF_sig_imp[,names(crisprhostRF_sig_imp) %in% paste0(unique(jumbophage_fun_df$crispr_host), ".pval")]

crisprhostRF_sig_imp_gen[crisprhostRF_sig_imp_pval >= 0.05] <- 0
crisprhostRF_sig_imp_gen <- crisprhostRF_sig_imp_gen[,colSums(crisprhostRF_sig_imp_gen) > 0]
crisprhostRF_sig_imp_gen <- crisprhostRF_sig_imp_gen[rowSums(crisprhostRF_sig_imp_gen) > 0,]
crisprhostRF_top10 <- as.matrix(crisprhostRF_sig_imp_gen)
crisprhostRF_top10[-order(crisprhostRF_sig_imp_gen, na.last=TRUE, decreasing=TRUE)[1:10]] <- NA
crisprhostRF_top10 <- apply(signif(crisprhostRF_top10, 3), 2, as.character)

rowLabels <- sapply(rownames(crisprhostRF_sig_imp_gen), function(x) unique(jumbophage_gff$protein_description[jumbophage_gff$protein %in% x]))
rowLabels <- paste(rowLabels, "(", rownames(crisprhostRF_sig_imp_gen), ")")
#rowColorLabels <- sapply(rownames(crisprhostRF_sig_imp_gen), function(x) unique(jumbophage_gff$Parent[jumbophage_gff$protein %in% x]))

# Heatmap of importance scores
tiff("figures/randomForest_importance.tiff", width = 1700, height = 2000, res = 100)
heatmap.2(as.matrix(crisprhostRF_sig_imp_gen),
          trace = "none",
          scale = "none",
          density.info = "none",
          hclustfun = function(x) {hclust(x, method = "ward.D2")},
          dendrogram = "both",
          col =  colorRampPalette(c("white", "darkblue"))(99),
          breaks = seq(min(crisprhostRF_sig_imp_gen), max(crisprhostRF_sig_imp_gen), length.out = 100),
          symbreaks = FALSE,
          keysize = 1,
          lhei = c(1,10),
          key.title = NA,
          key.xlab = "Importance score",
          key.ylab = NA,
          labRow = rowLabels,
          #RowSideColors = rowColors[rowColorLabels],
          labCol = names(crisprhostRF_sig_imp_gen),
          cexCol = 3,
          cexRow = 0.7,
          xlab = "",
          ylab = "",
          margin = c(30, 25),
          cellnote = crisprhostRF_top10
          )
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
jumbophages_trna <- read.table("data/jumbophage_trna_clean.tab", stringsAsFactors = FALSE)
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

# Write GFF file for largest three circular phages
circular_jumbophages <- unique(jumbophage_contigs_meta_circ$name[order(jumbophage_contigs_meta_circ$size, decreasing = TRUE)])[c(1,3)]
for (i in 1:length(circular_jumbophages)) {
  write.table(comb_gff[comb_gff$seqid == circular_jumbophages[i],], paste0("data/", circular_jumbophages[i], ".gff"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}

# Jumbo phages and args
args_data <- read.delim("data/jumbophage_proteins_card.out", stringsAsFactors = FALSE, header = FALSE)
args_data <- args_data[args_data$V3 >= 80,]