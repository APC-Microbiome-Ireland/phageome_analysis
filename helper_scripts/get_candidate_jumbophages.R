library(data.table)

# Read contig catalogue
contig_data <- readRDS("../data/contig_data.RDS")

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

# Output jumbo phage headers
write.table(jumbophage_contigs$name, file = "../data/jumbophage_contigs_for_annot.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)