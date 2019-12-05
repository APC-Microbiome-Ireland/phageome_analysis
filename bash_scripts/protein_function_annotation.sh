#!/bin/bash

# Predict proteins with Prodigal (v2.6.3)
prodigal -i ../data/megaphage_contigs.fasta -f gff -o ../data/megaphage_contigs_prodigal.gff -a ../data/megaphage_contigs_translations.faa -d ../data/megaphage_contigs_genes.fna -p meta

# remove trailing white space
sed 's/\ .*//g' megaphage_contigs_genes.fna > megaphage_contigs_genes2.fna

# DIAMOND (v0.38.5) against RefSeq (downloaded 3rd September 2019)
diamond makedb --in nr.gz -d nr
diamond blastp -q ../data/megaphage_contigs_translations.faa -d ../db/REFSEQ_NR_PROT/nr -o ../data/megaphage_proteins_diamond.out --outfmt 6 --threads 2 --evalue 1e-5

# HMM search (v3.2.1) of pVOGs (downloaded 1st November 2019)
hmmscan -E 0.00001 --tblout ../data/megaphage_pvogs_hmm_tblout.txt -o ../data/megaphage_pvogs_out.txt --cpu 2 ../db/AllvogHMMprofiles/all_vogs.hmm ../data/megaphage_contigs_translations.faa

# Select first description from pVOG metadata (downloaded 8th November 2019)
awk -F'[' '{print $1}' ../db/pVOGs/VOGTable.txt > ../db/pVOGs/VOGTable_singleprot.txt

# HMM search (v3.2.1) of pFAMs (downloaded 2nd September 2019)
hmmsearch --tblout megaphage_proteins_pfam.out -E 1e-5 --cpu 2 ../db/PFAM/Pfam-A.hmm ../data/megaphage_contigs_translations.faa > /dev/null

# HMM search (v3.2.1) of TIGRFAMs (downloaded 3rd September 2019)
hmmsearch --tblout megaphage_proteins_tigrfams.out -E 1e-5 --cpu $NSLOTS ../db/TIGRFAMS/TIGRFAMS.hmm ../data/megaphage_contigs_translations.faa > /dev/null

