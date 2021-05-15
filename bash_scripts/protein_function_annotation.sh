#!/bin/bash

# Predict proteins with Prodigal (v2.6.3)
prodigal -i ../data/jumbophage_contigs.fasta -f gff -o ../data/jumbophage_contigs_prodigal.gff -a ../data/jumbophage_contigs_translations.faa -d ../data/jumbophage_contigs_genes.fna -p meta

# HMM search (v3.2.1) of pVOGs (downloaded 1st November 2019)
hmmscan -E 0.00001 --tblout ../data/jumbophage_pvogs_hmm_tblout.txt -o ../data/jumbophage_pvogs_out.txt --cpu 2 ../db/AllvogHMMprofiles/all_vogs.hmm ../data/jumbophage_contigs_translations.faa

# HMM search (v3.2.1) of pFAMs (downloaded 2nd September 2019)
hmmsearch --tblout ../data/jumbophage_proteins_pfam.out -E 1e-5 --cpu 2 ../db/PFAM/Pfam-A.hmm ../data/jumbophage_contigs_translations.faa > /dev/null

# HMM search (v3.2.1) of TIGRFAMs (downloaded 3rd September 2019)
hmmsearch --tblout ../data/jumbophage_proteins_tigrfams.out -E 1e-5 --cpu 2 ../db/TIGRFAMS/TIGRFAMS.hmm ../data/jumbophage_contigs_translations.faa > /dev/null
