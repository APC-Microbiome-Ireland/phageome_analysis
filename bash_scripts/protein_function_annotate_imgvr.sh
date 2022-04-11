#!/bin/bash

# Predict proteins with Prodigal (v2.6.3)
prodigal -i ../data/imgvr_uvigs_scaffolds_200414.fasta -f gff -o ../data/imgvr_uvigs_prodigal.gff -a ../data/imgvr_uvigs_translations.faa -d ../data/imgvr_uvigs_genes.fna -p meta

# HMM search (v3.2.1) of pVOGs (downloaded 1st November 2019)
hmmscan -E 1e-5 --tblout ../data/imgvr_uvigs_pvogs_hmm_tblout.txt -o ../data/imgvr_uvigs_pvogs_out.txt --cpu 2 ../db/AllvogHMMprofiles/all_vogs.hmm ../data/imgvr_uvigs_translations.faa

# HMM search (v3.2.1) of pFAMs (downloaded 2nd September 2019)
hmmsearch --tblout ../data/imgvr_uvigs_pfam.out -E 1e-5 --cpu 2 ../db/PFAM/Pfam-A.hmm ../data/imgvr_uvigs_translations.faa > /dev/null

# HMM search (v3.2.1) of TIGRFAMs (downloaded 3rd September 2019)
hmmsearch --tblout ../data/imgvr_uvigs_tigrfams.out -E 1e-5 --cpu 2 ../db/TIGRFAMS/TIGRFAMS.hmm ../data/imgvr_uvigs_translations.faa > /dev/null
