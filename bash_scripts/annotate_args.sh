#!/bin/bash
/Users/vcarr/SOFTWARE/ncbi-blast-2.10.0+/bin/blastp -db ../db/CARD_3.0.0/protein_fasta_protein_homolog_model.fasta -query ../data/jumbophage_contigs_translations.faa -out ../data/jumbophage_proteins_card.out -outfmt 6 -num_threads 2 -evalue 1e-5
