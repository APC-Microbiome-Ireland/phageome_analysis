#!/bin/bash

# Identify trna from megaphages using aragorn (v1.2.36)
aragorn -t -i -c -d -w -o ../data/megaphage_trna.tab ../data/megaphage_contigs.fasta

# Clean output
python3 ../helper_scripts/clean_aragorn_output.py
