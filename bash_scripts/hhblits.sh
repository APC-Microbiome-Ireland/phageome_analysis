#!/bin/sh

python3 ../helper_scripts/extract_proteins_from_largest_phages.py

ls -1 ../data/proteins/*faa > prot_list.txt
num_files=$(cat prot_list.txt | wc -l)

for ((i=1;i<=$num_files;i++)); do
	file=$(sed -n "${i}p" prot_list.txt)
	hhblits -i $file -blasttab ${file}.hhblit.tab -n 2 -e 0.00001 -d ../db/UNICLUST30/uniclust30_2018_08/uniclust30_2018_08
done

cat ../data/proteins/SRS016200_NODE_2_length_551943_cov_20.3899*.tab > ../data/proteins/SRS016200_NODE_2_length_551943_cov_20.3899.tab
cat ../data/proteins/SRS017007_NODE_1_length_547686_cov_137.51*.tab > ../data/proteins/SRS017007_NODE_1_length_547686_cov_137.51.tab
cat ../data/proteins/SRS047100_NODE_1_length_449263_cov_16.0097*.tab > ../data/proteins/SRS047100_NODE_1_length_449263_cov_16.0097.tab

awk '{ if ($11 <= 1E-05) { print } }' ../data/proteins/SRS016200_NODE_2_length_551943_cov_20.3899.tab > ../data/proteins/SRS016200_NODE_2_length_551943_cov_20.3899.1E-5.tab
awk '{ if ($11 <= 1E-05) { print } }' ../data/proteins/SRS017007_NODE_1_length_547686_cov_137.51.tab > ../data/proteins/SRS017007_NODE_1_length_547686_cov_137.51.1E-5.tab
awk '{ if ($11 <= 1E-05) { print } }' ../data/proteins/SRS047100_NODE_1_length_449263_cov_16.0097.tab > ../data/proteins/SRS047100_NODE_1_length_449263_cov_16.0097.1E-5.tab


