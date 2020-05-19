#!/bin/sh

file=../data/alshayeb/gbk_files.txt
num=$(cat $file | wc -l)

for ((i=1;i<=$num;i++)); do
  gbk_file=$(sed -n "${i}p" $file)
  if grep -q 'Major capsid protein\|major capsid protein' ../data/alshayeb/phage-plasmid-unknown-genbank/${gbk_file}; then
    cp ../data/alshayeb/phage-plasmid-unknown-genbank/${gbk_file} ../data/alshayeb/major_capsid_protein_phages/
  fi
done
