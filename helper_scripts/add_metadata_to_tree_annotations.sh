#!/bin/bash

cd ../db/itol

cp binary_template.txt circular.txt
cat ../../data/jumbophage_mcp_metadata_edit_cols.txt | sed '1d' | awk -F'\t' '$3=="circular" {print $1"\t"1}' >> circular.txt
cat ../../data/jumbophage_mcp_metadata_edit_cols.txt | sed '1d' | awk -F'\t' '$3=="linear" {print $1"\t""-1"}' >> circular.txt

cp simplebar_template.txt sizes_imgvr.txt
grep "IMG/VR" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\t"$2}' >> sizes_imgvr.txt

cp simplebar_template.txt sizes_devoto.txt
grep "Devoto" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\t"$2}' >> sizes_devoto.txt

cp simplebar_template.txt sizes_alshayeb.txt
grep "Al-Shayeb" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\t"$2}' >> sizes_alshayeb.txt

cp simplebar_template.txt sizes_this_study.txt
grep "This study" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\t"$2}' >> sizes_this_study.txt

cp simplebar_template.txt sizes_refseq.txt
grep "RefSeq" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\t"$2}' >> sizes_refseq.txt

cp colorstrip_template.txt environment.txt
cat environment_legend.txt >> environment.txt
echo "DATA" >> environment.txt
cat ../../data/jumbophage_mcp_metadata_edit_cols.txt | sed '1d' | awk -F'\t' '{print $1"\t"$9}' >> environment.txt

cp colorstrip_template.txt hosts.txt
cat host_legend.txt >> hosts.txt
echo "DATA" >> hosts.txt
cat ../../data/jumbophage_mcp_metadata_edit_cols.txt | sed '1d' | awk -F'\t' '{print $1"\t"$10}' >> hosts.txt

cp tree_colors_template.txt clusters.txt
grep "VC_" ../../data/jumbophage_mcp_metadata_edit_cols.txt | awk -F'\t' '{print $1"\trange\t"$13"\t"$12}' >> clusters.txt

cp labels_template.txt cluster_labels.txt
cat ../../data/jumbophage_mcp_metadata_edit_cols.txt | sed '1d' | awk -F'\t' '{print $1"\t"$12}' >> cluster_labels.txt
