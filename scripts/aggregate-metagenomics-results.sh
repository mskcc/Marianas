#! /bin/bash

# run it from the output folder

# create merged file
python ${mpa_dir}/utils/merge_metaphlan_tables.py *_profiled_metagenome.txt > merged-abundance-table.txt

# make heatmap
python ${mpa_dir}/utils/metaphlan_hclust_heatmap.py -c bbcry --top 50 --minv 0.1 -s log -d euclidean --in merged-abundance-table.txt --out abundance-heatmap.pdf








#
