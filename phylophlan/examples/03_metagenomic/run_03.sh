#!/usr/bin/bash


# Retrieve 369 Ethiopian MAGs
wget https://www.dropbox.com/s/fuafzwj67tguj31/ethiopian_mags.tar.bz2?dl=1 -O ethiopian_mags.tar.bz2;
mkdir -p input_metagenomic;
tar -xjf ethiopian_mags.tar.bz2 -C input_metagenomic/;

# Find executables
[[ -z $(which phylophlan_metagenomic) ]] && ppa_meta="phylophlan_metagenomic.py" || ppa_meta=$(which phylophlan_metagenomic)
[[ -z $(which phylophlan_draw_metagenomic) ]] && ppa_draw_meta="phylophlan_draw_metagenomic.py" || ppa_draw_meta=$(which phylophlan_draw_metagenomic)

# Execute phylophlan_metagenomic
# the parameter "-d SGB.Jan19" specify the first database released by MetaRefSGB and is just an example. The user can use the command: 
# "phylophlan_metagenomic.py --database_list" to see the list of available databases
$ppa_meta -i input_metagenomic \
    -d SGB.Jan19 \
    -o output_metagenomic \
    --nproc 4 \
    -n 1 \
    --verbose 2>&1 | tee logs/phylophlan_metagenomic.log

# Draw heatmaps
$ppa_draw_meta -i output_metagenomic.tsv \
    --map bin2meta.tsv \
    -o output_heatmap \
    --top 20 \
    --verbose 2>&1 | tee logs/phylophlan_draw_metagenomic.log
