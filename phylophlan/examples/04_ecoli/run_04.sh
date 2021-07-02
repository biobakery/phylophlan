#!/bin/bash


# Find executables
[[ -z $(which phylophlan_setup_database) ]] && ppa_setup_db="phylophlan_setup_database.py" || ppa_setup_db=$(which phylophlan_setup_database)
[[ -z $(which phylophlan_get_reference) ]] && ppa_get_ref="phylophlan_get_reference.py" || ppa_get_ref=$(which phylophlan_get_reference)
[[ -z $(which phylophlan) ]] && ppa="phylophlan.py" || ppa=$(which phylophlan)

# Generate the phylogenetic markers database for E. coli based on the core set of UniRef90 proteins
$ppa_setup_db -o s__Escherichia_coli \
    -g s__Escherichia_coli \
    --verbose 2>&1 | tee logs/phylophlan_setup_database.log

# Retrieve 200 E. coli reference genomes
$ppa_get_ref -g s__Escherichia_coli \
    -o inputs/ \
    -n 200 \
    --verbose 2>&1 | tee logs/phylophlan_get_reference.log

# Retrieve 8 Ethiopian MAGs asssigned to E. coli
for i in $(grep kSGB_10068 ../03_metagenomic/output_metagenomic.tsv | cut -f1); do
    cp -a ../03_metagenomic/input_metagenomic/$i.fna inputs/;
done;

# Build the phylogeny
$ppa -i inputs \
    -o output_references \
    -d s__Escherichia_coli \
    -t a \
    -f references_config.cfg \
    --force_nucleotides \
    --nproc 4 \
    --subsample twentyfivepercent \
    --diversity low \
    --fast \
    --verbose 2>&1 |tee logs/phylophlan__inputs.log
