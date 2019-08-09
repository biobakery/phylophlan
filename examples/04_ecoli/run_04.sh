#!/bin/bash


# Generate the phylogenetic markers database for E. coli based on the core set of UniRef90 proteins
phylophlan2_setup_database.py \
    -g s__Escherichia_coli \
    --verbose 2>&1 | tee logs/phylophlan2_setup_database.log

# Retrieve 200 E. coli reference genomes
phylophlan2_get_reference.py \
    -g s__Escherichia_coli \
    -o inputs/ \
    -n 200 \
    --verbose 2>&1 | tee logs/phylophlan2_get_reference.log

# Retrieve 8 Ethiopian MAGs asssigned to E. coli
for i in $(grep kSGB_10068 ../03_metagenomic/output_metagenomic.tsv | cut -f1); do
    cp -a ../03_metagenomic/input_metagenomic/$i.fna inputs/;
done;

# Build the phylogeny
phylophlan2.py \
    -i inputs \
    -o output_references \
    -d s__Escherichia_coli \
    -t a \
    -f references_config.cfg \
    --force_nucleotides \
    --nproc 4 \
    --subsample twentyfivepercent \
    --diversity low \
    --fast \
    --verbose 2>&1 |tee logs/phylophlan2__inputs.log
