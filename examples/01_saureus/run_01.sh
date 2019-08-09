#!/bin/bash


# Download the 135 S. aureus isolate genomes
mkdir -p input_isolates;

for i in $(cut -f4 135_saureus_isolates.tsv | sed 1d); do
    o=`basename $i | cut -f1 -d'.'`;
    wget $i -O input_isolates/${o}.fna.gz;
done;

# Generate S. aureus database based on UniRef90
phylophlan2_setup_database.py \
    -g s__Staphylococcus_aureus \
    --verbose 2>&1 | tee logs/phylophlan2_setup_database.log

# Build the phylogeny of the 135 S. aureus strains
phylophlan2.py \
    -i input_isolates \
    -o output_isolates \
    -d s__Staphylococcus_aureus \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    --min_num_proteins 1064 \
    --min_num_entries 135 \
    -t a \
    -f isolates_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc 4 \
    --verbose 2>&1 | tee logs/phylophlan2__output_isolates.log

# Add 1,000 S. aureus reference genomes
phylophlan2_get_reference.py \
    -g s__Staphylococcus_aureus \
    -o input_references \
    -n 1000 \
    --verbose 2>&1 | tee logs/phylophlan2_get_reference.log

cp -a input_isolates/* input_references/

# Build the phylogeny of the 1,135 S. aureus genomes
phylophlan2.py \
    -i input_references \
    -o output_references \
    -d s__Staphylococcus_aureus \
    -t a \
    -f references_config.cfg \
    --nproc 4 \
    --subsample twentyfivepercent \
    --diversity low \
    --fast \
    --verbose 2>&1 |tee logs/phylophlan2__output_references.log

# Visualize the phylogenetic tree with GraPhlAn
# GraPhlAn is Python2-based and have different requirements than PhyloPhlAn2
graphlan/run_graphlan_isolates.sh
