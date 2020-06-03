#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $SCRIPT_DIR || exit 1

echo "# Downloading 10 S. aureus isolate genomes"
rm -rf input_isolates
mkdir -p input_isolates

for i in $(cut -f4 135_saureus_isolates.tsv | sed 1d | head -n 11); do
    o=`basename $i | cut -f1 -d'.'`;
    wget $i -O input_isolates/${o}.fna.gz;
done;

echo "# Generating S. aureus database based on UniRef90"
phylophlan_setup_database \
    -g s__Staphylococcus_aureus \
    --max_proteins 100 \
    --overwrite \
    --verbose || exit 1

echo "# Writing default config"
CONFIG_DIR='configs'
phylophlan_write_default_configs.sh $CONFIG_DIR

echo "# Writing isolates config file"
phylophlan_write_config_file -o isolates_config.cfg \
    --overwrite \
    -d a \
    --force_nucleotides \
    --db_aa diamond \
    --map_aa diamond \
    --map_dna diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml || exit 1

echo "# Building the phylogeny of the 10 S. aureus strains"
phylophlan \
    -i input_isolates \
    -o output_isolates \
    -d s__Staphylococcus_aureus \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    --min_num_proteins 1 \
    --min_num_entries 4 \
    --min_num_markers 1 \
    -t a \
    --configs_folder $CONFIG_DIR \
    -f isolates_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc 2 \
    --verbose || exit 1

echo "# Adding 5 S. aureus reference genomes"
phylophlan_get_reference \
    -g s__Staphylococcus_aureus \
    -o input_references \
    -n 5 \
    --verbose || exit 1

cp -a input_isolates/* input_references/ || exit 1

echo "# Building the phylogeny of the 15 S. aureus genomes"
phylophlan \
    -i input_references \
    -o output_references \
    -d s__Staphylococcus_aureus \
    -t a \
    --min_num_proteins 1 \
    --min_num_entries 4 \
    --min_num_markers 1 \
    -f references_config.cfg \
    --nproc 2 \
    --subsample twentyfivepercent \
    --diversity low \
    --fast \
    --verbose || exit 1

