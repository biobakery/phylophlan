#!/bin/bash


mkdir -p input_genomes;

# Find executables
[[ -z $(which phylophlan_get_reference) ]] && ppa_get_ref="phylophlan_get_reference.py" || ppa_get_ref=$(which phylophlan_get_reference)
[[ -z $(which phylophlan) ]] && ppa="phylophlan.py" || ppa=$(which phylophlan)

# Get reference genomes for all species
$ppa_get_ref -g all \
    -o input_genomes \
    -n 1 \
    --verbose 2>&1 | tee logs/phylophlan_get_reference.log

# Run phylophlan
$ppa -i input_genomes \
    -d phylophlan \
    -f tol_config.cfg \
    --diversity high \
    --fast \
    -o output_tol \
    --nproc 4 \
    --verbose 2>&1 | tee logs/phylophlan.log
