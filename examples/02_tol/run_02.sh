#!/bin/bash


#Get reference genomes for all species
phylophlan2_get_reference.py \
    -g all \
    -o input_genomes \
    -n 1 \
    --verbose 2>&1 | tee logs/phylophlan2_get_reference.log

#Run phylophlan
phylophlan2.py \
    -i input_genomes \
    -d phylophlan \
    -f tol_config.cfg \
    --diversity high \
    --fast \
    -o output_tol \
    --nproc 4 \
    --verbose 2>&1 | tee logs/phylophlan2.log

