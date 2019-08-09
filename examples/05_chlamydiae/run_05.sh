#!/bin/bash


# Download all MAGs of the uSGB 19436
./download_files.sh sequence_url_opendata_19436.txt input_bins/

# Download all reference genomes for all species in the Chlamydiae phylum
phylophlan2_get_reference.py \
    -g p__Chlamydiae \
    -n -1 \
    -o input_bins \
    --verbose 2>&1 | tee logs/phylophlan2_get_reference__chlamydiae.log

# Download all reference genomes for all species in the Planctomycetes phylum (outgroup for rooting)
phylophlan2_get_reference.py \
    -g p__Planctomycetes \
    -n -1 \
    -o input_bins \
    --verbose 2>&1 | tee logs/phylophlan2_get_reference__planctomycetes.log

# Retrieve 10 Ethiopian MAGs assigned to the uSGB 19436 (Chalmydiae phylum)
for i in $(grep uSGB_19436 ../03_metagenomic/output_metagenomic.tsv | cut -f1); do
    cp -a ../03_metagenomic/input_metagenomic/$i.fna input_bins/;
done;

# Build the phylogeny
phylophlan2.py \
    -i input_bins \
    -d phylophlan \
    --diversity high \
    --accurate \
    -f chlamydiae_config.cfg \
    -o output_chlamydiae \
    --force_nucleotides \
    --nproc 4 \
    -t a \
    --verbose 2>&1 | tee logs/phylophlan2__input_bins.log 


