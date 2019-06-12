#!/bin/bash


mkdir -p examples/04_coli/input_bins

for i in $(grep kSGB_10068 examples/03_metagenomic/output_metagenomic.tsv | cut -f1); do
    cp -a examples/03_metagenomic/input_metagenomic/$i.fna examples/04_coli/input_bins/;
done;
