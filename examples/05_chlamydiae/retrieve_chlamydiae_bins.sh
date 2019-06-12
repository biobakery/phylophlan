#!/bin/bash
mkdir -p examples/05_chlamydiae/input_bins

for i in $(grep uSGB_19436 examples/03_metagenomic/output_metagenomic.tsv | cut -f1); do
    cp -a examples/03_metagenomic/input_metagenomic/$i.fna examples/05_chlamydiae/input_bins/;
done;
