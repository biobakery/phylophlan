#!/usr/bin/bash

#Execute phylophlan_metagenomic
phylophlan_metagenomic.py -i examples/03_metagenomic/input_metagenomic \
-o examples/03_metagenomic/output_metagenomic --nproc 4 -n 1 --verbose \
2>&1 | tee examples/03_metagenomic/logs/phylophlan_metagenomic.log

#Draw heatmaps
phylophlan_draw_heatmaps.py -i examples/03_metagenomic/output_metagenomic.tsv \
 --map examples/03_metagenomic/bin2meta.tsv --top 20 --verbose \
2>&1 | tee examples/03_metagenomic/logs/phylophlan_draw_heatmaps.log
