#!/bin/bash

#Retrieve E.Coli genomes bins
./examples/04_coli/retrieve_coli_bins.sh

#Generate the phylogenetic markers database for E. coli based on the core set of UniRef90 proteins
phylophlan_setup_database.py -g s__Escherichia_coli -o examples/04_coli/ \
--verbose 2>&1 | tee examples/04_coli/logs/phylophlan_setup_database.log

# Add E. coli reference genomes
phylophlan_get_reference.py -g s__Escherichia_coli -o examples/04_coli/input_references/ \
-n 200 --verbose 2>&1 | tee examples/04_coli/logs/phylophlan_get_reference.log

cp -a examples/04_coli/input_bins/* examples/04_coli/input_references/

#Build the phylogeny
phylophlan2.py \
-i examples/04_coli/input_references \
-o output_references --output_folder examples/04_coli/ \
-d s__Escherichia_coli --databases_folder examples/04_coli/ \
-t a -f examples/04_coli/references_config.cfg --nproc 4 \
--subsample twentyfivepercent --diversity low --fast \
--verbose 2>&1 |tee examples/04_coli/logs/phylophlan2__s__Escherichia_coli.log
