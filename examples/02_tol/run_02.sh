#!usr/bin/bash

#Get reference genomes for all species
phylophlan_get_reference.py -g all -o examples/02_tol/input_genomes/ \
-n 1 --verbose 2>&1 | tee examples/02_tol/logs/phylophlan_get_reference.log

#Run phylophlan
phylophlan2.py -i examples/02_tol/input_genomes \
-d phylophlan -f 02_tol.cfg --diversity high --fast \
-o output_tol --output_folder output_genomes \
--verbose 2>&1 | tee examples/02_tol/logs/phylophlan2.log

