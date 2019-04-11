#usr/bin/bash!

#Download the 135 S. aureus isolate genomes
./examples/01_saureus/download_saureus_genomes.sh

#Generate S. aureus database based on UniRef90
phylophlan_setup_database.py -g s__Staphylococcus_aureus -o examples/01_saureus/ \
--verbose 2>&1 | tee examples/01_saureus/logs/phylophlan_setup_database.log

#Build the phylogeny of the 135 S. aureus strains
phylophlan2.py \
-i examples/01_saureus/input_isolates \
-o output_isolates --output_folder examples/01_saureus/ \
-d s__Staphylococcus_aureus --databases_folder examples/01_saureus/ \
--trim greedy --not_variant_threshold 0.99 --remove_fragmentary_entries \
--fragmentary_threshold 0.67 --min_num_proteins 1064 --min_num_entries 135 \
-t a -f examples/01_saureus/isolates_config.cfg --diversity low --force_nucleotides --nproc 4 \
--verbose 2>&1 | tee examples/01_saureus/logs/phylophlan2__output_isolates.log

#Add S. aureus reference genomes
phylophlan_get_reference.py -g s__Staphylococcus_aureus -o examples/01_saureus/input_references/ \
-n 1000 --verbose 2>&1 | tee examples/01_saureus/logs/phylophlan_get_reference.log

cp -a examples/01_saureus/input_isolates/* examples/01_saureus/input_references/

phylophlan2.py \
-i examples/01_saureus/input_references \
-o output_references --output_folder examples/01_saureus/ \
-d s__Staphylococcus_aureus --databases_folder examples/01_saureus/ \
-t a -f examples/01_saureus/references_config.cfg --nproc 4 \
--subsample twentyfivepercent --diversity low --fast \
2>&1 |tee examples/01_saureus/logs/phylophlan2__reference_genomes__s__Staphylococcus_aureus.log

#Visualize the phylogenetic tree with GraPhlAn
./examples/01_saureus/graphlan/run_graphlan_isolates.sh
