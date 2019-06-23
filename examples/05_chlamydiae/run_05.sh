# !/bin/bash

#Retrieve the genome bins in uSGB 19436 found in 03.Metagenomic
./examples/05_chlamydiae/retrieve_chlamydiae_bins.sh 

#Download other genome bins in uSGB 19436 
./examples/05_chlamydiae/download_files.sh
examples/05_chlamydiae/sequence_url_opendata_19436.txt \
examples/05_chlamydiae/input_bins

#Download other genome bins in uSGB 19435
./examples/05_chlamydiae/download_files.sh \
examples/05_chlamydiae/sequence_url_opendata_19435.txt \
examples/05_chlamydiae/input_bins

#DOwnload reference genomes for each species in Chlamydiae phylum
phylophlan_get_reference.py -g p__Chlamydiae -n 2 \
-o examples/05_chlamydiae/input_bins --verbose 2>&1 | \
tee examples/05_chlamyidiae/logs/phylophlan_get_reference.log

#Download genomes from Actinobacteria as outgroup
./examples/05_chlamydiae/retrieve_actinobacteria_genomes.sh

#Build the phylogeny
phylophlan2.py -i examples/05_chlamydiae/input_bins \
-d phylophlan --diversity medium --accurate -f examples/05_chlamydiae/chlamydiae_config.cfg \
-o output_chlamydiae --output_folder examples/05_chlamydiae/ --nproc 4 -t a \
--verbose 2>&1 | tee examples/05_chlamydiae/logs/phylophlan2_chlamydiae.log 


