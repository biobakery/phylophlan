#!/bin/bash

f=$(sort -k3 -n -r <(phylophlan_get_reference.py -l | grep p__Actinobacteria | grep s__)| head -2| cut -f7 -d'|'| cut -f1 -d'	')

one=$(echo $f | cut -f1 -d' ')
two=$(echo $f | cut -f2 -d' ')

phylophlan_get_reference.py -g $one -n 2 -o examples/05_chlamydiae/input_bins --verbose 2>&1 | tee examples/05_chlamydiae/logs/phylophlan_get_reference_${one}_actinobacteria.log
phylophlan_get_reference.py -g $two -n 2 -o examples/05_chlamydiae/input_bins --verbose 2>&1 | tee examples/05_chlamydiae/logs/phylophlan_get_reference_${two}_actinobacteria.log

