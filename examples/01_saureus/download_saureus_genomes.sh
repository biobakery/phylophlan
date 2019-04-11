#!/bin/bash

mkdir examples/01_saureus/input_isolates

for i in $(cut -f4 examples/01_saureus/135_saureus_isolates.tsv | sed 1d); do
    o=`basename $i | cut -f1 -d'.'`; 
    wget $i -O examples/01_saureus/input_isolates/${o}.fna.gz;
done
