#!/usr/bin/bash


echo "GraPhlAn annotate"
graphlan_annotate.py \
    --annot graphlan/isolates_annotation.txt \
    output_isolates/RAxML_bestTree.input_isolates_refined.tre \
    graphlan/isolates_annotated.xml

echo "GraPhlAn draw"
graphlan.py \
    graphlan/isolates_annotated.xml \
    graphlan/saureus_isolates.png \
    --dpi 300
