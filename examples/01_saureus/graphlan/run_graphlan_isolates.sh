#/usr/bin/bash!

graphlan_annotate.py --annot examples/01_saureus/graphlan/isolates_annotation.txt \
examples/01_saureus/output_isolates/RAxML_bestTree.input_isolates_refined.tre \
examples/01_saureus/graphlan/isolates_annotated.xml

graphlan.py --dpi 300 \
examples/01_saureus/graphlan/isolates_annotated.xml \
examples/01_saureus/graphlan/saureus_isolates.png
