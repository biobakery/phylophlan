#usr/bin/bash!

mkdir -p examples/04_coli/input_bins

for i in $(cat examples/03_metagenomic/output_metagenomic.tsv | grep Escherichia_coli | cut -f1); do
cp -a examples/03_metagenomic/input_metagenomic/$i.fna  examples/04_coli/input_bins
done;
