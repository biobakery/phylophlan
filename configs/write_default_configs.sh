#!/bin/bash


# supermatrix.cfg
python3 ../write_config_file.py -o supermatrix.cfg --map_dna --map_aa --trim --tree2 --overwrite
# supertree.cfg
python3 ../write_config_file.py -o supertree.cfg --map_dna --map_aa --trim --gene_tree1 --gene_tree2 --overwrite
