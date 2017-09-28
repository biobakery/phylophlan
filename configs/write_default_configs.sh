#!/bin/bash


# supermatrix_nt.cfg
python3 ../write_config_file.py \
    -o supermatrix_nt.cfg \
    -d n \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml \
    --overwrite
# supertree_nt.cfg
python3 ../write_config_file.py \
    -o supertree_nt.cfg \
    -d n \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --gene_tree1 fasttree \
    --gene_tree2 raxml \
    --tree1 astral \
    --overwrite
# supermatrix_aa.cfg
python3 ../write_config_file.py \
    -o supermatrix_aa.cfg \
    -d a \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml \
    --overwrite
# supertree_aa.cfg
python3 ../write_config_file.py \
    -o supertree_aa.cfg \
    -d a \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --gene_tree1 fasttree \
    --gene_tree2 raxml \
    --tree1 astral \
    --overwrite
# ppa_run_1.cfg
# python3 ../write_config_file.py -o ppa_run_1.cfg --map_dna --map_aa --msa muscle --tree1 fasttree --overwrite
# ppa_run_2.cfg
# python3 ../write_config_file.py -o ppa_run_2.cfg --map_dna --map_aa --msa muscle --tree1 fasttree --overwrite
# ppa_run_4.cfg
# python3 ../write_config_file.py -o ppa_run_4.cfg --map_dna --map_aa --msa mafft --tree1 fasttree --overwrite
# ppa_run_5.cfg
# python3 ../write_config_file.py -o ppa_run_5.cfg --map_dna --map_aa --msa mafft --tree1 fasttree --overwrite
# ppa_run_7.cfg
# python3 ../write_config_file.py -o ppa_run_7.cfg --map_dna --map_aa --msa muscle --tree1 fasttree --tree2 raxml --overwrite
# ppa_run_8.cfg
# python3 ../write_config_file.py -o ppa_run_8.cfg --map_dna --map_aa --msa muscle --tree1 fasttree --tree2 raxml --overwrite
# ppa_run_10.cfg
# python3 ../write_config_file.py -o ppa_run_10.cfg --map_dna --map_aa --msa mafft --tree1 fasttree --tree2 raxml --overwrite
# ppa_run_11.cfg
# python3 ../write_config_file.py -o ppa_run_11.cfg --map_dna --map_aa --msa mafft --tree1 fasttree --tree2 raxml --overwrite
# ppa_run_13.cfg
# python3 ../write_config_file.py -o ppa_run_13.cfg --map_dna --map_aa --msa mafft --gene_tree1 fasttree --tree1 astral --overwrite
# ppa_run_14.cfg
# python3 ../write_config_file.py -o ppa_run_14.cfg --map_dna --map_aa --msa mafft --gene_tree1 fasttree --tree1 astral --overwrite
# ppa_run_16.cfg
# python3 ../write_config_file.py -o ppa_run_16.cfg --map_dna --map_aa --msa mafft --gene_tree1 fasttree --tree1 astral --gene_tree2 raxml --overwrite
# ppa_run_17.cfg
# python3 ../write_config_file.py -o ppa_run_17.cfg --map_dna --map_aa --msa mafft --gene_tree1 fasttree --tree1 astral --gene_tree2 raxml --overwrite
