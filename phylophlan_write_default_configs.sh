#!/bin/bash


# supermatrix_nt.cfg
# python3 ../phylophlan_write_config_file.py -o supermatrix_nt.cfg \
phylophlan_write_config_file.py -o phylophlan_configs/supermatrix_nt.cfg \
    -d n \
    --db_dna makeblastdb \
    --map_dna blastn \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml \
    --overwrite
# supertree_nt.cfg
# python3 ../phylophlan_write_config_file.py -o supertree_nt.cfg \
phylophlan_write_config_file.py -o phylophlan_configs/supertree_nt.cfg \
    -d n \
    --db_dna makeblastdb \
    --map_dna blastn \
    --msa mafft \
    --trim trimal \
    --gene_tree1 fasttree \
    --gene_tree2 raxml \
    --tree1 astral \
    --overwrite
# supermatrix_aa.cfg
# python3 ../phylophlan_write_config_file.py -o supermatrix_aa.cfg \
phylophlan_write_config_file.py -o phylophlan_configs/supermatrix_aa.cfg \
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
# python3 ../phylophlan_write_config_file.py -o supertree_aa.cfg \
phylophlan_write_config_file.py -o phylophlan_configs/supertree_aa.cfg \
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
