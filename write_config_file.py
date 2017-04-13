#!/usr/bin/env python3.6


import configparser as cp


config = cp.ConfigParser()
config_folder_files = 'configs/'

# AVAILABLE OPTIONS:
# program_name: executable name to use
# program_name_parallel: executable name of the parallel (or multi-core) version (if available)
# params: list of params to use
# input: specify input option after which the input file will be placed
# database:
# output: specify output option after which the output file will be placed
# version: (if available) used to verify that the software is installed and running


# command: "@program_name|@program_name_parallel @input @output @database @params"


progs = {
    'markers_db': {'program_name': 'usearch9.2.64_i86linux32', 'program_name_parallel': '', 'params': '-quiet', 'input': '-makeudb_ublast', 'output': '-output', 'version': '-version'},
    'dna_map': {'program_name': 'tblastn', 'program_name_parallel': '', 'params': '', 'input': '', 'database': '', 'output': '', 'version': '-version'},
    'aa_map': {'program_name': 'usearch9.2.64_i86linux32', 'program_name_parallel': '', 'params': '-quiet -threads 1 -evalue 1e-10 -maxaccepts 8 -maxrejects 32', 'input': '-ublast', 'database': '-db', 'output': '-blast6out', 'version': '-version'},
    'msa': {'program_name': 'muscle3.8.1551', 'program_name_parallel': '', 'params': '', 'input': '', 'output': '', 'version': '-version'},
    # {'program_name': 'mafft', 'program_name_parallel': '', 'params': '', 'input': '', 'output': '', 'version': '--version'},
    'tree': {'program_name': 'FastTree-2.1.9-SSE3', 'program_name_parallel': 'FastTreeMP-2.1.9-SSE3', 'params': '', 'input': '', 'output': '', 'version': ''},
    # {'program_name': 'raxmlHPC', 'program_name_parallel': 'raxmlHPC-PTHREADS-SSE3', 'params': '', 'input': '', 'output': '', 'version': '-v'}
}

for prog, options in progs.items():
    config[prog] = {}

    for option, value in options.items():
        config[prog][option] = value

with open(config_folder_files+'software.config', 'w') as f:
    config.write(f)
