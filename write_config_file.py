#!/usr/bin/env python3


import configparser as cp


config = cp.ConfigParser()
config_folder_files = 'configs/'

# AVAILABLE OPTIONS:
# program_name: name of the executable to use
# program_name_parallel: name of the parallel (or multi-core) version of the executable
# params: params to use
# threads: specify the option to use to pass the number of threads
# input: specify the option to use for the input file
# database: specify the option to use for setting the database
# output_path: specify the to option to use to set the path of the folder that will contains the output file
# output: specify the option to use for the output file
# version: specify the option to use to get the version of the sotware, used to verify the software installation
# command_line: specify the command line to generate with the position of each argument

progs = {
    'markers_db': {'program_name': 'usearch9.2.64_i86linux32',
                   'params': '-quiet',
                   'input': '-makeudb_ublast',
                   'output': '-output',
                   'version': '-version',
                   'command_line': '#program_name# #params# #input# #output#'},
    'dna_map': {'program_name': 'tblastn',
                'params': '-outfmt 6 -evalue 1e-50',
                'input': '-subject',
                'database': '-query',
                'output': '-out',
                'version': '-version',
                'command_line': '#program_name# #params# #input# #output#'},
    'aa_map': {'program_name': 'usearch9.2.64_i86linux32',
               'params': '-quiet -evalue 1e-10 -maxaccepts 8 -maxrejects 32',
               'threads': '-threads',
               'input': '-ublast',
               'database': '-db',
               'output': '-blast6out',
               'version': '-version',
               'command_line': '#program_name# #params# #threads# #input# #database# #output#'},
    'msa': {'program_name': 'muscle3.8.1551',
            'params': '-quiet -maxiters 2',
            'input': '-in',
            'output': '-out',
            'version': '-version',
            'command_line': '#program_name# #params# #input# #output#'},
           # {'program_name': 'mafft', 'params': '--anysymbol --quiet', 'version': '--version'},
    'trim': {'program_name': 'trimal',
             'params': '-gappyout',
             'input': '-in',
             'output': '-out',
             'version': '--version',
             'command_line': '#program_name# #params# #input# #output#'},
    'tree': {'program_name': 'FastTree-2.1.9-SSE3',
             'program_name_parallel': 'FastTreeMP-2.1.9-SSE3',
             'params': '-quiet -fastest -mlnni 4 -no2nd',
             'output': '-out',
             'command_line': '#program_name# #params# #output# #input#'},
            # {'program_name': 'raxmlHPC', 'program_name_parallel': 'raxmlHPC-PTHREADS-SSE3', 'params': '-m PROTCATWAG -p 1989', 'threads': '-T', 'input': '-s', 'output_path':'-w', 'output': '-n', 'version': '-v'}
}

for prog, options in progs.items():
    config[prog] = {}

    for option, value in options.items():
        config[prog][option] = value

with open(config_folder_files+'software.config', 'w') as f:
    config.write(f)
