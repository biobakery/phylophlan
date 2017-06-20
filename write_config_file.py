#!/usr/bin/env python3


import os
import sys
import argparse as ap
import configparser as cp


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = ap.ArgumentParser(description="")

    p.add_argument('-o', '--output', type=str, required=True, help="Specify the output file where to write the configurations. If the file already exists the program won't overwrite the alredy existing file")
    p.add_argument('--map_dna', action='store_true', default=False, help="Specify to add the map_dna section. Note, at least one of map_dna and map_aa must be specified")
    p.add_argument('--map_aa', action='store_true', default=False, help="Specify to add the map_aa section. Note, at least one of map_dna and map_aa must be specified")
    p.add_argument('--trim', action='store_true', default=False, help="Specify to add the trim section")
    p.add_argument('--gene_tree1', action='store_true', default=False, help="Specify to add the gene_tree1 section")
    p.add_argument('--gene_tree2', action='store_true', default=False, help="Specify to add the gene_tree2 section")
    p.add_argument('--tree2', action='store_true', default=False, help="Specify to add the tree2 section")

    return p.parse_args()


def check_params(args):
    if (not args.map_dna) and (not args.map_aa):
        error('at least one of --map_dna and --map_aa must be specified', exit=True)

    if os.path.isfile(args.output):
        error('cannot write ouptut file {} because it already exists'.format(args.output), exit=True)


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


if __name__ == '__main__':
    args = read_params()
    check_params(args)
    config = cp.ConfigParser()

    progs = {
        'db_aa': {'program_name': 'usearch9.2.64_i86linux32',
                  'params': '-quiet',
                  'input': '-makeudb_ublast',
                  'output': '-output',
                  'version': '-version',
                  'command_line': '#program_name# #params# #input# #output#'},
        'msa': {'program_name': 'muscle3.8.1551',
                'params': '-quiet -maxiters 2',
                'input': '-in',
                'output': '-out',
                'version': '-version',
                'command_line': '#program_name# #params# #input# #output#'}
               # {'program_name': 'mafft',
               #  'params': '--anysymbol --quiet',
               #  'input': '',
               #  'output': '',
               #  'version': '--version'},
               #  'command_line': '#program_name# #params# #input# #output#'}
    }

    if args.map_dna:
        progs['map_dna'] = {'program_name': 'tblastn',
                            'params': '-outfmt "6 saccver qaccver pident length mismatch gapopen sstart send qstart qend evalue bitscore" -evalue 1e-50',
                            'input': '-subject',
                            'database': '-query',
                            'output': '-out',
                            'version': '-version',
                            'command_line': '#program_name# #params# #input# #database# #output#'}

    if args.map_aa:
        progs['map_aa'] = {'program_name': 'usearch9.2.64_i86linux32',
                           'params': '-quiet -evalue 1e-10 -maxaccepts 8 -maxrejects 32',
                           'threads': '-threads',
                           'input': '-ublast',
                           'database': '-db',
                           'output': '-blast6out',
                           'version': '-version',
                           'command_line': '#program_name# #params# #threads# #input# #database# #output#'}

    if args.trim:
        progs['trim'] = {'program_name': 'trimal',
                         'params': '-gappyout',
                         'input': '-in',
                         'output': '-out',
                         'version': '--version',
                         'command_line': '#program_name# #params# #input# #output#'}

    if args.gene_tree1:
        progs['tree1'] = {'program_name': 'java -jar /CM/tools/astral-4.11.1/astral.4.11.1.jar',
                          'input': '-i',
                          'output': '-o',
                          'version': '--help',
                          'command_line': '#program_name# #input# #output#'}
        progs['gene_tree1'] = {'program_name': 'FastTree-2.1.9-SSE3',
                               'params': '-quiet -mlacc 2 -slownni -spr 4 -fastest -mlnni 4 -no2nd',
                               'output': '-out',
                               'command_line': '#program_name# #params# #output# #input#'}
    else:
        progs['tree1'] = {'program_name_parallel': 'FastTreeMP-2.1.9-SSE3',
                          'params': '-quiet -mlacc 2 -slownni -spr 4 -fastest -mlnni 4 -no2nd',
                          'output': '-out',
                          'command_line': '#program_name_parallel# #params# #output# #input#'}

    if args.gene_tree2:
        progs['gene_tree2'] = {'program_name': 'raxmlHPC',
                               'params': '-p 1989',
                               'model': '-m',
                               'database': '-g', # starting tree
                               'input': '-s',
                               'output_path':'-w',
                               'output': '-n',
                               'version': '-v',
                               'command_line': '#program_name# #model# #params# #database# #output_path# #input# #output#'}

    if args.tree2:
        progs['tree2'] = {'program_name_parallel': 'raxmlHPC-PTHREADS-SSE3',
                          'params': '-m PROTGAMMAAUTO -p 1989',
                          'threads': '-T',
                          'database': '-g', # starting tree
                          'input': '-s',
                          'output_path':'-w',
                          'output': '-n',
                          'version': '-v',
                          'command_line': '#program_name_parallel# #params# #threads# #database# #output_path# #input# #output#'}

    for prog, options in progs.items():
        config[prog] = {}

        for option, value in options.items():
            config[prog][option] = value

    with open(args.output, 'w') as f:
        config.write(f)

    sys.exit(0)
