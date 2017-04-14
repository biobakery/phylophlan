#!/usr/bin/env python3.6


__author__ = 'Nicola Segata (nsegata@hsph.harvard.edu), Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.01'
__date__ = '10 Apr 2017'


import os
import sys
import glob
import shutil
import argparse as ap
import configparser as cp
import subprocess as sb
import multiprocessing as mp
from Bio import SeqIO # Biopython (v 1.69) require NumPy (v 1.12.1)
import tempfile as tf
import bz2
# import tarfile
# import collections
# try:
#     import cPickle as pickle
# except:
#     import pickle
# import urllib2
# from contextlib import closing
# from StringIO import StringIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# import taxcuration as taxc
# import time
# from itertools import chain
# import hashlib
# try:
#     from collections import Counter # works only with python >= 2.7
#     collections_counter = True
# except:
#     collections_counter = False


config_sections_all = ['markers_db', 'dna_map', 'aa_map', 'msa', 'tree']
config_options_all = ['program_name', 'program_name_parallel', 'params', 'threads', 'input', 'database', 'output_path', 'output', 'version']
config_options_mandatory = ['program_name']
# download = ""
# ppa_fna = "data/ppa.seeds.faa"
# ppa_fna_40 = "data/ppa.seeds.40.faa"
# ppa_aln = "data/ppafull.aln.faa"
# ppa_up2prots = "data/ppafull.up2prots.txt"
# ppa_ors2prots = "data/ppafull.orgs2prots.txt"
# ppa_tax = "data/ppafull.tax.txt"
# ppa_alns = ("data/ppaalns/list.txt", "data/ppaalns/ppa.aln.tar.bz2")
# ppa_alns_fol = "data/ppaalns/"
# ppa_xml = "data/ppafull.xml"
# ppa_udb = "data/ppa.udb"
# up2prots = "up2prots.txt"
# ors2prots = "orgs2prots.txt"
# aln_tot = "aln.fna"
# aln_int_tot = "aln.int.fna"
# ups2faa_pkl = "ups2faa.pkl"
# o_tree = "tree.nwk"
# o_inttree = "tree.int.nwk"
# p2t_map = "p2t_map.txt"
# cleanfaa_fld = "clean_faa/"
# old2newprots = "old2newprots.txt"
# NOT_ENOUGH_MAPPINGS = "Not enough mappings"
# few_maps = "few_maps.txt"
# f_protein_ids_map = "protein_ids_map.pkl"


def info(s, init_new_line=False, exit=False, exit_value=0):
    inf = '\n' if init_new_line else ''
    sys.stdout.write('{}{}'.format(inf, s))
    sys.stdout.flush()
    sys.exit(exit_value) if exit else None


def error(s, init_new_line=False, exit=False, exit_value=1):
    err = '\n' if init_new_line else ''
    err += '[e] '
    sys.stderr.write('{}{}\n'.format(err, s))
    sys.stderr.flush()
    sys.exit(exit_value) if exit else None


def read_params():
    p = ap.ArgumentParser(description="")

    group = p.add_mutually_exclusive_group()
    group.add_argument('-i', '--integrate', metavar='PROJECT_NAME', type=str, default=None, help="Integrate user genomes into the PhyloPhlAn tree")
    group.add_argument('-u', '--user_tree', metavar='PROJECT_NAME', type=str, default=None, help="Build a phylogenetic tree using only genomes provided by the user")
    group.add_argument('-c', '--clean', metavar='PROJECT_NAME', type=str, default=None, help="Clean the final and partial data produced for the specified project")

    p.add_argument('-d', '--database', type=str, default=None, required=False, help="The name of the database to use")

    group = p.add_mutually_exclusive_group()
    group.add_argument('--strain', action='store_true', default=False, help="")
    group.add_argument('--tol', action='store_true', default=False, help="")
    group.add_argument('--meta', action='store_true', default=False, help="")

    group = p.add_argument_group(title="Folders location", description="Parameters for setting different folders location")
    group.add_argument('--input_folder', type=str, default='input/', help="Path to the folder containing the folder with the input data, default input/")
    group.add_argument('--configs_folder', type=str, default='configs/', help="Path to the folder containing the configuration files for the external software, default configs/")
    group.add_argument('--data_folder', type=str, default='data/', help="Path to the folder where to store the intermediate files, default data/")
    group.add_argument('--databases_folder', type=str, default='databases/', help="Path to the folder where to store the intermediate files, default databases/")
    group.add_argument('--output_folder', type=str, default='output/', help="Path to the output folder where to save the results, default output/")

    p.add_argument('--clean_all', action='store_true', default=False, help="Remove all instalation and database file that is automatically extracted and formatted at the first pipeline run")
    p.add_argument('--database_list', action='store_true', default=False, help="If specified lists the available databases that can be specified with the --database option")

    p.add_argument('--nproc', type=int, default=1, help="The number of CPUs to use, default 1")

    # # decide which database to use with blast
    # p.add_argument('--blast_full', action='store_true', default=False, help=
    #     "If specified, tells blast to use the full dataset of universal proteins. As default the small "
    #     "dataset of universal proteins is used.\n")
    # # decide if you want perform the faa cleaning of the .faa
    # p.add_argument('--faa_cleaning', action='store_true', default=False, help=
    #     "When specified perform a cleaning on the number and the length of proteins, changing also the "
    #     "proteins id such that are unique among all the genomes. As default this step will be skipped.\n")

    p.add_argument('--min_num_proteins', type=int, default=100, help="Proteomes (.faa) with less than this number of proteins will be discarded, default is 100")
    p.add_argument('--min_len_protein', type=int, default=50, help="Proteins in proteomes (.faa) shorter than this value will be discarded, default is 50\n")
    p.add_argument('--min_num_markers', type=int, default=0, help="Inputs that map less than this number of markers will be discarded, default is 0, i.e., no input will be discarded")

    group = p.add_argument_group(title="Filename extensions", description="Parameters for setting the extensions of the input files")
    group.add_argument('--genome_extension', type=str, default='.fna', help="Set the extension for the genomes in your inputs, default .fna")
    group.add_argument('--proteome_extension', type=str, default='.faa', help="Set the extension for the proteomes in your inputs, default .faa")

    p.add_argument('--verbose', action='store_true', default=False, help="Makes PhyloPhlAn2 verbose")
    p.add_argument('-v', '--version', action='store_true', default=False, help="Prints the current PhyloPhlAn2 version")

    return p.parse_args()


def read_configs(configs_folder, verbose=False):
    configs = {}

    for config_file in glob.iglob(configs_folder+"*.config"):
        info('Reading configuration file {}\n'.format(config_file)) if verbose else None
        config = cp.ConfigParser()
        config.read(config_file)

        for section in config_sections_all:
            if section in config.sections(): # "DEFAULT" section not included!
                configs[section] = {}

                for option in config_options_all:
                    if option in config[section]:
                        configs[section][option] = config[section][option]
                    else:
                        configs[section][option] = None
            else:
                configs[section] = None

    return configs


def check_args(args, verbose):
    if args.version:
        info('PhyloPhlAn2 version {} ({})\n'.format(__version__, __date__), exit=True)

    if args.database_list:
        info('Available databases:\n    {}\n'.format('\n    '.join(database_list(args.databases_folder))), exit=True)

    if args.clean_all:
        return

    args.input_folder += '/' if not args.input_folder.endswith('/') else ''
    args.configs_folder += '/' if not args.configs_folder.endswith('/') else ''
    args.data_folder += '/' if not args.data_folder.endswith('/') else ''
    args.databases_folder += '/' if not args.databases_folder.endswith('/') else ''
    args.output_folder += '/' if not args.output_folder.endswith('/') else ''

    if not os.path.isdir(args.input_folder):
        error('{} folder does not exists'.format(args.input_folder), exit=True)

    if not os.path.isdir(args.configs_folder):
        error('{} folder does not exists'.format(args.configs_folder), exit=True)

    if not os.path.isdir(args.data_folder):
        error('{} folder does not exists'.format(args.data_folder), exit=True)

    if not os.path.isdir(args.databases_folder):
        error('{} folder does not exists'.format(args.databases_folder), exit=True)

    if not os.path.isdir(args.output_folder):
        error('{} folder does not exists'.format(args.output_folder), exit=True)

    if (not args.integrate) and (not args.user_tree) and (not args.clean):
        error('Either -i (--integrate), or -u (--user_tree), or -c (--clean) must be specified', exit=True)

    project_name = args.integrate if args.integrate else None
    project_name = args.user_tree if args.user_tree else project_name
    project_name = args.clean if args.clean else project_name
    project_name += '/' if not project_name.endswith('/') else ''

    args.input_folder += project_name
    args.data_folder += project_name
    args.output_folder += project_name

    if not os.path.isdir(args.input_folder):
        error('Project folder does not exists in input directory ({})'.format(args.input_folder), exit=True)

    if not os.path.isdir(args.data_folder):
        info('Creating folder {}\n'.format(args.data_folder)) if verbose else None
        os.mkdir(args.data_folder, mode=0o775)

    if not os.path.isdir(args.output_folder):
        info('Creating folder {}\n'.format(args.output_folder)) if verbose else None
        os.mkdir(args.output_folder, mode=0o775)

    if args.database:
        if not os.path.isdir(args.databases_folder+args.database):
            error('{} folder does not exists\nAvailable databases are: {}'.format(args.databases_folder+args.database, ', '.join(database_list(args.databases_folder))), exit=True)
    else:
        error('Database not selected\nAvailable databases are: {}'.format(', '.join(database_list(args.databases_folder))), exit=True)

    args.genome_extension = '.'+args.genome_extension if not args.genome_extension.startswith('.') else args.genome_extension
    args.genome_extension = args.genome_extension[:-1] if args.genome_extension.endswith('.') else args.genome_extension

    args.proteome_extension = '.'+args.proteome_extension if not args.proteome_extension.startswith('.') else args.proteome_extension
    args.proteome_extension = args.proteome_extension[:-1] if args.proteome_extension.endswith('.') else args.proteome_extension

    if args.strain: # params for strain-level phylogenies
        print('>  ARGS.STRAIN  <')
        pass
    elif args.tol: # params for tree-of-life phylogenies
        print('>  ARGS.TOL  <')
        pass
    elif args.meta: # params for phylogenetic placement of metagenomic contigs
        print('>  ARGS.META  <')
        pass
    else: # default params
        print('>  DEFAULT  <')
        pass


def check_configs(configs, verbose=False):
    for section in config_sections_all:
        info('Checking {} section in configuration file\n'.format(section)) if verbose else None
        error('Could not find {} section in configuration file'.format(section), exit=True) if section not in configs else None
        for option in config_options_mandatory:
            error('Could not find {} mandatory option in section {} in configuration file'.format(option, section), exit=True) if option not in configs[section] else None


def check_dependencies(configs, nproc, verbose=False):
    for prog in [compose_command(configs[params], check=True, nproc=nproc) for params in configs]:
        try:
            info('Checking {}\n'.format(' '.join(prog))) if verbose else None
            # sb.run(prog, stdout=sb.DEVNULL, stderr=sb.DEVNULL, check=True)
            print('>  RICORDA DI DECOMMENTARE IL RUN  <')
        except sb.CalledProcessError as cpe:
            error('{} returned the following error: {}'.format(' '.join(prog), cpe), exit=True)
        except:
            error('{} not installed or not present in system path'.format(' '.join(prog)), exit=True)


def database_list(databases_folder):
    return os.listdir(databases_folder)


def compose_command(params, check=False, input_file=None, database=None, output_file=None, nproc=1):
    command = [params['program_name']]

    if (nproc > 1) and params['program_name_parallel']:
        command = [params['program_name_parallel']]

    if params['threads']:
        command.append(params['threads'])
        command.append(nproc)

    if check and params['version']:
        command.append(params['version'])
    elif params['params']:
        command += params['params'].split(' ')

    if input_file:
        if params['input']:
            command.append(params['input'])

        command.append(input_file)

    if database:
        if params['database']:
            command.append(params['database'])

        command.append(database)

    if output_file:
        if params['output']:
            command.append(params['output'])

        command.append(output_file)

    return [str(a) for a in command]


def init_database(database, databases_folder, command, verbose=False):
    db_fasta, db_output = None, None

    if os.path.isfile(databases_folder+database): # assumed to be a fasta file containing the markers
        db_fasta = databases_folder+database
        db_output = databases_folder+database[:database.rfind('.')]+'.udb'
    elif os.path.isdir(databases_folder+database): # assumed to be a folder with a fasta file for each marker
        db_fasta = databases_folder+database+'/'+database+'.fna'
        db_output = databases_folder+database+'/'+database+'.udb'

        if not os.path.isfile(db_fasta):
            with open(db_fasta, 'w') as f:
                for i in glob.iglob(databases_folder+database+'/*'):
                    with open(i) as g:
                        f.write(g.read())
        elif verbose:
            info('File {} already present\n'.format(db_fasta))
    else: # what's that??
        error("Custom set of markers not recognize", exit=True)

    if db_fasta and db_output:
        if not os.path.isfile(db_output):
            try:
                info("Generating usearch indexed DB {}\n".format(db_output))
                cmd = compose_command(command, input_file=db_fasta, output_file=db_output)
                # sb.run(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL, check=True))
                print('>  RICORDA DI DECOMMENTARE IL RUN  <')
            except sb.CalledProcessError as cpe:
                error('{} returned the following error: {}'.format(' '.join(cmd), cpe), exit=True)
            except:
                error('{} not installed or not present in system path'.format(' '.join(cmd)), exit=True)
        elif verbose:
            info('usearch database {} already present\n'.format(db_output))

    return db_output


def clean_all(loc_dat, verbose=False):
    print('>  DA FARE, BISOGNA CAPIRE QUALE SARA` LO STATO INIZIALE DEL PROGETTO, DATABASES, ETC...  <')
    sys.exit(0)


def clean_project(data_folder, output_folder, verbose=False):
    if os.path.exists(data_folder):
        info('Cleaning prject folder {}\n'.format(data_folder)) if verbose else None
        shutil.rmtree(data_folder)

    if os.path.exists(output_folder):
        info('Cleaning prject folder {}\n'.format(output_folder)) if verbose else None
        shutil.rmtree(output_folder)

    sys.exit(0)


def load_input_files(input_folder, extension, verbose=False):
    inputs, inputs_compressed = [], []
    files = glob.iglob(input_folder+'*'+extension+'*')

    for file in files:
        if file.endswith('.bz2'):
            tmp_faa = tf.NamedTemporaryFile(suffix=extension, prefix=file[file.rfind('/')+1:file.find('.')], dir=file[:file.rfind('/')+1])

            with open(tmp_faa.name, 'w') as f:
                with bz2.open(file, 'rt') as g:
                    records = list(SeqIO.parse(g, "fasta"))

                SeqIO.write(records, f, "fasta")

            inputs_compressed.append(tmp_faa)
            inputs.append(tmp_faa.name)
        elif file.endswith(extension):
            inputs.append(file)
        else:
            info('input file {} not recognized\n'.format(file))

    return (inputs, inputs_compressed)


def check_input_proteomes(inputs, min_num_proteins, min_len_protein, verbose=False):
    good_inputs = []

    for inp in inputs:
        info('Checking input file {}\n'.format(inp)) if verbose else None
        good_proteins = 0

        for seq_record in SeqIO.parse(inp, "fasta"):
            if len(seq_record) >= min_len_protein:
                good_proteins += 1

        if good_proteins >= min_num_proteins:
            good_inputs.append(inp)
        elif verbose:
            info('Input {} dscarded, not enough good proteins\n'.format(inp))

    return good_inputs


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def gene_markers_identification(configs, key, inputs, output_folder, database_name, database, nproc=1, verbose=False):
    commands = []

    if not os.path.isdir(output_folder):
        info('Creating folder {}\n'.format(output_folder)) if verbose else None
        os.mkdir(output_folder)
    elif verbose:
        info('Folder {} already exists\n'.format(output_folder))

    for inp in inputs:
        out = inp[inp.rfind('/')+1:]
        out = output_folder+out[:out.rfind('.')]+'.b6o'

        if not os.path.isfile(out):
            commands.append((configs[key], inp, database, out))

    if commands:
        info('Mapping {} on {} inputs (key: {})\n'.format(database_name, len(commands), key))
        pool_error = False
        terminating = mp.Event()
        pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

        try:
            pool.imap_unordered(gene_markers_identification_rec, commands, chunksize=nproc) # if there are a lot of cpus it is reasonable assume that there will be a lot of inputs too!
            pool.close()
        except:
            pool.terminate()
            pool_error = True

        pool.join()
        error('gene_markers_identification crashes', exit=True) if pool_error else None
    else:
        info('{} markers already mapped (key: {})\n'.format(database_name, key))


def gene_markers_identification_rec(x):
    if not terminating.is_set():
        try:
            params, inp, db, out = x
            cmd = compose_command(params, input_file=inp, database=db, output_file=out)
            info('Starting {}\n'.format(inp))
            # sb.run(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL, check=True))
            print('>  RICORDA DI DECOMMENTARE IL RUN  <')
            # shutil.copy2(out, out+'.bkp') # save the original results
            print('>  RICORDA DI DECOMMENTARE LA COPIA  <')

            print('>  SCREENING OF THE RESULTS OF USEARCH  <')
        except sb.CalledProcessError as cpe:
            error('{} returned the following error: {}'.format(' '.join(cmd), cpe))
            terminating.set()
            return
        except:
            error('Cannot execute command {}'.format(' '.join(cmd)))
            terminating.set()
            return

        info('{} generated\n'.format(out))
    else:
        terminating.set()


def gene_markers_extraction(input_folder, output_folder, nproc=1, verbose=False):
    if not os.path.isdir(output_folder):
        info('Creating folder {}\n'.format(output_folder)) if verbose else None
        os.mkdir(output_folder)
    elif verbose:
        info('Folder {} already exists\n'.format(output_folder))

    info('Extracting markers from inputs\n')
    pool_error = False
    terminating = mp.Event()
    pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

    try:
        pool.imap_unordered(gene_markers_extraction_rec, glob.iglob(output_folder+'*.b6o'), chunksize=nproc) # if there are a lot of cpus it is reasonable assume that there will be a lot of inputs too!
        pool.close()
    except:
        pool.terminate()
        pool_error = True

    pool.join()
    error('gene_markers_extraction crashes', exit=True) if pool_error else None


def gene_markers_extraction_rec(x):
    if not terminating.is_set():
        try:
            info('Extracting {}\n'.format(x))
        except:
            error('Cannot execute command {}'.format(' '.join(cmd)))
            terminating.set()
            return

        info('{} generated\n'.format(out))
    else:
        terminating.set()


# def get_inputs(proj, params):
#     if not os.path.isdir(inp_fol):
#         error("Folder not found: '"+inp_fol+"'", init_new_line=True, exit=True)

#     faa_in = []
#     fna_in = []
#     txt_in = []
#     tax_in = []
#     mdt_in = []

#     for f in iglob(inp_fol+'*'):
#         if 'faa' in f.split('.'):
#             cut = -2 if f.endswith('.bz2') else -1
#             faa_in.append('.'.join(f[f.rfind('/')+1:].split('.')[:cut]))

#         if 'fna' in f.split('.'):
#             cut = -2 if f.endswith('.bz2') else -1
#             fna_in.append('.'.join(f[f.rfind('/')+1:].split('.')[:cut]))

#         if os.path.splitext(f)[1] == ".txt":
#             txt_in.append(f[f.rfind('/')+1:])

#         if os.path.splitext(f)[1] == ".tax":
#             tax_in.append(f[f.rfind('/')+1:])

#         if os.path.splitext(f)[1] == ".metadata":
#             mdt_in.append(f[f.rfind('/')+1:])

#     if not (len(faa_in) + len(fna_in)):
#         error("No '.faa' or '.fna' input files found in '"+inp_fol+"'", init_new_line=True, exit=True)

#     if len(txt_in) > 1:
#         error("More than one '.txt' input files found in input/\n"
#               "[No more than one txt file (the taxonomy) allowed", init_new_line=True, exit=True)

#     if len(tax_in) > 1:
#         error("More than one '.tax' input files found in input/\n"
#               "[No more than one tax file (the taxonomy for taxonomic analysis) allowed", init_new_line=True, exit=True)

#     if len(tax_in) > 1:
#         error("More than one '.metadata' input files found in input/\n"
#               "[No more than one metadata file allowed", init_new_line=True, exit=True)

#     if not os.path.isdir(dat_fol):
#         os.mkdir(dat_fol)

#     # clean the faa files by filter on the number of proteins and their minimum length
#     if faa_in and params['faa_cleaning']:
#         t0 = time.time()
#         faa_in = clean_faa(proj, faa_in, params['min_num_proteins'], params['min_len_protein']) # , verbose=True)
#         info("\nCleaning faa inputs took "+str(int(time.time()-t0))+" s\n")

#     fna_in = uniq_filenames(faa_in, fna_in) # make uniq filename for fna inputs
#     if not os.path.isfile(dat_fol+ors2prots):
#         with open(dat_fol+ors2prots, 'w') as outf:
#             for f in faa_in:
#                 fld = inp_fol if not params['faa_cleaning'] else dat_cln_fol

#                 if os.path.isfile(fld+f+'.faa'):
#                     prots = sorted([l.split()[0][1:] for l in open(fld+f+'.faa') if l.startswith('>')])
#                 else:
#                     with BZ2File(fld+f+'.faa.bz2', 'rU') as inp:
#                         prots = sorted([l.split()[0][1:] for l in inp if l.startswith('>')])

#                 outf.write('\t'.join([f] + prots) + '\n')

#             for f, fo in fna_in:
#                 if os.path.isfile(inp_fol+f+'.fna'):
#                     prots = sorted([l.split()[0][1:] for l in open(inp_fol+f+'.fna') if l.startswith('>')])
#                 else:
#                     with BZ2File(inp_fol+f+'.fna.bz2', 'rU') as inp:
#                         prots = sorted([l.split()[0][1:] for l in inp if l.startswith('>')])

#                 outf.write('\t'.join([fo] + prots) + '\n')

#     return faa_in, fna_in, txt_in[0] if txt_in else None, tax_in[0] if tax_in else None, mdt_in[0] if mdt_in else None


# def check_inp(inps):
#     init_dig = [l for l in inps if l[0].isdigit()]

#     if init_dig:
#         error("The following genome IDs start with a digit:\n    "+"\n    ".join(init_dig) +
#               "\nPlease rename accordingly the corresponding input file names", init_new_line=True, exit=True)

#     ppa_ids = [l[1:] for l in open(ppa_aln) if l.startswith('>')]
#     init_dup = [l for l in inps if l in ppa_ids]

#     if init_dup:
#         error("The following genome IDs are already in PhyloPhlAn:\n    "+"\n    ".join(init_dig),
#               init_new_line=True, exit=True)


# def clean_faa(proj, faa_in, min_num_proteins, min_len_protein, verbose=False):
#     if not os.path.isdir(dat_cln_fol): # create the directory if it does not exists
#         os.mkdir(dat_cln_fol)
#     elif os.path.isfile(dat_fol+f_protein_ids_map): # load the protein_ids_map mapping and return it
#         with open(dat_fol+f_protein_ids_map) as f:
#             return pickle.load(f).keys()

#     protein_ids_map = dict() # will contains the old protein id and the new protein id assigned

#     for f in faa_in: # parallelize?
#         records = []
#         out = []
#         seqid = 0
#         msg = ""
#         old_id = dict()

#         if os.path.isfile(inp_fol+f+'.faa'):
#             with open(inp_fol+f+'.faa', 'rU') as h_in:
#                 records = list(SeqIO.parse(h_in, "fasta"))
#         elif os.path.isfile(inp_fol+f+'.faa.bz2'):
#             with BZ2File(inp_fol+f+'.faa.bz2', 'rU') as h_in:
#                 records = list(SeqIO.parse(h_in, "fasta"))
#         else:
#             info(' '.join(["[clean_faa]", "File:", f, "not found (.faa or .faa.bz2)"]) + "\n")
#             continue

#         if len(records) >= min_num_proteins: # check the number of proteins of the genome
#             for protein in records:
#                 if len(protein.seq) >= min_len_protein: # check the length of the protein
#                     pid = f+'_'+str(seqid)
#                     out.append(SeqRecord(protein.seq, id=pid, description=''))
#                     old_id[pid] = protein.id
#                     seqid += 1
#         else:
#             msg = "not enough proteins!"

#         # if some proteins has been removed, checked that the good ones kept are at least min_num_proteins
#         if len(out) >= min_num_proteins:
#             if len(records) - len(out):
#                 msg = str(len(records)-len(out)) + " too short proteins removed!"

#             # save the old and new mapping
#             if f in protein_ids_map:
#                 app = protein_ids_map[f].append()
#                 for t in ((old_id[s.id], s.id) for s in out):
#                     app(t)
#             else:
#                 protein_ids_map[f] = [(old_id[s.id], s.id) for s in out]

#             # write the new .faa
#             with open(dat_cln_fol+f+'.faa', 'w') as h_out:
#                 SeqIO.write(out, h_out, "fasta")
#         else:
#             msg = "not enough good proteins!"

#         if msg and verbose:
#             info(' '.join(["[clean_faa]", f, msg]) + "\n")

#     # write the mapping of the old and new proteins id of all genomes
#     with open(dat_fol+old2newprots, 'w') as f:
#         for k, v in protein_ids_map.iteritems():
#             f.write('\t'.join([k] + ['('+a+','+b+')' for a,b in v]) + '\n')

#     # serialize the protein_ids_map mapping
#     with open(dat_fol+f_protein_ids_map, 'w') as f:
#         pickle.dump(protein_ids_map, f, pickle.HIGHEST_PROTOCOL)

#     return protein_ids_map.keys()


# def initu(terminating_, loc_inp_, dat_fol_, fin_dict_, fake_proteomes_):
#     global terminating
#     global loc_inp
#     global dat_fol
#     global fin_dict
#     global fake_proteomes
#     terminating = terminating_
#     loc_inp = loc_inp_
#     dat_fol = dat_fol_
#     fin_dict = fin_dict_
#     fake_proteomes = fake_proteomes_


# def exe_usearch(x):

#     def screen_usearch_wdb(inpf):
#         tab = ((l[0].split(' ')[0], l[1].split('_')[1], float(l[-1])) for l in (ll.strip().split('\t') for ll in open(inpf)))
#         f2p, f2b, p2b = {}, {}, {}

#         for f,p,b in tab:
#             if f in f2p and b < f2b[f]:
#                 continue
#             if p in p2b and b < p2b[p]:
#                 continue
#             if p in p2b:
#                 for ff in f2p.keys():
#                     if p in f2p[ff]:
#                         del f2p[ff]
#                         del f2b[ff]

#             f2p[f], f2b[f], p2b[p] = p, b, b

#         with open(inpf, 'w') as outf:
#             # there should be at least min_uprots of universal proteins mapped
#             if len(f2p) >= min_uprots:
#                 for k, v in f2p.iteritems():
#                     outf.write("\t".join([str(k), str(v)]) + "\n")
#             else:
#                 outf.write("\t".join([inpf[inpf.rfind('/')+1:], NOT_ENOUGH_MAPPINGS]) + "\n")

#     #
#     if not terminating.is_set():
#         try:
#             # info("Starting "+x[7][x[7].rfind('/')+1:]+"\n")
#             info("Starting "+x[7]+"\n")
#             tmp_faa = None
#             t = None

#             if not os.path.isfile(x[3]):
#                 tmp_faa = NamedTemporaryFile(suffix='.faa', prefix=x[3][x[3].rfind('/')+1:x[3].find('.')], dir=x[3][:x[3].rfind('/')+1])

#                 with open(tmp_faa.name, 'w') as h_out:
#                     with BZ2File(x[3]+'.bz2', 'rU') as h_in:
#                         records = list(SeqIO.parse(h_in, "fasta"))

#                     SeqIO.write(records, h_out, "fasta")

#                 cmd = x[:3] + [tmp_faa.name] + x[4:]
#             else:
#                 cmd = x

#             t = sb.call(cmd)

#             if t:
#                 t.wait()

#             shutil.copy2(x[7], x[7]+'.bkp') # copy the results of usearch
#             screen_usearch_wdb(x[7])

#             if tmp_faa:
#                 tmp_faa.close()

#             # info(x[7][x[7].rfind('/')+1:]+" generated!\n")
#             info(x[7]+" generated!\n")
#         except:
#             if t:
#                 t.kill()

#             if tmp_faa:
#                 tmp_faa.close()

#             terminating.set()
#             error(' '.join(x))
#             return
#     else:
#         terminating.set()

#     return


# def faa2ppafaa(inps, nproc, proj, faa_cleaning):
#     if not os.path.isdir(dat_fol): # create the data folder if does not exists
#         os.mkdir(dat_fol)

#     loc_inp = dat_cln_fol if faa_cleaning else inp_fol
#     mmap = ((loc_inp+i+'.faa', dat_fol+i+'.b6o') for i in inps if not os.path.isfile(dat_fol+i+'.b6o'))
#     us_cmd = [["usearch", "-quiet", "-ublast", i, "-db", cf_udb if cf_udb else ppa_udb, "-blast6out", o, "-threads", "1",
#             "-evalue", "1e-10", "-maxaccepts", "8", "-maxrejects", "32"] for i, o in mmap] # usearch8.0.1517

#     if not us_cmd:
#         info("All usearch runs already performed!\n")
#     else:
#         sw = 'your custom set of' if cf_up else 'PhyloPhlAn'
#         info('Looking for ' + sw + ' proteins in input .faa files\n')

#         terminating = mp.Event()
#         pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

#         try:
#             pool.map(exe_usearch, us_cmd)
#             pool.close()
#         except:
#             pool.terminate()
#             pool.join()
#             error("Quitting faa2ppafaa()")
#             return

#         pool.join()
#         info("All usearch runs performed!\n")

#     if not os.path.isfile(dat_fol+up2prots):
#         up2p = collections.defaultdict(list)
#         too_few_maps = set()

#         for i in inps:
#             for p, up in (l.strip().split('\t') for l in open(dat_fol+i+'.b6o').readlines()):
#                 if not up.startswith(NOT_ENOUGH_MAPPINGS):
#                     up2p[up].append(p)
#                 else:
#                     too_few_maps |= set([i])

#         # write the skipped genomes
#         if too_few_maps:
#             with open(dat_fol+few_maps, 'w') as fout:
#                 fout.write('\n'.join(too_few_maps) + '\n')

#         # write the mapping between universal proteins and proteins mapped
#         if up2p:
#             with open(dat_fol+up2prots, 'w') as outf:
#                 for k, v in up2p.iteritems():
#                     outf.write("\t".join([k]+v) + "\n")


# def blastx_exe(x):
#     if not terminating.is_set():
#         try:
#             # info("Starting "+x[6][x[6].rfind('/')+1:]+"\n")
#             info("Starting "+x[6]+"\n")
#             tmp_fna = None
#             t = None

#             if not os.path.isfile(x[2]):
#                 tmp_fna = NamedTemporaryFile(suffix='.fna', prefix=x[2][x[2].rfind('/')+1:x[2].find('.')], dir=x[2][:x[2].rfind('/')+1])

#                 with open(tmp_fna.name, 'w') as h_out:
#                     with BZ2File(x[2]+'.bz2', 'rU') as h_in:
#                         records = list(SeqIO.parse(h_in, "fasta"))

#                     SeqIO.write(records, h_out, "fasta")

#                 cmd = x[:2] + [tmp_fna.name] + x[3:]
#             else:
#                 cmd = x

#             with open(os.devnull, 'w') as devnull:
#                 t = sb.call(cmd, stderr=devnull) # tblastn quiet homemade!

#             if t:
#                 t.wait()

#             if tmp_fna:
#                 tmp_fna.close()

#             # info(x[6][x[6].rfind('/')+1:]+" generated!\n")
#             info(x[6]+" generated!\n")
#         except:
#             if t:
#                 t.kill()

#             if tmp_fna:
#                 tmp_fna.close()

#             terminating.set()
#             error(' '.join(x))
#             return
#     else:
#         terminating.set()

#     return


# def blast(inps, nproc, proj, blast_full=False):
#     loc_dat = dat_fol+"tblastn/"
#     terminating = mp.Event()
#     pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)
#     mmap = [(inp_fol+i+'.fna', loc_dat+io+'.b6o') for i, io in inps.iteritems() if (not os.path.isfile(loc_dat+io+'.b6o'))]

#     if not os.path.isdir(loc_dat): # create the tmp directory if does not exists
#         os.mkdir(loc_dat)

#     if not mmap:
#         info('All tblastn runs already performed!\n')
#     else:
#         sw = 'your custom set of' if cf_up else 'PhyloPhlAn'
#         info('Looking for ' + sw + ' proteins in input .fna files ('+str(len(mmap))+')\n')

#         # small dataset
#         dataset = ppa_fna_40
#         evalue = '1e-40'

#         if cf_up: # custom set of universal proteins
#             dataset = cf_fna
#         elif blast_full: # full dataset
#             dataset = ppa_fna
#             evalue = '1e-60'

#         us_cmd = [['tblastn', '-subject', i, '-query', dataset, '-out', o, '-outfmt', '6', '-evalue', evalue]
#                   for i, o in mmap] # tblastn2.2.28+

#         try:
#             pool.map(blastx_exe, us_cmd)
#             pool.close()
#         except:
#             pool.terminate()
#             pool.join()
#             error("Quitting blast()")
#             return

#         pool.join()
#         info('All tblastn runs performed!\n')


# def gens2prots(inps, proj, faa_cleaning):
#     loc_inp = dat_cln_fol if faa_cleaning else inp_fol

#     if not os.path.isdir(dat_fol): # create the tmp directory if does not exists
#         os.mkdir(dat_fol)

#     if os.path.isfile(dat_fol+ups2faa_pkl):
#         return

#     ups2prots = dict([(ll[0], set(ll[1:])) for ll in (l.strip().split('\t') for l in open(dat_fol+up2prots))])
#     prots2ups = {}

#     for k, v in ups2prots.iteritems():
#         for vv in v:
#             if vv in prots2ups:
#                 error(str(vv) + " already in dict!")
#             prots2ups[vv] = k

#     ups2faa = collections.defaultdict(set)
#     e = None

#     for i in inps:
#         inp = None

#         if os.path.isfile(loc_inp+i+'.faa'):
#             inp = open(loc_inp+i+'.faa', 'rU')
#         else:
#             inp = BZ2File(loc_inp+i+'.faa.bz2', 'rU')

#         for l in inp:
#             if l.startswith(">"):
#                 if e in prots2ups:
#                     ups2faa[prots2ups[e]].add(s+"\n")

#                 e, s = l.strip().split()[0][1:], ""

#             s += l

#         if e in prots2ups:
#             ups2faa[prots2ups[e]].add(s+"\n")

#         if inp:
#             inp.close()

#     not_enough_seqs = []

#     for k, v in ups2faa.iteritems():
#         if len(v) > 1: # write only ups with at least 2 sequences!
#             if not os.path.isfile(dat_fol+k+'.faa'):
#                 with open(dat_fol+k+'.faa', 'w') as outf:
#                     for vv in v:
#                         outf.write("".join(vv))
#         else:
#             not_enough_seqs.append(k)

#     if not_enough_seqs:
#         for k in not_enough_seqs:
#             ups2faa.pop(k, None)

#     with open(dat_fol+ups2faa_pkl, 'w') as outf:
#         pickle.dump(ups2faa, outf, pickle.HIGHEST_PROTOCOL)


# def screen(stat, cols, sf=None, unknown_fraction=0.5, n=10):
#     lena, nsc = len(cols[0]), []

#     if sf and os.path.isfile(sf):
#         with open(sf) as sfinp:
#             scores = [(ll[:12], ll[12:]) for ll in [l for l in sfinp.readlines()]]
#     else:
#         scores = [("1.0", c) for c in cols]

#     for sc, st in zip(scores, stat):
#         ssc, lsc = sc

#         try:
#             float(ssc)
#         except:
#             continue

#         if (len(st) == 1) or \
#            (len(st) == 2 and ("-" in st or "X" in st)) or \
#            (len(st) == 3 and "X" in st and "-" in st):
#             continue

#         nn = 0
#         for k in ["-","X"]:
#             if k in st:
#                 nn += st[k]

#         if nn > lena*unknown_fraction:
#             continue

#         nsc.append((float(ssc), lsc.strip()))

#     nsc = sorted(nsc, key=lambda x: x[0], reverse=True)
#     ret = [v for _, v in nsc][-n:]

#     return ret if ret else [list(['-']*lena)]


# def aln_subsample(inp_f, out_f, scores, unknown_fraction, namn):
#     fnas = list(SeqIO.parse(inp_f, "fasta"))
#     ids = [f.id for f in fnas]
#     cols = zip(*[f.seq for f in fnas])
#     col_alf = [set(x) for x in cols]
#     col_sta = [dict([(a, seq.count(a)) for a in s]) for s, seq in zip(col_alf, cols)]
#     col_sta_ok = screen(col_sta, cols, sf=scores, unknown_fraction=unknown_fraction, n=namn)
#     raws = zip(*col_sta_ok)
#     recs = [SeqRecord(seq=Seq("".join(r)), id=nid, description=nid) for r, nid in zip(raws, ids)]
#     SeqIO.write(recs, out_f, "fasta")


# def exe_muscle(x):
#     if not terminating.is_set():
#         i1, i2, i3, i4, t = None, None, None, None, None

#         try:
#             if x[0] == 'muscle':
#                 i4 = x[3]
#                 info("Running muscle on "+i4+"\n")

#                 with open(os.devnull, 'w') as devnull:
#                     t = sb.call(x[:-3], stderr=devnull) # quiet mode

#                 if t:
#                     t.wait()

#                 i1, i3 = x[5], x[-6]
#             elif x[0] == 'mafft':
#                 i4 = x[-6]
#                 o4 = x[-5]
#                 info("Running mafft on "+i4+"\n")


#                 with open(os.devnull, 'w') as devnull:
#                     with open(o4+'.refine', 'w') as f:
#                         t = sb.call(x[:-5], stdout=f, stderr=devnull) # quiet mode

#                 if t:
#                     t.wait()
#                     wait_for(o4+'.refine')
#                 else:
#                     error(' '.join(x[:-5])+', t=None')

#                 t = None

#                 # compute the score file with muscle (in any case, for the moment)
#                 with open(os.devnull, 'w') as devnull:
#                     t = sb.call(['muscle', '-in', o4+'.refine', '-out', o4, '-refine', '-scorefile', x[-4]], stderr=devnull) # quiet mode

#                 if t:
#                     t.wait()
#                     wait_for(x[-4])
#                 else:
#                     error(' '.join(['muscle', '-in', o4+'.refine', '-out', o4, '-refine', '-scorefile', x[-4]])+', t=None')

#                 i1, i3 = o4, x[-4]
#             else:
#                 error('Invalid program: "'+x[0]+'"')

#             if x[-1]:
#                 pn = 80
#             else:
#                 pn = max(int(max(int((400.0-int(x[-2][1:]))*30.0/400.0), 1)**2 /30.0), 3)

#             aln_subsample(i1, x[-3], i3, 0.1, pn)
#             info(x[-3] + " generated (from "+i4+")\n")
#         except:
#             if t:
#                 t.kill()

#             terminating.set()
#             error(' '.join(x))
#             return
#     else:
#         terminating.set()

#     return


# def faa2aln(nproc, proj, integrate, mafft):
#     ppa_fol = cf_up if cf_up else ppa_alns_fol

#     with open(dat_fol+ups2faa_pkl) as inp:
#         prots = pickle.load(inp)

#     if cf_up:
#         up_list = (i[i.rfind('/')+1:i.rfind('.')] for i in iglob(cf_up+'*')) # assumes the folder contains all and only the needed files
#     else:
#         up_list = ('p{num:04d}'.format(num=aa) for aa in range(400))

#     mmap = [(dat_fol+i+".faa", dat_fol+i+".aln", dat_fol+i+".sc",
#              ppa_fol+i+".aln", ppa_fol+i+".aln.score",
#              dat_fol+i+".int.aln", dat_fol+i+".int.sc",
#              dat_fol+i+".sub.aln", dat_fol+i+".int.sub.aln",
#              cf_up, i in prots, i) for i in up_list]

#     if mafft:
#         us_cmd = [["mafft", "--anysymbol", "--quiet", i, o, s, so, pn, up] for i,o,s,_,_,_,_,so,_,up,pres,pn in mmap if not os.path.isfile(o) and pres]
#     else:
#         us_cmd = [["muscle", "-quiet",
#                    "-in", i,
#                    "-out", o,
#                    "-scorefile", s,
#                    "-maxiters", "2",
#                    so, pn, up] for i,o,s,_,_,_,_,so,_,up,pres,pn in mmap if not os.path.isfile(o) and pres]

#     if us_cmd:
#         info("There are "+str(len(us_cmd))+" proteins to align\n")

#         if integrate:
#             for pn in [pn for _, _, _, _, _, _, _, _, _, _, pres, pn in mmap if pres]:
#                 with open(dat_fol+pn+".faa", 'a') as f:
#                     with open(ppa_fol+pn+".faa", 'r') as g:
#                         f.write(g.read())

#         terminating = mp.Event()
#         pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

#         try:
#             pool.map(exe_muscle, us_cmd)
#             pool.close()
#         except:
#             pool.terminate()
#             pool.join()
#             error("Quitting faa2aln()")
#             return

#         pool.join()
#         info("All alignments performed!\n")
#     else:
#         info("All alignments already performed!\n")

#     # if os.path.exists(dat_fol+up2prots):
#     #     info("All alignments already computed!\n")
#     #     return

#     # up2p = collections.defaultdict( list )
#     # for i in inps:
#     #     for p,up in (l.strip().split('\t')
#     #             for l in open(dat_fol+i+".b6o").readlines()):
#     #         up2p[up].append( p )
#     # if integrate:
#     #     with open(ppa_up2prots) as inpf:
#     #         for l in (ll.strip().split('\t') for ll in inpf):
#     #             up2p[l[0]] += l[1:]

#     # with open(dat_fol+up2prots,"w") as outf:
#     #     for k,v in up2p.items():
#     #         outf.write( "\t".join([k]+v) + "\n" )


# def aln_merge(proj, integrate):
#     loc_dat = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot

#     if os.path.isfile(loc_dat):
#         info("All alignments already merged!\n")
#         return

#     up2p = dict([(l[0], set(l[1:])) for l in (ll.strip().split('\t') for ll in open(dat_fol+up2prots))])

#     if integrate:
#         for l in (ll.strip().split('\t') for ll in open(ppa_up2prots)):
#             if l[0] in up2p:
#                 up2p[l[0]] |= set(l[1:])
#             else:
#                 up2p[l[0]] = set(l[1:])

#     # check if there are proteins id duplicated and if yes warn the user and exit
#     if collections_counter:
#         print "TEST", [k for k, v in Counter(sum([list(x) for x in up2p.values()], [])).iteritems() if v > 1]

#     all_prots = set.union(*up2p.values()) # all proteins id
#     genomes_to_skip = [r.strip() for r in open(dat_fol+few_maps, 'r')] if os.path.isfile(dat_fol+few_maps) else []
#     t2p = dict([(l[0], set(l[1:])) for l in (ll.strip().split('\t') for ll in open(dat_fol+ors2prots)) if l[0] not in genomes_to_skip])

#     if integrate:
#         for l in (ll.strip().split('\t') for ll in open(ppa_ors2prots)):
#             t2p[l[0]] = set(l[1:])

#     # map proteins id to genomes id
#     p2t = dict()
#     if os.path.isfile(dat_fol+p2t_map):
#         p2t = dict([(l[0], l[1]) for l in (row.strip().split('\t') for row in open(dat_fol+p2t_map, 'r'))])

#     for t, p in t2p.iteritems():
#         for pp in p:
#             if pp in all_prots: # if the protein is present in the current alignments
#                 p2t[pp] = t  # map proteins to their genome

#     all_g = set(t2p.keys()) # all genomes id
#     aln = dict([(t, "") for t in t2p.keys()]) # dictionary that will contains all alignments
#     up = dict([(t, 0) for t in t2p.keys()]) # dictionary that will contains the number of universal proteins mapped

#     for f in sorted((dat_fol+ff for ff in os.listdir(dat_fol) if ff.endswith(".sub.aln"))):
#         g2aln = dict()
#         lenal = 0 # keep the lenght of the current alignments

#         for seq in SeqIO.parse(f, "fasta"):
#             if not lenal:
#                 lenal = len(seq.seq)

#             g2aln[p2t[seq.id]] = seq

#         for g in all_g: # for all genomes
#             if g in g2aln: # if there is alignment
#                 up[g] += 1
#                 aln[g] += g2aln[g]
#             else: # otherwise add gaps
#                 aln[g] += SeqRecord(Seq('-'*lenal), id=str(g))

#     out_faas = []
#     for k, v in aln.iteritems():
#         v.id = k
#         v.description = ""
#         out_faas.append(v)

#     SeqIO.write(out_faas, loc_dat, "fasta")
#     info("All alignments merged into "+loc_dat+"!\n")

#     if not os.path.isdir(out_fol):
#         os.mkdir(out_fol) # create output directory if does not exists

#     # write statistics file
#     with open(out_fol+'aln_stats.tsv', 'w') as f:
#         f.write('\t'.join(['id', 'tot_length', 'aln_length', 'tot_gaps', 'up_mapped'])+'\n')

#         for k, v in aln.iteritems():
#             tot_len = len(v)
#             tot_gaps = str(v.seq).count('-')
#             f.write('\t'.join([str(k), str(tot_len), str(tot_len-tot_gaps), str(tot_gaps), str(up[k])])+'\n')

#     info(out_fol+'aln_stats.tsv generated!\n')


# def build_phylo_tree(proj, integrate, nproc, raxml):
#     loc_out = proj+"."+ (o_inttree if integrate else o_tree)

#     if os.path.isfile(loc_out):
#         info("Final tree already built ("+loc_out+")!\n")
#         return

#     if not os.path.isdir(out_fol):
#         os.mkdir(out_fol)

#     aln_in = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot
#     info("Building the phylogenetic tree with ")

#     if raxml:
#         if nproc > 1:
#             info("raxmlHPC-PTHREADS-SSE3\n")
#             cmd = ["raxmlHPC-PTHREADS-SSE3", "-T", str(nproc)]
#         else:
#             info("raxmlHPC\n")
#             cmd = ["raxmlHPC"]

#         cmd += ["-m", "PROTCATWAG", "-n", loc_out, "-s", aln_in, "-w", os.getcwd()+'/'+out_fol, "-p", "1989"]
#     else:
#         info("FastTree\n")
#         cmd = ["FastTree", "-quiet", "-fastest", "-mlnni", "4", "-spr", "4", "-mlacc", "2", "-slownni", "-no2nd", "-out", out_fol+loc_out, aln_in]

#     t = None

#     with open(os.devnull, 'w') as devnull:
#         t = sb.call(cmd, stdout=devnull, stderr=devnull) # homemade quiet mode

#     if t:
#         t.wait()

#     info("Tree built! The output newick file is in "+out_fol+"\n")


# def circlader(proj, integrate, tax=None):
#     # inp_fol = cf_input if cf_input else "input/"
#     # inp_fol += proj+"/"
#     # dat_fol = cf_data if cf_data else "data/"
#     # dat_fol += proj+"/"
#     # out_fol = cf_output if cf_output else "output/"
#     # out_fol += proj+"/"
#     # out_img = "output/"+proj+"/imgs/"

#     tree = out_fol+proj+"."+ (o_inttree if integrate else o_tree)
#     xtree_rr = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"xml"
#     xtree_an = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
#     pngtree = out_fol+proj+".tree"+ (".int" if integrate else "")+".png"
#     annotation_f = dat_fol+proj+(".int" if integrate else "")+".annot.txt"
#     archaea_f = dat_fol+proj+(".int" if integrate else "")+".archaea.txt"
#     tax_d = [l.strip().split('\t') for l in open(inp_fol+tax)] if tax else []
#     # o2t = dict([(o,t.split('.')) for t,o in tax_d])
#     t2o = dict(tax_d)

#     if os.path.isfile(pngtree):
#         info("Output tree images already generated!\n")
#         return

#     # info("Generating tree output images ...")
#     c_circ = "circlader/circlader.py"
#     c_circ_an = "circlader/circlader_annotate.py"
#     c_reroot = "./taxcuration/tree_reroot.py"
#     names = [l.strip().split('\t')[0] for l in open(dat_fol+ors2prots)]

#     if integrate:
#         return

#         int_names = [l.strip().split('\t')[0] for l in open(ppa_ors2prots)]
#         mt4tax_d = [l.strip().split('\t')[::-1] for l in open(ppa_tax)]
#         # mt4o2t = dict([(o,t.split('.')) for t,o in mt4tax_d])
#         mt4t2o = dict(mt4tax_d)
#         mt4archaea = [o for t,o in mt4t2o.iteritems() if t.split(".")[0] == "d__Archaea"]

#         with open(archaea_f, 'w') as out:
#             out.write("\n".join(mt4archaea) +"\n")

#         sb.call([c_reroot, "-s", "lca", "-f", archaea_f, tree, xtree_rr])
#         annotation_f = dat_fol+proj+".annot.txt"

#         with open(annotation_f, 'w') as outf:
#             for v1,v2 in [('clade_fill_color','r')]:
#                 outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n")

#             for v1,v2 in [('clade_size','7.0'), ]:
#                 outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names+int_names]) + "\n")

#             for v1,v2 in [('clade_size','2.0')]:
#                 outf.write("\t".join(["*",v1,v2]) + "\n")

#             for v1,v2 in [('sub_branches_angle_reduction','0.0'), ('clade_edge_line_width','0.2'), ('branch_line_width','0.45'), ]:
#                 outf.write("\t".join([v1,v2]) + "\n")
#     else:
#         n = 1 if len(names) < 4 else int(round(min(len(names)*0.05,100)))
#         archaea = [o for t,o in t2o.iteritems() if t.split(".")[0] == "d__Archaea"]

#         if t2o and archaea:
#             with open(archaea_f, 'w') as out:
#                 out.write("\n".join(archaea) +"\n")

#             sb.call([c_reroot, "-s", "lca", "-f", archaea_f, tree, xtree_rr])
#         elif n == 1 and len(names) < 8:
#             sb.call([c_reroot, "-s","longest_edge", tree, xtree_rr])
#         elif n == 1:
#             sb.call([c_reroot, "-s", "longest_internal_edge", tree, xtree_rr])
#         else:
#             sb.call([c_reroot, "-s", "longest_internal_edge_n", "-n", str(n), tree, xtree_rr])

#         return

#         with open(annotation_f, 'w') as outf:
#             for v1,v2 in [('clade_fill_color','r'), ('clade_size','7.0'), ]:
#                 outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n")

#             for v1,v2 in [('clade_size','2.0')]:
#                 outf.write("\t".join(["*",v1,v2]) + "\n")

#             for v1,v2 in [('sub_branches_angle_reduction','0.0'), ('clade_edge_line_width','0.2'), ('branch_line_width','0.45'), ]:
#                 outf.write("\t".join([v1,v2]) + "\n")

#     sb.call([c_circ_an, "--annot", annotation_f, xtree_rr, xtree_an])
#     sb.call([c_circ, "--dpi", "300", xtree_an, pngtree])
#     info("Done!\n")


# def tax_curation(proj, tax, integrate, mtdt=None, inps=None):
#     # inp_fol = cf_input if cf_input else "input/"
#     # inp_fol += proj+"/"
#     # out_fol = cf_output if cf_output else "output/"
#     # out_fol += proj+"/"
#     taxin = ppa_tax if integrate else inp_fol+tax
#     taxout = None if integrate else out_fol+"new."+tax
#     taxdiffs = None if integrate else out_fol+"diff."+tax
#     taxboth = None if integrate else out_fol+"comp."+tax
#     intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
#     operation = "imputation" if integrate else "curation"

#     if not os.path.isfile(intree):
#         intree = out_fol+proj+".tree."+ ("int." if integrate else "")+"nwk"

#     if tax:
#         info("Starting taxonomic curation: detecting suspicious taxonomic labels... ")
#         taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, to_integrate=inp_fol+tax)
#         taxCleaner.remove_tlabels()
#         info("Done!\n")
#         info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy... ")
#     else:
#         taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, inps=inps)
#         info("Trying to impute taxonomic labels for taxa newly integrated into the tree... ")

#     taxCleaner.infer_tlabels()
#     info("Done!\n")

#     info("Writing taxonomic "+operation+" outputs... ")
#     tid2img = taxCleaner.write_output(proj, images=False)
#     info("Done!\n")
#     info("Writing final pdf reports... ")
#     if integrate:
#         pass
#     elif mtdt:
#         taxCleaner.write_report(out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images=True, typ="refined")
#         taxCleaner.write_report(out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images=True, typ="corrected")

#     if not integrate:
#         taxCleaner.write_new_taxonomy(taxout, taxdiffs, taxboth)

#     info("Done!\n")


# def tax_imputation(proj, tax, integrate=False, mtdt=None, inps=None):
#     # inp_fol = cf_input if cf_input else "input/"
#     # inp_fol += proj+"/"
#     # out_fol = cf_output if cf_output else "output/"
#     # out_fol += proj+"/"
#     taxin = ppa_tax if integrate else inp_fol+tax
#     # taxout = None if integrate else out_fol+"new."+tax
#     # taxdiffs = None if integrate else out_fol+"diff."+tax
#     # taxboth = None if integrate else out_fol+"comp."+tax
#     intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
#     operation = "imputation" if integrate else "curation"

#     if not os.path.isfile(intree):
#         intree = out_fol+proj+".tree."+ ("int." if integrate else "")+"nwk"

#     info("Trying to impute taxonomic labels for taxa newly integrated into the tree... ")
#     taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, inps=inps)
#     taxCleaner.infer_tlabels()
#     info("Done!\n")
#     info("Writing taxonomic "+operation+" outputs... ")
#     taxCleaner.write_output(proj, images=False, imp_only=True)
#     info("Done!\n")


# def tax_curation_test(proj, tax, nerrorrs=10, error_type='wrong', taxl='s', tmin=3, tex=None, name='err.txt', integrate=False, descr=None):
#     err_tax = out_fol+name
#     intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
#     taxCleaner = taxc.TaxCleaner(inp_fol+tax, intree)
#     taxCleaner.introduce_errors(nerrorrs, err_tax, error_type, taxl, tmin, tex)
#     taxCleaner = taxc.TaxCleaner(err_tax, intree)
#     info("Starting taxonomic curation: detecting suspicious taxonomic labels... ")
#     taxCleaner.remove_tlabels()
#     info("Done!\n")
#     info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy... ")
#     taxCleaner.infer_tlabels()
#     info("Done!\n")
#     info("Writing taxonomic curation outputs and evaluations... ")
#     taxCleaner.write_output("log.txt", proj, outname=name) # , images=True )
#     taxCleaner.evaluate(proj, inp_fol+tax, err_tax, name, descr=descr)
#     info("Done!\n")


# def notinclusters(xx, yy, clusters):
#     reverse = True if xx > yy else False

#     for x, y in clusters:
#         if reverse and (x > y):
#             if ((x <= xx) and (xx <= y)) and ((x <= yy) and (yy <= y)):
#                 return False
#         elif not reverse and (x < y):
#             if ((y <= xx) and (xx <= x)) and ((y <= yy) and (yy <= x)):
#                 return False

#     return True


# def longest_not_overlapping(points):
#     spoints = sorted(points, key=lambda x: abs(x[1]-x[0]), reverse=True)
#     lno = [spoints[0]]

#     for a, b in spoints[1:]:
#         to_add = True

#         for x, y in lno:
#             if not (((a < x) and (b < x) and (a < y) and (b < y)) or
#                     ((a > x) and (b > x) and (a > y) and (b > y))):
#                 to_add = False
#                 break

#         if to_add:
#             lno.append((a, b))

#     return lno


# def find_clusters(reprs, points, f1, f2, reverse=True):
#     clusters = set()
#     changes = True

#     while changes:
#         changes = False
#         tmp_clusters = set()

#         for x, y in set(reprs).difference(clusters):
#             xyrev = True if x > y else False

#             for a, b in points:

#                 if reverse and xyrev:
#                     a, b = b, a

#                 if ((x != a) and (y != b)) and notinclusters(a, b, tmp_clusters):
#                     if ((a <= x) and (b >= x) and (b <= y)) or ((a >= x) and (a <= y) and (b >= y)):
#                         tmp_clusters = tmp_clusters.union(set([(f1(x, a), f2(y, b))]))

#             if tmp_clusters:
#                 break
#             else:
#                 if (x, y) not in clusters:
#                     clusters = clusters.union(set([(x, y)]))

#         if tmp_clusters:
#             changes = True
#             to_remove = set([(x, y)])

#             for a, b in reprs:
#                 for x, y in tmp_clusters:
#                     if (a >= x and a <= y and b >= x and b <= y) and ((a, b) not in to_remove):
#                         to_remove = to_remove.union(set([(a, b)]))

#             reprs = list(set(reprs).difference(to_remove))
#             reprs.append((f1([a for a, _ in tmp_clusters]), f2([b for _, b in tmp_clusters])))
#             tmp_clusters = set()

#     return clusters


# def fake_proteome_exe(f):
#     key = f[f.rfind('/')+1:f.rfind('.')]

#     if not key in fin_dict:
#         return key

#     mappings = {}

#     with open(f) as hf:
#         for r in hf:
#             row = r.strip().split('\t')
#             parsed_row = [row[0], float(row[2])]
#             parsed_row += [int(e) for e in row[3:10]]
#             parsed_row += [float(e) for e in row[10:]]

#             if row[1] in mappings:
#                 mappings[row[1]].append(parsed_row)
#             else:
#                 mappings[row[1]] = [parsed_row]

#     clusters = {}

#     for c, d in mappings.iteritems():
#         clusters[c] = set()
#         yerev = []
#         norev = []

#         for dd in d:
#             if dd[7] > dd[8]:
#                 yerev.append((dd[7], dd[8]))
#             else:
#                 norev.append((dd[7], dd[8]))

#         if norev:
#             lno = longest_not_overlapping(norev) # find longest and not-overlapping segments
#             clusters[c] |= find_clusters(lno, norev, min, max, reverse=False) # search for clusters

#         if yerev:
#             lno = longest_not_overlapping(yerev) # find longest and not-overlapping segments
#             clusters[c] |= find_clusters(lno, yerev, max, min, reverse=True) # search for clusters

#     # open input file and extract and translate the clusters
#     fin = fin_dict[key]
#     ff = None

#     if os.path.isfile(inp_fol+fin+'.fna'):
#         ff = open(inp_fol+fin+'.fna')
#         aaa = inp_fol+fin+'.fna'
#     elif os.path.isfile(inp_fol+fin+'.fna.bz2'):
#         ff = BZ2File(inp_fol+fin+'.fna.bz2')
#         aaa = inp_fol+fin+'.fna.bz2'
#     else:
#         info('File: '+inp_fol+fin+'.fna(.bz2) not found!\n')

#     proteome = []
#     p2t = dict()

#     if ff:
#         for record in SeqIO.parse(ff, 'fasta'):
#             for contig, cpoints in clusters.iteritems():
#                 if record.id in contig:
#                     for s, e in cpoints:
#                         reverse = False

#                         if s > e:
#                             s, e = e, s
#                             reverse = True

#                     sequence1 = Seq(str(record.seq)[s:e])

#                     if s-1 < 0:
#                         sequence2 = Seq('N' + str(record.seq)[s:e])
#                     else:
#                         sequence2 = Seq(str(record.seq)[s-1:e])

#                     if s-2 < 0:
#                         if s-1 < 0:
#                             sequence3 = Seq('NN' + str(record.seq)[s:e])
#                         else:
#                             sequence3 = Seq('N' + str(record.seq)[s-1:e])
#                     else:
#                         sequence3 = Seq(str(record.seq)[s-2:e])

#                     # if sequence no div by 3, add Ns
#                     while (len(sequence1) % 3) != 0:
#                         sequence1 += Seq('N')

#                     while (len(sequence2) % 3) != 0:
#                         sequence2 += Seq('N')

#                     while (len(sequence3) % 3) != 0:
#                         sequence3 += Seq('N')

#                     if reverse:
#                         sequence1 = sequence1.reverse_complement()
#                         sequence2 = sequence2.reverse_complement()
#                         sequence3 = sequence3.reverse_complement()

#                     rev = ':c' if reverse else ':' # reverse or not

#                     seqid1 = key+'_'+record.id+rev+str(s)+'-'+str(e)+'_s1'
#                     seqid2 = key+'_'+record.id+rev+str(s-1)+'-'+str(e)+'_s2'
#                     seqid3 = key+'_'+record.id+rev+str(s-2)+'-'+str(e)+'_s3'

#                     aminoacids1 = Seq.translate(sequence1)
#                     aminoacids2 = Seq.translate(sequence2)
#                     aminoacids3 = Seq.translate(sequence3)

#                     proteome.append(SeqRecord(aminoacids1, id=seqid1, description=''))
#                     proteome.append(SeqRecord(aminoacids2, id=seqid2, description=''))
#                     proteome.append(SeqRecord(aminoacids3, id=seqid3, description=''))

#                     # save map from seqid to genome
#                     p2t[seqid1] = key
#                     p2t[seqid2] = key
#                     p2t[seqid3] = key

#         ff.close()

#         # write output file
#         if proteome:
#             with open(dat_fol+ors2prots[:ors2prots.rfind('.')]+'/'+key+'.txt', 'w') as outf:
#                 outf.write('\t'.join([key] + [seq.id for seq in proteome]) + '\n')

#             with open(loc_inp+key+'.faa', 'w') as ff:
#                 SeqIO.write(proteome, ff, 'fasta')

#             info(loc_inp+key+'.faa generated (from '+aaa+' and '+f+')\n')

#         # write the partial mapping of proteins into genomes
#         if p2t:
#             with open(dat_fol+p2t_map[:p2t_map.rfind('.')]+'/'+key+'.txt', 'w') as f:
#                 for k, v in p2t.iteritems():
#                     f.write(str(k)+'\t'+str(v)+'\n')

#     if proteome and p2t:
#         return key

#     return


# def fake_proteome(proj, fna_out, faa_cleaning, nproc):
#     loc_inp = dat_cln_fol if faa_cleaning else inp_fol
#     b6o_files = iglob(dat_fol+'tblastn/*.b6o')
#     fake_proteomes = None
#     done = []

#     # create necessary directories if they do not exists
#     if not os.path.isdir(dat_fol+ors2prots[:ors2prots.rfind('.')]):
#         os.mkdir(dat_fol+ors2prots[:ors2prots.rfind('.')])

#     if not os.path.isdir(dat_fol+p2t_map[:p2t_map.rfind('.')]):
#         os.mkdir(dat_fol+p2t_map[:p2t_map.rfind('.')])

#     if not os.path.isdir(loc_inp):
#         os.mkdir(loc_inp)
#     else: # check already computed result files
#         for fi, fo in fna_out:
#             if os.path.isfile(loc_inp+fo+'.faa'):
#                 done.append((fi,fo))

#     fin_dict = dict([(fo, fi) for fi, fo in list(set(fna_out) - set(done))])

#     if not fin_dict:
#         info('All fake proteomes already computed!\n')
#         if os.path.isfile(dat_fol+ors2prots):
#             return [row.strip().split('\t')[0] for row in open(dat_fol+ors2prots)]
#         else:
#             return None

#     info('Reading tblastn output .b6o files ('+str(len(fin_dict))+')\n')
#     terminating = mp.Event()
#     pool = mp.Pool(initializer=initu,
#                    initargs=(terminating, loc_inp, dat_fol, fin_dict, fake_proteomes, ),
#                    processes=nproc)

#     try:
#         fake_proteomes = pool.map(fake_proteome_exe, b6o_files)
#         pool.close()
#     except:
#         pool.terminate()
#         pool.join()
#         error("Quitting fake_proteome_exe()")
#         return

#     pool.join()

#     # merge ors2prots
#     with open(dat_fol+ors2prots, 'w') as f:
#         for i in iglob(dat_fol+ors2prots[:ors2prots.rfind('.')]+'/*'):
#             with open(i) as g:
#                 f.write(g.read())

#     # merge p2t_map
#     with open(dat_fol+p2t_map, 'w') as f:
#         for i in iglob(dat_fol+p2t_map[:p2t_map.rfind('.')]+'/*'):
#             with open(i) as g:
#                 f.write(g.read())

#     return fake_proteomes


# def uniq_filenames(faa_in, fna_in):
#     new_fna_in = []

#     for fna in fna_in:
#         if fna in faa_in:
#             hashh = hashlib.sha1()
#             hashh.update(fna)
#             new_fna_in.append((fna, fna+'_'+str(hashh.hexdigest()[:4])))
#         else:
#             new_fna_in.append((fna, fna))

#     return new_fna_in


if __name__ == '__main__':
    args = read_params()
    check_args(args, args.verbose)
    clean_all(loc_dat, args.verbose) if args.clean_all else None
    clean_project(args.data_folder, args.output_folder, args.verbose) if args.clean else None

    if not args.meta: # standard phylogeny reconstruction
        configs = read_configs(args.configs_folder, args.verbose)
        check_configs(configs, args.verbose)
        check_dependencies(configs, args.nproc, args.verbose)
        database = init_database(args.database, args.databases_folder, configs['markers_db'], args.verbose)
        input_fna, input_fna_compressed = load_input_files(args.input_folder, args.genome_extension, args.verbose)

        if input_fna:
            gene_markers_identification(configs, 'dna_map', input_fna, args.data_folder, args.database, database, args.nproc, args.verbose)
            gene_markers_extraction(args.data_folder+'dna_map/', args.data_folder+'dna_markers/', args.nproc, args.verbose)
            fake_proteomes(args.data_folder+'dna_markers/', args.data_folder+'fake_proteomes/', args.proteome_extension, args.nrpoc, args.verbose)

        input_faa, input_faa_compressed = load_input_files(args.input_folder, args.proteome_extension, args.verbose)
        input_faa_fake, input_faa_compressed_fake = load_input_files(args.data_folder+'fake_proteomes/', args.proteome_extension, args.verbose)
        # merge input_faa & input_faa_fake?
        input_faa = check_input_proteomes(input_faa, args.min_num_proteins, args.min_len_protein, args.verbose)

        if input_faa:
            gene_markers_identification(configs, 'aa_map', input_faa, args.data_folder+'aa_map/', args.database, database, args.nproc, args.verbose)
            gene_markers_extraction(args.data_folder+'aa_map/', args.data_folder+'aa_markers/', args.nproc, args.verbose)

        for i in input_fna_compressed+input_faa_compressed+input_faa_compressed_fake: # close all NamedTemporaryFile
            print('Removing temporary file {}'.format(i.name)) if args.verbose else None
            i.close()
    else: # metagenomic application
        pass

#     global cf_data
#     cf_data = pars['c_dat']
#     if cf_data:
#         cf_data += '/' if not cf_data.endswith('/') else ''

#     if pars['cleanall']:
#         if pars['taxonomic_analysis'] or pars['user_tree'] or pars['integrate']:
#             error("--cleanall is in conflict with -t, -u, and -i", exit=True)
#         else:
#             info("Cleaning \"data/\" and \"data/ppaalns/\" folders... ")
#             clean_all(cf_data if cf_data else 'data/')
#             info("Done!\n")

#         sys.exit(0)

#     projn = pars['inp']
#     if not projn:
#         error("Project name not provided.", exit=True)

#     global cf_input
#     cf_input = pars['c_in']
#     if cf_input:
#         cf_input += '/' if not cf_input.endswith('/') else ''

#     global cf_output
#     cf_output = pars['c_out']
#     if cf_output:
#         cf_output += '/' if not cf_output.endswith('/') else ''

#     global dat_fol
#     global dat_cln_fol
#     dat_fol = cf_data if cf_data else 'data/'
#     dat_fol += projn+'/'
#     dat_cln_fol = dat_fol+cleanfaa_fld

#     global inp_fol
#     global inp_cln_fol
#     inp_fol = cf_input if cf_input else 'input/'
#     inp_fol += projn+'/'
#     inp_cln_fol = inp_fol+cleanfaa_fld

#     global out_fol
#     out_fol = cf_output if cf_output else 'output/'
#     out_fol += projn+'/'

#     if pars['clean']:
#         if pars['taxonomic_analysis'] or pars['user_tree'] or pars['integrate']:
#             error("-c/--clean is in conflict with -t, -u, and -i", exit=True)
#         else:
#             info("Cleaning project \""+projn+"\"... ")
#             clean_project()
#             info("Done!\n")

#         sys.exit(0)

#     global cf_up
#     cf_up = pars['c_up']
#     if cf_up:
#         cf_up += '/' if not cf_up.endswith('/') else ''

#     dep_checks(pars['mafft'], pars['raxml'], pars['nproc'])

#     global cf_fna
#     global cf_udb
#     cf_fna, cf_udb = init()

#     info("Loading and checking input files: ")
#     t0 = time.time()
#     faa_in, fna_in, tax, rtax, mtdt = get_inputs(projn, pars)
#     check_inp(faa_in + [f for _, f in fna_in])
#     info(str(len(faa_in)+len(fna_in))+" input files loaded in "+str(int(time.time()-t0))+" s!\n")

#     if not tax and pars['taxonomic_analysis'] and not pars['integrate']:
#         error("No taxonomy file found for the taxonomic analysis", exit=True)

#     global min_uprots
#     min_uprots = pars['min_uprots']
#     t0 = time.time()
#     torm_sol = None

#     if fna_in:
#         try:
#             blast(dict(fna_in), pars['nproc'], projn, pars['blast_full'])
#         except:
#             error("Quitting PhyloPhlAn [blast error]", exit=True)

#         fake = fake_proteome(projn, fna_in, pars['faa_cleaning'], pars['nproc'])

#         if fake:
#             faa_in += [i for i in fake if i]
#             loc_inp = dat_cln_fol if pars['faa_cleaning'] else inp_fol
#             torm_sol = [loc_inp+i+'.faa' for i in fake if i]

#     if faa_in:
#         try:
#             faa2ppafaa(faa_in, pars['nproc'], projn, pars['faa_cleaning'])
#             gens2prots(faa_in, projn, pars['faa_cleaning'])
#         except:
#             error("Quitting PhyloPhlAn [usearch error]", exit=True)

#     t1 = time.time()
#     info("Mapping finished in "+str(int(t1-t0))+" secs.\n")

#     try:
#         faa2aln(pars['nproc'], projn, pars['integrate'], pars['mafft'])
#     except:
#         error("Quitting PhyloPhlAn [muscle error]", exit=True)

#     t2 = time.time()
#     info("Aligning finished in "+str(int(t2-t1))+" secs ("+str(int(t2-t0))+" total time).\n")
#     aln_merge(projn, pars['integrate'])
#     build_phylo_tree(projn, pars['integrate'], pars['nproc'], pars['raxml'])
#     t4 = time.time()
#     info("Tree building finished in "+str(int(t4-t2))+" secs ("+str(int(t4-t0))+" total time).\n")

#     # circlader(projn, pars['integrate'], tax)

#     # if pars['integrate'] and tax and not pars['taxonomic_analysis']:
#     #     tax_imputation( projn, tax, mtdt = mtdt, integrate = pars['integrate'], inps = inps )

#     if pars['integrate'] and pars['taxonomic_analysis']:
#         tax_imputation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)
#         sys.exit()

#     if pars['taxonomic_analysis']:
#         if pars['tax_test']:
#             nerrorrs, error_type, taxl, tmin, tex, name = pars['tax_test'].split(":")
#             tax_curation_test(projn, tax,
#                               nerrorrs=int(nerrorrs),
#                               error_type=error_type,
#                               taxl=taxl,
#                               tmin=int(tmin),
#                               tex=int(tex) if tex else None,
#                               name=name,
#                               descr=pars['tax_test'])
#         else:
#             tax_curation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)

#     if torm_sol:
#         info("Removing "+str(len(torm_sol))+" input files... ")
#         delete_fake_proteomes(torm_sol)
#         info("Done!\n")
