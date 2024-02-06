#!/usr/bin/env python


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Katarina Mladenovic (katarina.mladenovic@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Fabio Cumbo (fabio.cumbo@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Paolo Manghi (paolo.manghi@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '3.0.40'
__date__ = '5 February 2024'


import sys
import glob
import os
import argparse as ap
from urllib.request import urlretrieve
import time
import subprocess as sb
import multiprocessing as mp
import bz2
import hashlib
import numpy as np
import tarfile
import datetime
import itertools
import copy
import pandas as pd
import csv


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn {} requires Python 3, your current Python version is {}.{}.{}"
                    .format(__version__, sys.version_info[0], sys.version_info[1], sys.version_info[2]))

HOW_MANY = "8"
#DOWNLOAD_URL = "https://www.dropbox.com/s/xdqm836d2w22npb/phylophlan_metagenomic.txt?dl=1"
DOWNLOAD_URL = "http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan_SGB_databases.txt"
DATABASE_FOLDER = 'phylophlan_databases/'


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
    p = ap.ArgumentParser(description=("The phylophlan_assign_sgbs.py script assigns the SGB and taxonomy to a given set of input genomes. "
                                       "Outputs can be of three types: (1) for each input genomes returns the list of the closest "
                                       "-n/--how_many SGBs sorted by average Mash distance; (2) for each input genomes returns the "
                                       "closest SGB, GGB, FGB, and reference genomes; (3) returns a all vs. all matrix with all the "
                                       "pairwise mash distances"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str,
                   help="Input folder containing genomes and/or metagenome-assembled genomes (MAGs) to be indexed")
    p.add_argument('-o', '--output_prefix', type=str, default=None,
                   help=("Prefix used for the output folders: indexed bins, distance estimations. If not specified, "
                         "the input folder will be used"))
    p.add_argument('-d', '--database', type=str, default=None,
                   help="Database name, available options can be listed using the --database_list parameter")
    p.add_argument('--database_list', action='store_true', default=False,
                   help="List of all the available databases that can be specified with the -d/--database option")
    p.add_argument('--database_update', action='store_true', default=False, help="Update the databases file")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help=("Specify the extension of the input file(s) specified via -i/--input. If not specified will "
                         "try to infer it from the input files"))
    p.add_argument('-n', '--how_many', type=str, default=HOW_MANY,
                   help=('Specify the number of SGBs to report in the output; "all" is a special value to report all the SGBs;'
                         ' this param is not used when "--only_input" is specified'))
    p.add_argument('--nproc', type=int, default=1, help="The number of CPUs to use")
    p.add_argument('--database_folder', type=str, default=DATABASE_FOLDER,
                   help="Path to the folder that contains the database file")
    p.add_argument('--only_input', action='store_true', default=False,
                   help="If specified provides a distance matrix between only the input genomes provided")
    p.add_argument('--add_ggb_fgb', action='store_true', default=False,
                   help=("If specified adds GGB and FGB assignments. It will be adding a column for each "
                         " that reports the closest reference genome, -n/--how_many will be set to 1"))
    p.add_argument('--overwrite', action='store_true', default=False, help="If specified overwrite the output file if exists")
    p.add_argument('--citation', action='version',
                   version=('Asnicar, F., Thomas, A.M., Beghini, F. et al. '
                            'Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using PhyloPhlAn 3.0. '
                            'Nat Commun 11, 2500 (2020). '
                            'https://doi.org/10.1038/s41467-020-16366-7'),
                   help="Show citation")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_assign_sgbs.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan_assign_sgbs.py version and exit")

    return p.parse_args()


def database_list(databases_folder, update=False, exit=False, exit_value=0):
    sgbs_url = os.path.basename(DOWNLOAD_URL).replace('?dl=1', '')
    download(DOWNLOAD_URL, sgbs_url, overwrite=update, verbose=False)
    urls = set([r.strip().split('\t')[-1].replace('.md5', '').replace('.tar', '')
                for r in open(sgbs_url)
                if not r.startswith('#') and (len(r.split('\t')) == 2) and ('tutorial' not in r)])

    if not update:
        info('\nAvailable databases that can be specified with -d/--database:\n    ')
        info('\n    '.join(urls) + '\n', exit=exit, exit_value=exit_value)


def check_params(args, verbose=False):
    # database folder
    if os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), args.database_folder)):
        args.database_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.database_folder)

        if verbose:
            info('Setting --database_folder to "{}"\n'.format(args.database_folder))

    if args.database_update:
        database_list(args.database_folder, update=args.database_update)
        args.database_update = False

    if args.database_list:
        database_list(args.database_folder, update=args.database_update, exit=True)

    if not os.path.isdir(args.database_folder):
        create_folder(args.database_folder, verbose=verbose)

    if not args.only_input:
        if (not args.input) or (not args.database):
            error('both -i/--input and -d/--database must be specified', init_new_line=True)
            database_list(args.database_folder, update=args.database_update, exit=True)

        args.mapping = args.database

        if not args.mapping.endswith('.txt.bz2'):
            args.mapping += '.txt.bz2'

    if not os.path.isdir(args.input):
        error('"{}" folder not found, -i/--input must be a folder'.format(args.input), exit=True)

    if not args.input_extension:
        exts = set([os.path.splitext(i)[1] for i in glob.iglob(args.input + '/*') if os.path.splitext(i)[1]])

        if len(exts) > 1:
            error('Could not automatically infer the input extension (extensions detected: "{}"), please specify '
                  'using the -e/--input_extension param'.format('", "'.join(exts)), exit=True)

        args.input_extension = list(exts)[0]

        if verbose:
            info('Setting input extension to "{}"\n'.format(args.input_extension))

    if not args.input_extension.startswith('.'):
        args.input_extension = '.' + args.input_extension

    if not args.output_prefix:
        args.output_prefix = os.path.abspath(args.input)

        if verbose:
            info('Setting output prefix to "{}"\n'.format(args.output_prefix))

    if os.path.isdir(args.output_prefix):
        if args.output_prefix.endswith('/'):
           args.output_prefix = args.output_prefix[:-1]

        args.output_prefix = os.path.join(args.output_prefix, os.path.basename(args.output_prefix))

        if verbose:
            info('Output prefix is a folder, setting it to "{}"\n'.format(args.output_prefix))

    if args.how_many != 'all':
        try:
            how_many = int(args.how_many)
        except Exception as _:
            if verbose:
                info('Unrecognized value "{}", setting -n/--how_many to default value "{}"'.format(args.how_many, HOW_MANY))

            args.how_many = HOW_MANY

        args.how_many = how_many

    create_folder(args.output_prefix + '_sketches', verbose=args.verbose)
    create_folder(args.output_prefix + '_sketches/inputs', verbose=args.verbose)
    create_folder(args.output_prefix + '_dists', verbose=args.verbose)
    create_folder(args.output_prefix + '_prefiltering', verbose=True)
    create_folder(args.output_prefix + '_prefiltering/pref_dbs', verbose=True)

    if args.add_ggb_fgb:
        args.how_many = 1

        if verbose:
            info('--add_ggb_fgb specified, setting -n/--how_many to "{}"'.format(args.how_many))

    if verbose:
        info('Arguments: {}\n'.format(vars(args)), init_new_line=True)


def check_dependencies(verbose=False):
    if verbose:
        info('Checking "mash"\n', init_new_line=True)

    try:
        sb.check_call(['mash'], stdout=sb.DEVNULL, stderr=sb.DEVNULL)
    except Exception as e:
        error(str(e), init_new_line=True)
        error('mash is not installed or not present in the system path\n', init_new_line=True, exit=True)


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return (byte / 1048576)


class ReportHook():

    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
v        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, overwrite=False, verbose=False):
    """
    Download a file from a url
    """

    if (not os.path.isfile(download_file)) or overwrite:
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            error('unable to download "{}"'.format(url), exit=True)
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Folder "{}" already present\n'.format(output))


def remove_file(filename, path=None, verbose=False):
    to_remove = ''

    if os.path.isfile(filename):
        to_remove = filename
    elif os.path.isfile(os.path.join(path, filename)):
        to_remove = os.path.join(path, filename)

    if os.path.isfile(to_remove):
        if verbose:
            info('Removing "{}"\n'.format(to_remove))

        os.remove(to_remove)
    elif verbose:
        error('cannot remove "{}", file not found'.format(to_remove))


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def sketching_inputs_for_input_input_dist(input_folder, input_extension, output_prefix, nproc=1, verbose=False):
    commands = []

    for i in glob.iglob(os.path.join(input_folder, '*' + input_extension)):
        out = os.path.splitext(os.path.basename(i))[0]
        out_sketch = os.path.join(output_prefix + "_sketches/inputs", out)
        commands.append((i, out_sketch, verbose))

    if commands:
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(sketching_inputs_for_input_input_dist_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('sketching crashed', init_new_line=True, exit=True)
    else:
        info('No inputs found!\n')


def sketching_inputs_for_input_input_dist_rec(x):
    if not terminating.is_set():
        try:
            inp_bin, out_sketch, verbose = x

            if verbose:
                t0 = time.time()
                info('Sketching "{}"\n'.format(inp_bin))

            # sketch
            if not os.path.isfile(out_sketch + ".msh"):
                cmd = ['mash', 'sketch', '-k', '21', '-s', '10000', '-o', out_sketch, inp_bin]

                try:
                    sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
                except Exception as e:
                    terminating.set()
                    remove_file(out_sketch + ".msh", verbose=verbose)
                    error(str(e), init_new_line=True)
                    error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                    raise

            if verbose:
                t1 = time.time()
                info('Sketch for "{}" computed in {}s\n'.format(inp_bin, int(t1 - t0)))

        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while sketching_inputs_for_input_input_dist_rec\n    {}'
                  .format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def sketching(input_folder, input_extension, output_prefix, nproc=1, verbose=False):
    commands = []

    for i in glob.iglob(os.path.join(input_folder, '*' + input_extension)):
        out = os.path.splitext(os.path.basename(i))[0]
        out_sketch = os.path.join(output_prefix + "_sketches/inputs", out)
        commands.append((i, out_sketch, verbose))

    if commands:
        terminating = mp.Event()
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(sketching_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('sketching crashed', init_new_line=True, exit=True)
    else:
        info('No inputs found!\n')


def sketching_rec(x):
    if not terminating.is_set():
        try:
            inp_bin, out_sketch, verbose = x

            if verbose:
                t0 = time.time()
                info('Sketching "{}"\n'.format(inp_bin))

            if not os.path.isfile(out_sketch + ".msh"):
                cmd = ['mash', 'sketch', '-k', '21', '-s', '10000', '-o', out_sketch, inp_bin]

                try:
                    sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
                except Exception as e:
                    terminating.set()
                    remove_file(out_sketch + ".msh", verbose=verbose)
                    error(str(e), init_new_line=True)
                    error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                    raise

            if verbose:
                t1 = time.time()
                info('Sketch for "{}" computed in {}s\n'.format(inp_bin, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while sketching\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def pasting(output_prefix, prj_name, verbose=False):
    outf = output_prefix + "_sketches/" + prj_name + "_paste_{}"
    inpf = output_prefix + "_sketches/{}_inputs_list_{{}}.txt".format(prj_name)

    if verbose:
        t0 = time.time()
        info('Pasting inputs\n')

    # mash paste crashes with >70k
    istart = 0
    iend = 0
    step = 50000
    chunk = 1
    msh_idx = glob.glob(output_prefix + "_sketches/inputs/*.msh")

    while iend <= len(msh_idx):
        iend = min(iend + step, len(msh_idx))

        with open(inpf.format(chunk), 'w') as f:
            f.write('\n'.join(msh_idx[istart:iend]) + '\n')

        if iend == len(msh_idx):
            break

        istart = iend
        chunk += 1

    for inpc in glob.iglob(output_prefix + "_sketches/{}_inputs_list_*.txt".format(prj_name)):
        chunk = int(inpc[::-1].split('.')[1].split('_')[0][::-1])
        outc = outf.format(chunk)

        if os.path.isfile('{}.msh'.format(outc)):
            if verbose:
                info('"{}.msh" already exists\n'.format(outc))
                remove_file(inpc, output_prefix + "_sketches")
            continue

        if verbose:
            info('Pasting "{}.msh"\n'.format(outc))

        cmd = ['mash', 'paste', '-l', outc, inpc]

        try:
            sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
            remove_file(inpc, output_prefix + "_sketches")
        except Exception as e:
            error(str(e), init_new_line=True)
            error('cannot execute command\n {}'.format(' '.join(cmd)), init_new_line=True)
            raise

    if verbose:
        t1 = time.time()
        info('Inputs pasted in {}s\n'.format(int(t1 - t0)))


def prefiltering(output_prefix, prj_name, db_pref, nproc=1, verbose=True):
    """
    disting with prefiltering dbs ( dbs have to be created from representative genomes containing SGB ID in filename )
    """
    commands = []

    for inps in glob.glob(output_prefix + "_sketches/" + prj_name + "_paste_*.msh"):
        dist_file = os.path.join(output_prefix + "_prefiltering/pref_dbs", os.path.basename(inps).replace('.msh', '.tsv'))
        commands.append((inps, db_pref, dist_file, verbose, nproc))

    if commands:
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(prefiltering_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('prefiltering crashed', init_new_line=False, exit=False)
    else:
        info('Prefiltering mash dist already computed!\n')


def prefiltering_rec(x):
    if not terminating.is_set():
        try:
            msh_idx, db_pref, dist_file, verbose, nproc = x

            if not os.path.isfile(dist_file):
                fout = open(dist_file, 'w')

                if verbose:
                    t0 = time.time()
                    info('Disting "{}" with prefiltering database\n'.format(msh_idx))

                for db in os.listdir(db_pref):
                    cmd = ['mash', 'dist', '-p', str(nproc), os.path.join(db_pref, db), msh_idx]

                    try:
                        sb.check_call(cmd, stdout=fout, stderr=sb.DEVNULL, env=os.environ.copy().update({'OMP_NUM_THREADS': '1'}))
                    except Exception as e:
                        terminating.set()
                        fout.close()
                        remove_file(dist_file, verbose=verbose)
                        error(str(e), init_new_line=True)
                        error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                        raise

                if verbose:
                    t1 = time.time()
                    info('Prefiltering dist for "{}" computed in {}s\n'.format(db_pref, int(t1 - t0)))

                fout.close()

            elif verbose:
                    info('"{}" already present\n'.format(dist_file))

#            remove_file(msh_idx)

        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while disting\n    {} \n with prefiltering database'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def prefiltering_pasting(output_prefix, input_extension, chocophlan_list, database, verbose=True):
    """
    pasting inputs per sgb
    """
    outf = output_prefix + "_prefiltering/" + "{}_paste"
    inpf = output_prefix + "_prefiltering/{}_genomes_inputs_list.txt"

    for i in os.listdir(os.path.join(output_prefix + "_prefiltering/pref_dbs")):
        table = pd.read_csv(os.path.join(output_prefix + "_prefiltering/pref_dbs"+'/'+i), sep='\t',names=['sgb','mag','dist','pvalue','hits'])
        mags_full = list(table['mag'].unique())

        table=table[table['dist'] < 0.2]
        table['sgb'] = table['sgb'].map(lambda x: x.split("__")[0].replace('genomes/SGB', "") )
        table['mag'] = table['mag'].map(lambda x: x.split('/')[-1].rsplit('.', 1)[0] + '.msh')
        input_mags = list(table['mag'].unique())
        diff = list(filter(lambda x: x.split('/')[-1].rsplit('.', 1)[0] + '.msh' not in input_mags, mags_full))

	# filter out mags who's closest SGB is not in Chocophlan and put them in unassigned.txt
        for mag in input_mags:
            mag_tbl = table[table['mag']==mag].sort_values('dist')
            if mag_tbl.iloc[0]['sgb'] not in chocophlan_list:
                table.drop(table[table['mag'] == mag].index, inplace=True)

        for sgb in chocophlan_list:
            genomes = [output_prefix+'_sketches/inputs/'+ g for g in table['mag'][table['sgb']==sgb].values]
            input_mags = [x for x in input_mags if output_prefix+'_sketches/inputs/' + x not in genomes]
            outs = outf.format(sgb)
            inpg = inpf.format(sgb)

            with open(inpg, 'w') as f:
                f.write('\n'.join(genomes))

            if os.path.isfile('{}.msh'.format(outs)):
                if verbose:
                    info('"{}.msh" already exists\n'.format(outs))
                    remove_file(inpg, output_prefix + "_prefiltering")
                continue

            if verbose:
                info('Pasting "{}.msh"\n'.format(outs))

            cmd = ['mash', 'paste', '-l', outs, inpg]

            try:
                sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
                remove_file(inpg, output_prefix + "_prefiltering")
            except Exception as e:
                error(str(e), init_new_line=True)
                error('cannot execute command\n {}'.format(' '.join(cmd)), init_new_line=True)
                raise

        with open('unassigned.txt', 'w') as f:
            for mag in input_mags:
                f.write(f"{mag.strip('.msh')} is close to an SGB not present in the {database} database\n")
            for mag in diff:
                f.write(f"{mag.split('/')[-1].rsplit('.', 1)[0].strip('.msh')} is not close to any SGB present in the {database} database\n")

        if verbose:
            t1 = time.time()
            info('Inputs pasted in {}s\n'.format(int(t1 - t0)))


def disting(output_prefix, db, nproc=1, verbose=True): 
    commands = []

    for inps in glob.glob(output_prefix + "_prefiltering/" + "*_paste.msh"):
        sgb_msh_idx=os.path.join(db, '{}.msh')
        dist_file = os.path.join(output_prefix + "_dists", os.path.basename(inps).replace('_paste.msh', '.tsv'))
        commands.append((inps, sgb_msh_idx, dist_file, verbose))

    if commands:
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(disting_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('disting crashed', init_new_line=True, exit=False)
    else:
        info('Mash dist already computed!\n')


def disting_rec(x):
    if not terminating.is_set():
        try:
            msh_idx, sgb_msh_idx, dist_file, verbose = x
            sgb_msh_idx=sgb_msh_idx.format(msh_idx.split('/')[-1].split('_')[0].strip('SGB'))

            if not os.path.isfile(dist_file):
                fout = open(dist_file, 'a')

                if verbose:
                    t0 = time.time()
                    info('Disting "{}"\n'.format(sgb_msh_idx))

                cmd = ['mash', 'dist', '-p', '1', sgb_msh_idx, msh_idx]

                try:
                    sb.check_call(cmd, stdout=fout, stderr=sb.DEVNULL, env=os.environ.copy().update({'OMP_NUM_THREADS': '1'}))
                except Exception as e:
                    terminating.set()
                    fout.close()
                    remove_file(dist_file, verbose=verbose)
                    error(str(e), init_new_line=True)
                    error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                    raise

                if verbose:
                    t1 = time.time()
                    info('Dist for "{}" computed in {}s\n'.format(sgb_msh_idx, int(t1 - t0)))

                fout.close()

            elif verbose:
                info('"{}" already present\n'.format(dist_file))

            remove_file(msh_idx)

        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while disting\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def disting_input_vs_input(output_prefix, prj_name, output_file, nproc=1, verbose=False):
    cmd = []
    inps = sorted(glob.glob(output_prefix + "_sketches/" + prj_name + "_paste_*.msh"))
    inps_combs = itertools.combinations_with_replacement(inps, r=2)
    _, out_ext = os.path.splitext(output_file)
    out_f = os.path.join(output_prefix + "_dists", prj_name)

    for inp_a, inp_b in inps_combs:
        a = inp_a.split('_')[-1].split('.')[0]
        b = inp_b.split('_')[-1].split('.')[0]
        out_t = '{}_{}vs{}{}'.format(out_f, a, b, out_ext)

        if not os.path.isfile(out_t):
            cmd = ['mash', 'dist', '-t', '-p', str(nproc), inp_a, inp_b]

            try:
                sb.check_call(cmd, stdout=open(out_t, 'w'), stderr=sb.DEVNULL, env=os.environ.copy().update({'OMP_NUM_THREADS': str(nproc)}))
            except Exception as e:
                remove_file(out_t, verbose=verbose)
                error(str(e), init_new_line=True)
                error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                raise


def group_assignment(output_prefix, db, groups_map, nproc=1, verbose=True):
    """
    `groups_map` is a tab-separated file of the format:
       SGB_ID<tab>SGB_ID[<tab>SGB_ID]*

    the fisrt SGB ID is the representative SGB ID of the group
    the list of the row describes all members of that group (representative included)

    Example:
        1397<tab>1395
        17162<tab>17159<tab>17164<tab>17158

    this function will return a dictionary like the following:
        {'SGB1397': 'SGB1397_group',
         'SGB1395': 'SGB1397_group',
         'SGB17162': 'SGB17162_group',
         'SGB17159': 'SGB17162_group',
         'SGB17164': 'SGB17162_group',
         'SGB17158': 'SGB17162_group'}

    mapping each SGB to the group its belong to.
    """
    sgb_2_group = {}

    with bz2.open(groups_map, 'rt') as f:
        for row in f:
            sgb_repr = row.strip().split('\t')[0]

            for sgb in row.strip().split('\t'):
                sgb_2_group[sgb] = f'{sgb_repr}_group'

    return sgb_2_group


def check_md5(tar_file, md5_file, verbose=False):
    md5_md5 = None
    md5_tar = None

    if os.path.isfile(md5_file):
        with open(md5_file) as f:
            for row in f:
                md5_md5 = row.strip().split(' ')[0]
    else:
        error('file "{}" not found!\n'.format(md5_file))

    # compute MD5 of .tar
    if os.path.isfile(tar_file):
        hash_md5 = hashlib.md5()

        with open(tar_file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)

        md5_tar = hash_md5.hexdigest()[:32]
    else:
        error('file "{}" not found!\n'.format(tar_file))

    if (md5_tar is None) or (md5_md5 is None):
        error("MD5 checksums not found, something went wrong!", exit=True)

    # compare checksums
    if md5_tar != md5_md5:
        error("MD5 checksums do not correspond! Try removing the database files and rerun PhyloPhlAn "
              "to re-download them. If this happens again, please report this error.", exit=True)


def untar_and_decompress(tar_file, database, database_folder, nproc=1, verbose=False):
    # untar
    if not os.path.isdir(database):
        try:
            if verbose:
                info('Untar {} mash database\n'.format(tar_file))
            tarfile_handle = tarfile.open(tar_file)
            tarfile_handle.extractall(path=database_folder)
            tarfile_handle.close()
        except EnvironmentError:
            error('Warning: Unable to extract "{}".\n'.format(tar_file))
    elif verbose:
        info('Mash database already untarred\n')

    # uncompress mash indexes
    commands = [(os.path.join(database, f), os.path.join(database, f.replace('.bz2', '')), verbose)
                for f in os.listdir(database)
                if not os.path.isfile(os.path.join(database, f.replace('.bz2', '')))]

    if commands:
        if verbose:
            info('Decompressing {} mash indexes\n'.format(len(commands)))

        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(decompress_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('untar_and_decompress crashed', init_new_line=True, exit=True)
    elif verbose:
        info('Mash indexes already decompressed\n')


def decompress_rec(x):
    if not terminating.is_set():
        try:
            bz2_file, msh_file, verbose = x

            if verbose:
                info('Decompressing "{}"\n'.format(bz2_file))

            with open(msh_file, 'wb') as msh_h, bz2.open(bz2_file, 'rb') as bz2_h:
                for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                    msh_h.write(data)
            remove_file(bz2_file)
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while decompress_rec\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def merging(output_prefix, prj_name, output_file, verbose=False):
    out_f = os.path.join(output_prefix + "_dists", prj_name)
    _, out_ext = os.path.splitext(output_file)
    to_be_merged = glob.glob('{}_*vs*{}'.format(out_f, out_ext))

    if len(to_be_merged) == 1:
        os.rename(to_be_merged[0], output_file)
    else:
        if verbose:
            info('[w] this operation can take several hours and hundreds of Giga of RAM\n')

        # Map pair-wise distance matrices
        matrix_map = {}

        for filepath in to_be_merged:
            indexes = sorted([ int(i) for i in os.path.splitext(os.path.basename(filepath))[0].split("_")[-1].split("vs") ])

            if indexes[0] not in matrix_map:
                matrix_map[indexes[0]] = {}

            matrix_map[indexes[0]][indexes[1]] = filepath

            if indexes[1] not in matrix_map:
                matrix_map[indexes[1]] = {}

            matrix_map[indexes[1]][indexes[0]] = filepath

        # Concat pair-wise distance matrices
        matrix = None
        matrix_empty = True

        for position1 in sorted(matrix_map.keys()):
            row = None
            row_empty = True

            for position2 in sorted(matrix_map[position1].keys()):
                if verbose:
                    info('Loading "{}"\n'.format(matrix_map[position1][position2]))

                partial = pd.read_csv( matrix_map[position1][position2], sep="\t", index_col=0 )

                if position2 > position1:
                    partial = partial.T

                if row_empty:
                    row = partial
                    row_empty = False
                else:
                    row = pd.concat( [ row, partial ], axis=1 )

            if matrix_empty:
                matrix = row
                matrix_empty = False
            else:
                matrix = pd.concat( [ matrix, row ] )

        if verbose:
            info('Writing "{}" output file\n'.format(output_file))

        matrix.to_csv( output_file, sep="\t" )


def phylophlan_metagenomic():
    args = read_params()

    if args.verbose:
        info('phylophlan_assign_sgbs.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    check_dependencies(verbose=args.verbose)

    prj_name=os.path.basename(args.output_prefix)

    if not args.only_input:  # if mashing vs. the SGBs
        if (    not os.path.exists(os.path.join(args.database_folder, args.database)) or
                not os.path.exists(os.path.join(args.database_folder, args.mapping)) or
                not os.path.exists(os.path.join(args.database_folder, f'{args.database}_groups.tsv.bz2')) ):
               # not os.path.exists(os.path.join(args.database_folder, args.database + '.md5'))    ):
            sgbs_url = os.path.basename(DOWNLOAD_URL).replace('?dl=1', '')
            download(DOWNLOAD_URL, sgbs_url, verbose=args.verbose)
            urls = [tuple(r.strip().split('\t')) for r in open(sgbs_url)
                    if not r.startswith('#') and (args.database in r) and (len(r.split('\t')) == 2)]
            if len(urls) != 2:
                error('invalid number of URLs for "{}" in the downloaded file'.format(args.database),
                      exit=True)

            for url, filename in urls:
                download(url, os.path.join(args.database_folder, filename), verbose=args.verbose)

        args.mapping = os.path.join(args.database_folder, args.mapping)
        groups_map = os.path.join(args.database_folder, f'{args.database}_groups.tsv.bz2')
        args.db_pref = os.path.join(args.database_folder, args.database + '_pref')
        args.database = os.path.join(args.database_folder, args.database)  # keep this at the end otherwise args.database will be the full path

        if ( os.path.exists(os.path.join(args.database_folder, args.database + '.md5')) or not os.path.exists(os.path.join(args.database_folder, args.database)) ):
            check_md5(args.database + '.tar', args.database + '.md5', verbose=args.verbose)
        untar_and_decompress(args.database + '.tar', args.database, args.database_folder, nproc=args.nproc, verbose=args.verbose)

        chocophlan_list = [i.replace('.msh', '') for i in os.listdir(args.database)]
    else:  # mashing inputs against themselves
        sketching_inputs_for_input_input_dist(args.input, args.input_extension, args.output_prefix,
                                              nproc=args.nproc, verbose=args.verbose)
        args.database = args.output_prefix + '_sketches'

    sketching(args.input, args.input_extension, args.output_prefix, nproc=args.nproc, verbose=args.verbose)
    pasting(args.output_prefix, os.path.basename(args.output_prefix), verbose=args.verbose)
    prefiltering(args.output_prefix, prj_name, args.db_pref , args.nproc, verbose=True)
    prefiltering_pasting(args.output_prefix, args.input_extension, chocophlan_list, os.path.basename(args.database), verbose=True)

    output_file = args.output_prefix + ('.tsv' if not args.only_input else '_distmat.tsv')

    if os.path.isfile(output_file) and (not args.overwrite):
        timestamp = str(datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        output_file = output_file.replace(".tsv", "_" + timestamp + ".tsv")

    if not args.only_input:  # if mashing vs. the SGBs
        disting(args.output_prefix, args.database, nproc=args.nproc, verbose=args.verbose)
        sgb_2_group = group_assignment(args.output_prefix, args.database_folder, groups_map, nproc=args.nproc, verbose=args.verbose)

        # SGBs mapping file
        if args.verbose:
            info('Loading SGB mapping file\n')

        metadata_rows = []

        for r in bz2.open(args.mapping, 'rt'):
            if r.startswith('#'):
                metadata_rows.append(r.strip().split('\t'))

        sketches_folder = args.output_prefix + "_sketches/inputs"
        dists_folder = args.output_prefix + "_dists"
        mdidx = dict([(m, i) for i, m in enumerate(metadata_rows[-1][2:])])
        sgb_2_info = dict([(r.strip().split('\t')[1], r.strip().split('\t')[2:])
                           for r in bz2.open(args.mapping, 'rt')
                           if (r.strip().split('\t')[0] == 'SGB') and (not r.startswith('#'))])
        ggb_2_info = None
        fgb_2_info = None
        sgb_2_ggb = {}
        ggb_2_sgb = {}
        sgb_2_fgb = {}
        binn_2_sgb = dict([(b.replace('.msh', ''), dict()) for b in os.listdir(sketches_folder)])
        binn_2_ggb = None
        binn_2_fgb = None
        binn_2_refgen = None
        refgen_list = None

        if args.add_ggb_fgb:
            if args.verbose:
                info('Loading GGB mapping file\n')

            ggb_2_info = dict([(r.strip().split('\t')[1], r.strip().split('\t')[2:])
                               for r in bz2.open(args.mapping, 'rt')
                               if (r.strip().split('\t')[0] == 'GGB') and (not r.startswith('#'))])

            # SGB <--> GGB
            for ggb, ggb_info in ggb_2_info.items():
                ggb_2_sgb[ggb] = [sgb.replace('SGB', '') for sgb in ggb_info[mdidx['List of reconstructed genomes']].split(',')]

                for sgb in ggb_info[mdidx['List of reconstructed genomes']].split(','):
                    sgb_2_ggb[sgb.replace('SGB', '')] = ggb

            binn_2_ggb = copy.deepcopy(binn_2_sgb)

            if args.verbose:
                info('Loading FGB mapping file\n')

            fgb_2_info = dict([(r.strip().split('\t')[1], r.strip().split('\t')[2:])
                               for r in bz2.open(args.mapping, 'rt')
                               if (r.strip().split('\t')[0] == 'FGB') and (not r.startswith('#'))])

            # SGB --> FGB
            for fgb, fgb_info in fgb_2_info.items():
                for ggb in fgb_info[mdidx['List of reconstructed genomes']].split(','):
                    for sgb in ggb_2_sgb[ggb.replace('GGB', '')]:
                        sgb_2_fgb[sgb] = fgb

            binn_2_fgb = copy.deepcopy(binn_2_sgb)

        if args.add_ggb_fgb:
            if args.verbose:
                info('Loading reference genomes list\n')

            binn_2_refgen = dict([(binn, []) for binn in binn_2_sgb])
            refgen_list = set(itertools.chain.from_iterable([x[mdidx['List of reference genomes']].strip().split(',')
                                                             for x in sgb_2_info.values()
                                                             if x[mdidx['List of reference genomes']].strip() != '-']))

        if args.verbose:
            info('Loading mash dist files\n')

        for sgb in os.listdir(dists_folder):
            sgb_id = sgb.replace('.tsv', '')
            binn_2_dists = dict([(binn, []) for binn in binn_2_sgb])

            with open(os.path.join(dists_folder, sgb), 'rt') as f:
                for r in f:
                    rc = r.strip().split('\t')
                    binn = os.path.splitext(os.path.basename(rc[1]))[0]
                    binn_2_dists[binn].append(float(rc[2]))

                    if args.add_ggb_fgb:
                        sgb_member = os.path.splitext(os.path.basename(rc[0]))[0]

                        #if sgb_member in refgen_list:
                        binn_2_refgen[binn].append((sgb_member, float(rc[2])))

            for binn, dists in binn_2_dists.items():
                if dists:
                    binn_2_sgb[binn][sgb_id] = np.mean(dists)

        if args.add_ggb_fgb:
            for binn, dists in copy.deepcopy(binn_2_sgb).items():
                if dists == {}:
                    binn_2_sgb.pop(binn)
                    if args.add_ggb_fgb:
                        binn_2_ggb.pop(binn)
                        binn_2_fgb.pop(binn)
            for binn, dists in copy.deepcopy(binn_2_refgen).items():
                if dists == []:
                    binn_2_refgen.pop(binn)

        if args.add_ggb_fgb:
            if args.verbose:
                info('Computing GGB distances\n')

            for ggb_id, info_list in ggb_2_info.items():
                sgb_2_count = dict([(sgb_id, (int(sgb_2_info[sgb_id][mdidx['Number of reconstructed genomes']]) +
                                              int(sgb_2_info[sgb_id][mdidx['Number of reference genomes']])))
                                    for sgb_id in [x.replace('SGB', '') for x in info_list[mdidx['List of reconstructed genomes']].split(',')
                                    if info_list[mdidx['List of reconstructed genomes']].strip() != '-']])
                tot_sgbs = sum(sgb_2_count.values())

                for binn in binn_2_sgb: 
                    binn_2_ggb[binn][ggb_id] = sum([((sgb_sum / tot_sgbs) * binn_2_sgb[binn][sgb_id]) for sgb_id, sgb_sum in sgb_2_count.items()
                                                    if sgb_id in binn_2_sgb[binn].keys()])

            if args.verbose:
                info('Computing FGB distances\n')

            for fgb_id, info_list in fgb_2_info.items():
                ggb_2_count = dict([(ggb_id, (int(ggb_2_info[ggb_id][mdidx['Number of reconstructed genomes']]) +
                                              int(ggb_2_info[ggb_id][mdidx['Number of reference genomes']])))
                                    for ggb_id in [x.replace('GGB', '') for x in info_list[mdidx['List of reconstructed genomes']].split(',')
                                    if info_list[mdidx['List of reconstructed genomes']].strip() != '-']])
                tot_sgbs = sum(ggb_2_count.values())

                for binn in binn_2_ggb:
                    binn_2_fgb[binn][fgb_id] = sum([((ggb_sum / tot_sgbs) * binn_2_ggb[binn][ggb_id])
                                                    for ggb_id, ggb_sum in ggb_2_count.items()])

        if args.how_many == 'all':
            args.how_many = len(glob.glob(os.path.join(args.database, '*.msh')))

        with open(output_file, 'w') as f:
            if args.add_ggb_fgb:
                f.write('#last SGB id {}\n'.format(sorted(sgb_2_info, key=int, reverse=True)[0]))
                f.write('#last GGB id {}\n'.format(sorted(ggb_2_info, key=int, reverse=True)[0]))
                f.write('#last FGB id {}\n'.format(sorted(fgb_2_info, key=int, reverse=True)[0]))
                f.write('\t'.join(['#input_bin',
                                   '[u|k]_SGBid:taxa_level:taxonomy:avg_dist',
                                   '[u|k]_GGBid:taxa_level:taxonomy:avg_dist',
                                   '[u|k]_FGBid:taxa_level:taxonomy:avg_dist',
                                   'reference_genome:dist[|reference_genome:dist]']) + '\n')

                for binn in binn_2_sgb:
                    sgb_id, sgb_dist = sorted(binn_2_sgb[binn].items(), key=lambda x: x[1])[:args.how_many][0]

                    if sgb_dist <= 0.05:  # if SGB is assigned then output the GGB and FGB of the assigned SGB
                        ggb_id = sgb_2_ggb[sgb_id.split('_')[0]]
                        ggb_dist = binn_2_ggb[binn][ggb_id]
                        fgb_id = sgb_2_fgb[sgb_id.split('_')[0]]
                        fgb_dist = binn_2_fgb[binn][fgb_id]
                    else:  # otherwise report the closest GGB and FGB found
                        ggb_id, ggb_dist = sorted(binn_2_ggb[binn].items(), key=lambda x: x[1])[:args.how_many][0]
                        fgb_id, fgb_dist = sorted(binn_2_fgb[binn].items(), key=lambda x: x[1])[:args.how_many][0]

                    # closest ref.gens within 5% Mash distance of the closest ref.gen and no more than 100 in any case
                    refgen_sorted = list(sorted(binn_2_refgen[binn], key=lambda x: x[1]))
                    refgen_closest_thr = refgen_sorted[0][1] + (refgen_sorted[0][1] * 0.05)
                    refgens_closest = [i for i in refgen_sorted if i[1] <= refgen_closest_thr][:100]

                    # group_assignment
                    if sgb_id in sgb_2_group:
                        sgb_id = sgb_2_group[sgb_id]

                    f.write('\t'.join([binn,
                                       "{}_{}:{}:{}:{}".format(sgb_2_info[sgb_id.split('_')[0]][5], sgb_id, sgb_2_info[sgb_id.split('_')[0]][6],
                                                               sgb_2_info[sgb_id.split('_')[0]][7], sgb_dist),
                                       "{}_{}:{}:{}:{}".format(ggb_2_info[ggb_id][5], ggb_id, ggb_2_info[ggb_id][6],
                                                               ggb_2_info[ggb_id][7], ggb_dist),
                                       "{}_{}:{}:{}:{}".format(fgb_2_info[fgb_id][5], fgb_id, fgb_2_info[fgb_id][6],
                                                               fgb_2_info[fgb_id][7], fgb_dist),
                                       "|".join(["{}:{}".format(a[0], a[1]) for a in refgens_closest])]) + '\n')
            else:
                f.write('\t'.join(['#input_bin'] + ['[u|k]_[S|G|F]GBid:taxa_level:taxonomy:avg_dist'] * args.how_many) + '\n')

                # group_assignment
                for binn, sgb_dists in binn_2_sgb.items():
                    if sgb_dists:
                        for i in sorted(sgb_dists.items(), key=lambda x: x[1]):
                            sgb_id=i[0]
                            if sgb_id in sgb_2_group:
                                sgb_id=sgb_2_group[sgb_id]
                                if sgb_id not in sgb_dists.keys():
                                    sgb_dists[sgb_id]=i[1]
                                    sgb_dists.pop(i[0])
                                else:
                                    if sgb_dists[sgb_id] > i[1]:
                                        sgb_dists[sgb_id] = i[1]
                                    sgb_dists.pop(i[0])


                        f.write('\t'.join([binn] + ["{}_{}:{}:{}:{}".format(sgb_2_info[i[0].split('_')[0]][5],
                                                                            i[0],
                                                                            sgb_2_info[i[0].split('_')[0]][6],
                                                                            sgb_2_info[i[0].split('_')[0]][7],
                                                                            i[1])
                                                    for i in sorted(sgb_dists.items(), key=lambda x: x[1])[:args.how_many]]) + '\n')

            if os.path.isfile('unassigned.txt'):
                with open('unassigned.txt') as infile:
                    f.write(infile.read())

                remove_file('unassigned.txt', verbose=True)


    else:  # input vs. input mode
        disting_input_vs_input(args.output_prefix, os.path.basename(args.output_prefix), output_file,
                               nproc=args.nproc, verbose=args.verbose)
        merging(args.output_prefix, os.path.basename(args.output_prefix), output_file, verbose=args.verbose)

    info('Results saved to "{}"\n'.format(output_file))


if __name__ == '__main__':
    t0 = time.time()
    phylophlan_metagenomic()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)



