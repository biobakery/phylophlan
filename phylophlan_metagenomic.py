#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.06'
__date__ = '5 March 2019'


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


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

HOW_MANY = "10"
DOWNLOAD_URL = ""
MAPPING_FILE = "db_sgbs_4930_sgb2taxa.tsv.bz2"
DATABASE_FILE = "db_sgbs_4930_k21_s10000"
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
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str, required=True,
                   help="Input folder containing the metagnomic bins to be indexed")
    p.add_argument('-o', '--output_prefix', type=str, default=None,
                   help=("Prefix used for the output folders: indexed bins, distance estimations. If not specified, "
                         "the input folder will be used"))
    p.add_argument('-d', '--database', type=str, default=DATABASE_FILE,
                   help="Specify the name of the database, if not found locally will be automatically downloaded")
    p.add_argument('-m', '--mapping', type=str, default=MAPPING_FILE,
                   help="Specify the name of the mapping file, if not found locally will be automatically downloaded")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help=("Specify the extension of the input file(s) specified via -i/--input. If not specified will "
                         "try to infer it from the input files"))
    p.add_argument('-n', '--how_many', type=str, default=HOW_MANY,
                   help='Specify the number of SGBs to report in the output; "all" is a special value to report all the SGBs')
    p.add_argument('--nproc', type=int, default=1, help="The number of CPUs to use")
    p.add_argument('--database_folder', type=str, default=DATABASE_FOLDER,
                   help="Path to the folder that contains the database file")
    p.add_argument('--only_input', action='store_true', default=False,
                   help="If specified")
    p.add_argument('--overwrite', action='store_true', default=False, help="If specified overwrites the output file if exists")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_metagenomic.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan_metagenomic.py version and exit")
    return p.parse_args()


def check_params(args, verbose=False):
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
            info('Setting prefix output folder to "{}"\n'.format(args.output_prefix))

    if not os.path.isdir(args.database_folder):
        if os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), args.database_folder)):
            args.database_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.database_folder)

            if verbose:
                info('Setting --database_folder to "{}"\n'.format(args.database_folder))
        else:
            error('"{}" folder not found, neither in "{}"'.format(args.database_folder, os.path.dirname(os.path.abspath(__file__))),
                  exit=True)

    if args.how_many != 'all':
        try:
            how_many = int(args.how_many)
        except Exception as e:
            if verbose:
                info('Unrecognized value "{}", setting -n/--how_many to default value "{}"'.format(args.how_many, HOW_MANY))

            args.how_many = HOW_MANY

        args.how_many = how_many

    args.database = os.path.join(args.database_folder, args.database)
    args.mapping = os.path.join(args.database_folder, args.mapping)

    create_folder(args.output_prefix + '_sketches', verbose=args.verbose)
    create_folder(args.output_prefix + '_dists', verbose=args.verbose)

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
        Print download progress message
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


def download(url, download_file, verbose=False):
    """
    Download a file from a url
    """

    if not os.path.isfile(download_file):
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
        except EnvironmentError:
            error('unable to download "{}"'.format(url), exit=True)
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Output "{}" folder present\n'.format(output))


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


def sketching_dists_filter(input_folder, input_extension, output_prefix, database, nproc=1, verbose=False):
    commands = []

    for i in glob.iglob(os.path.join(input_folder, '*' + input_extension)):
        out = os.path.splitext(os.path.basename(i))[0]
        out_sketch = os.path.join(output_prefix + "_sketches", out)
        out_dist = os.path.join(output_prefix + "_dists", out)
        commands.append((i, out_sketch, out_dist, database, verbose))

    if commands:
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(sketching_dists_filter_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('sketching crashed', init_new_line=True, exit=True)
    else:
        info('No inputs found!\n')


def sketching_dists_filter_rec(x):
    if not terminating.is_set():
        try:
            inp_bin, out_sketch, dist_fld, db, verbose = x

            if verbose:
                t0 = time.time()
                info('Analyzing "{}"\n'.format(inp_bin))

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

            # dist
            for sgb_msh_idx in glob.iglob(os.path.join(db, '*.msh')):
                dist_file = os.path.join(dist_fld, os.path.basename(sgb_msh_idx).replace('.msh', '.tsv'))

                if not os.path.isdir(dist_fld):
                    create_folder(dist_fld)

                if not os.path.isfile(dist_file):
                    cmd = ['mash', 'dist', sgb_msh_idx, out_sketch + ".msh"]

                    try:
                        sb.check_call(cmd, stdout=open(dist_file, 'w'), stderr=sb.DEVNULL)
                    except Exception as e:
                        terminating.set()
                        remove_file(dist_file, verbose=verbose)
                        error(str(e), init_new_line=True)
                        error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                        raise

            if verbose:
                t1 = time.time()
                info('Analysis for "{}" completed in {}s\n'.format(inp_bin, int(t1 - t0)))
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while sketching disting filtering\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


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
        error("MD5 checksums do not correspond! If this happens again, you should remove the database files and "
              "rerun MetaPhlAn2 so they are re-downloaded", exit=True)


def untar_and_decompress(tar_file, folder, nproc=1, verbose=False):
    # untar
    try:
        tarfile_handle = tarfile.open(tar_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except EnvironmentError:
        error('Warning: Unable to extract "{}".\n'.format(tar_file))

    # uncompress mash indexes
    commands = [(os.path.join(folder, f), os.path.join(folder, f.replace('.bz2', '')), verbose)
                for f in os.listdir(folder)
                if not os.path.isfile(os.path.join(folder, f.replace('.bz2', '')))]

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
    else:
        info('Mash indexes already decompressed\n')


def decompress_rec(x):
    if not terminating.is_set():
        try:
            bz2_file, msh_file, verbose = x

            with open(msh_file, 'wb') as msh_h, bz2.BZ2File(bz2_file, 'rb') as bz2_h:
                for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                    msh_h.write(data)
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while decompress_rec\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def phylophlan_metagenomic():
    args = read_params()

    if args.verbose:
        info('phylophlan_metagenomic.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    check_dependencies(verbose=args.verbose)

    # if mashing vs. the SGBs
    if not args.only_input:
        download(os.path.join(DOWNLOAD_URL, DATABASE_FILE + '.tar'), args.database + '.tar', verbose=args.verbose)
        download(os.path.join(DOWNLOAD_URL, DATABASE_FILE + '.md5'), args.database + '.md5', verbose=args.verbose)
        check_md5(args.database + '.tar', args.database + '.md5', verbose=args.verbose)
        untar_and_decompress(args.database + '.tar', args.database, nproc=args.nproc, verbose=args.verbose)
        download(os.path.join(DOWNLOAD_URL, MAPPING_FILE), args.mapping, verbose=args.verbose)
    else:  # mashing inputs against themselves
        pass

    if args.how_many == 'all':
        args.how_many = len(glob.glob(os.path.join(args.database, '*.msh')))

    sketching_dists_filter(args.input, args.input_extension, args.output_prefix, args.database, nproc=args.nproc, verbose=args.verbose)

    # SGBs mapping file
    sgb_2_info = dict([(r.strip().split('\t')[0], r.strip().split('\t')[1:]) for r in bz2.open(args.mapping, 'rt')])

    # output
    binn_2_sgb = {}
    dists_folder = args.output_prefix + "_dists"

    for binn in os.listdir(dists_folder):
        binn_2_sgb[binn] = []
        binn_folder = os.path.join(dists_folder, binn)

        for sgb in os.listdir(binn_folder):
            avg_dist = np.mean([float(r.strip().split('\t')[2]) for r in open(os.path.join(binn_folder, sgb))])
            binn_2_sgb[binn].append((sgb.replace('.tsv', ''), avg_dist))

    output_file = args.output_prefix + '.tsv'

    if os.path.isfile(output_file) and (not args.overwrite):
        timestamp = str(datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        output_file = args.output_prefix + '_' + timestamp + '.tsv'

    with open(output_file, 'w') as f:
        f.write('\t'.join(['#input_bin'] + ['[u|k]_SGBid(taxa_level):avg_dist'] * args.how_many) + '\n')

        for binn, sgb_dists in binn_2_sgb.items():
            f.write('\t'.join([binn] + ["{}SGB_{}({}):{}".format('u' if sgb_2_info[i[0]][2].upper() == 'YES' else 'k',
                                                                 i[0],
                                                                 sgb_2_info[i[0]][3],
                                                                 i[1])
                                        for i in sorted(sgb_dists, key=lambda x: x[1])[:args.how_many]]) + '\n')

    info('Results saved to "{}"'.format(output_file), init_new_line=True)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan_metagenomic()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
