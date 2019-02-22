#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.04'
__date__ = '13 February 2019'


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


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

THRESHOLD = 0.05
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
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


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __repr__(self):
        return "({},{})".format(self.start, self.end)


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
    p.add_argument('-t', '--threshold', type=float, choices=[Range(0.0, 1.0)], default=THRESHOLD,
                   help="Used to determine the closest hits from the distance estimation")

    p.add_argument('--nproc', type=int, default=1, help="The number of CPUs to use")
    p.add_argument('--database_folder', type=str, default=DATABASE_FOLDER,
                   help="Path to the folder that contains the database file")

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
            error('"{}" folder not found, neither in "{}"'.format(args.database_folder, os.path.dirname(os.path.abspath(__file__))), exit=True)

    args.database = os.path.join(args.database_folder, args.database)
    args.mapping = os.path.join(args.database_folder, args.mapping)

    create_folder(args.output_prefix + '_sketches', verbose=args.verbose)
    create_folder(args.output_prefix + '_dists', verbose=args.verbose)

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def check_dependencies(verbose=False):
    if verbose:
        info('Checking "mash"\n')

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


def sketching_dists_filter(input_folder, input_extension, output_prefix, database, threshold, nproc=1, verbose=False):
    commands = []
    out_filter_fld = os.path.join(output_prefix + "_filter_{}".format(threshold))

    if not os.path.exists(out_filter_fld):
        create_folder(out_filter_fld, verbose=verbose)

    for i in glob.iglob(os.path.join(input_folder, '*' + input_extension)):
        out = os.path.splitext(os.path.basename(i))[0]
        out_sketch = os.path.join(output_prefix + "_sketches", out)
        out_dist = os.path.join(output_prefix + "_dists", out + ".tsv")
        out_filter = os.path.join(out_filter_fld, out + ".tsv")

        if os.path.isfile(out_sketch + ".msh") and os.path.isfile(out_dist) and os.path.isfile(out_filter):
            continue

        commands.append((i, out_sketch, out_dist, out_filter, database, threshold, verbose))

    if commands:
        terminating = mp.Event()

        with mp.Pool(initializer=initt, initargs=(terminating,), processes=nproc) as pool:
            try:
                [_ for _ in pool.imap_unordered(sketching_dists_filter_rec, commands, chunksize=1)]
            except Exception as e:
                error(str(e), init_new_line=True)
                error('sketching crashed', init_new_line=True, exit=True)
    else:
        info('Inputs already sketched, dist, and filtered\n')


def sketching_dists_filter_rec(x):
    if not terminating.is_set():
        try:
            inp_bin, out_sketch, out_dist, out_filter, db, thr, verbose = x

            # sketching
            if not os.path.isfile(out_sketch + ".msh"):
                t0 = time.time()
                info('Sketching "{}"\n'.format(inp_bin))
                cmd = ['mash', 'sketch', '-k', '21', '-s', '10000', '-o', out_sketch, inp_bin]

                try:
                    sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
                except Exception as e:
                    terminating.set()
                    remove_file(out_sketch + ".msh", verbose=verbose)
                    error(str(e), init_new_line=True)
                    error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                    raise

                t1 = time.time()
                info('"{}.msh" generated in {}s\n'.format(out_sketch, int(t1 - t0)))

#######

            # dist
            if not os.path.isfile(out_dist):
                t0 = time.time()
                info('Computing distance for "{}"\n'.format(out_sketch + ".msh"))
                cmd = ['mash', 'dist', db, out_sketch + ".msh"]

                try:
                    sb.check_call(cmd, stdout=open(out_dist, 'w'), stderr=sb.DEVNULL)
                except Exception as e:
                    terminating.set()
                    remove_file(out_dist, verbose=verbose)
                    error(str(e), init_new_line=True)
                    error('cannot execute command\n    {}'.format(' '.join(cmd)), init_new_line=True)
                    raise

                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out_dist, int(t1 - t0)))

            # filter
            if not os.path.isfile(out_filter):
                t0 = time.time()
                info('Filtering "{}"\n'.format(out_dist))
                with open(out_dist) as fin, open(out_filter, 'w') as fout:
                    for r in fin:
                        rc = r.strip().split('\t')  # database-ID, input-ID, mash-distance, p-value, #-shared-hashes

                        if float(rc[2]) <= thr:
                            fout.write(rc[0] + '\n')

                t1 = time.time()
                info('"{}" generated in {}s\n'.format(out_filter, int(t1 - t0)))

#######

        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('error while sketching disting filtering\n    {}'.format('\n    '.join([str(a) for a in x])), init_new_line=True)
            raise
    else:
        terminating.set()


def check_md5(tar_file, md5_file):
    md5_md5 = None
    md5_tar = None

    if os.path.isfile(md5_file):
        with open(md5_file) as f:
            for row in f:
                md5_md5 = row.strip().split(' ')[0]
    else:
        error('file "{}" not found!\n'.format(md5_file))

    # compute MD5 of .tar.bz2
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


def phylophlan_metagenomic():
    args = read_params()

    if args.verbose:
        info('phylophlan_metagenomic.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, verbose=args.verbose)
    check_dependencies(verbose=args.verbose)
    download(os.path.join(DOWNLOAD_URL, DATABASE_FILE + '.tar.bz2'), args.database + '.tar.bz2', verbose=args.verbose)
    download(os.path.join(DOWNLOAD_URL, DATABASE_FILE + '.md5'), args.database + '.md5', verbose=args.verbose)
    check_md5(args.database + '.tar.bz2', args.database + '.md5')

    # untar
    try:
        tarfile_handle = tarfile.open(tar_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except EnvironmentError:
        sys.stderr.write("Warning: Unable to extract {}.\n".format(tar_file))

    # uncompress sequences
    bz2_file = os.path.join(folder, "mpa_" + download_file_name + ".fna.bz2")
    fna_file = os.path.join(folder, "mpa_" + download_file_name + ".fna")

    if not os.path.isfile(fna_file):
        sys.stderr.write('\n\nDecompressing {} into {}\n'.format(bz2_file, fna_file))

        with open(fna_file, 'wb') as fna_h, bz2.BZ2File(bz2_file, 'rb') as bz2_h:
            for data in iter(lambda: bz2_h.read(100 * 1024), b''):
                fna_h.write(data)



#######

    download(os.path.join(DOWNLOAD_URL, MAPPING_FILE), args.mapping, verbose=args.verbose)
    sketching_dists_filter(args.input, args.input_extension, args.output_prefix, args.database, args.threshold,
                           nproc=args.nproc, verbose=args.verbose)

    # load mapping file
    sgb2taxa = dict([r.strip().split('\t') for r in bz2.open(args.mapping, 'rt')])
    taxa2bin = {}

    for b in glob.iglob(args.output_prefix + '_filter_{}/*.tsv'.format(args.threshold)):
        bc = b.replace(args.output_prefix + '_filter_{}/'.format(args.threshold), '').replace('.tsv', '')

        for taxa in set([sgb2taxa[r.strip()] for r in open(b)]):
            if taxa in taxa2bin:
                taxa2bin[taxa].append(bc)
            else:
                taxa2bin[taxa] = [bc]

    with open(args.output_prefix + '_filter_{}.tsv'.format(args.threshold), 'w') as f:
        info('#taxnomy\tnum_bins\tbins_list_tab_separated\n', init_new_line=True)  # header
        f.write('#taxnomy\tnum_bin\tbins_list_tab_separated\n')  # header

        for taxa, bins in taxa2bin.items():
            info('{}\t{}\t{}\n'.format(taxa, len(bins), '\t'.join(bins)))
            f.write('{}\t{}\t{}\n'.format(taxa, len(bins), '\t'.join(bins)))

    info('Results saved to "{}_filter_{}.tsv"\n'.format(args.output_prefix, args.threshold), init_new_line=True)


if __name__ == '__main__':
    t0 = time.time()
    phylophlan_metagenomic()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
