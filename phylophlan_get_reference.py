#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.02'
__date__ = '17 May 2018'


import sys
import bz2
import os
import argparse as ap
from urllib.request import urlretrieve
import time


DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
TAXA2GENOMES_FILE = "taxa2genomes_latest.txt"


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

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', '--get', type=str,
                       help=('Specify the taxonomic label for which download the set of reference proteomes. '
                             'The label must represents a valid taxonomic level or the special case "all"'))
    group.add_argument('-l', '--list_clades', action='store_true', default=False,
                       help='Print for all taxa the total number of species and reference genomes available')

    p.add_argument('-e', '--output_file_extension', type=str, default='.faa.gz',
                   help="Specify path to the extension of the output files")
    p.add_argument('-o', '--output', type=str,
                   help="Specify path to the output folder where to save the files, required when -g/--get is specified")
    p.add_argument('-n', '--num_ref', type=int, default=4,
                   help='Specify how many reference proteomes to download, where -1 stands for "all available"')
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")

    return p.parse_args()


def check_params(args, verbose=False):
    if args.list_clades:
        return

    if args.get[:3] not in ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 'all']:
        error('Taxnomic label provided "{}" is not in the correct format'.format(args.get), exit=True)

    if not args.output:
        error('-o/--output is required', exit=True)

    if os.path.exists(args.output):
        if not os.path.isdir(args.output):
            error('output param is not a directory', exit=True)

    if args.num_ref < 0:
        args.num_ref = None

    if not args.output_file_extension.startswith('.'):
        args.output_file_extension = '.' + args.output_file_extension

    if args.output_file_extension.endswith('.'):
        args.output_file_extension = args.output_file_extension[:-1]


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info("Output folder already present\n")


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

            urlretrieve(url, filename=download_file, reporthook=ReportHook().report)
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" already present!\n'.format(download_file))


def list_available_clades(taxa2proteomes_file, verbose=False):
    clades = {}
    metadata = None

    for r in bz2.open(taxa2proteomes_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        taxa = r.strip().split('\t')[1].split('|')
        num_ref_gen = len(r.strip().split('\t')[-1].split(';'))

        for i, c in enumerate(taxa):
            cl = '|'.join(taxa[:i+1])

            if cl in clades:
                num_spp, num_ref = clades[cl]
                clades[cl] = (num_spp + 1, num_ref + num_ref_gen)
            else:
                clades[cl] = (1, num_ref_gen)

    info('#taxa\tspecies\treference_genomes\n')
    for k, v in sorted(clades.items(), key=lambda x: x[0]):
        info('\t'.join([k, str(v[0]), str(v[1])]) + '\n')


def get_reference_genomes(taxa2genomes_file, taxa_label, num_ref, out_file_ext, output, verbose=False):
    core_genomes = {}
    metadata = None
    url = None

    # taxa2genomes format
    #
    #   # taxid   taxonomy      UniRefURL       genomes_list
    #   12345     k__|...|s__   http://..{}..   UP123;UP456;...

    for r in bz2.open(taxa2genomes_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        r_clean = r.strip().split('\t')

        if (taxa_label in r_clean[1].split('|')) or (taxa_label == 'all'):
            url = r_clean[2]
            core_genomes[r_clean[1]] = r_clean[3].split(';')[:num_ref]

    if taxa_label != 'all':
        if not len(core_genomes):
            error('no entry found for "{}", please check the taxonomic label provided'.format(taxa_label), exit=True)

    for lbl, genomes in core_genomes.items():
        if verbose:
            info('Downloading {} reference genomes for {}\n'.format(len(genomes), lbl))

        for genome in genomes:
            download(url.format(genomes), os.path.join(output, genomes + out_file_ext), verbose=verbose)


if __name__ == '__main__':
    taxa2genomes_file_latest = None
    args = read_params()
    check_params(args, verbose=args.verbose)
    download(os.path.join(DOWNLOAD_URL, TAXA2GENOMES_FILE), TAXA2GENOMES_FILE, verbose=args.verbose)

    with open(TAXA2GENOMES_FILE) as f:
        for r in f:
            if not r.startswith('#'):
                taxa2genomes_file_latest = r.strip()
                break  # file should contains only one line, i.e., the name of the latest taxa2genomes file

    download(os.path.join(DOWNLOAD_URL, taxa2genomes_file_latest), taxa2genomes_file_latest, verbose=args.verbose)

    if args.list_clades:
        list_available_clades(taxa2genomes_file_latest, verbose=args.verbose)
        sys.exit(0)

    create_folder(os.path.join(args.output), verbose=args.verbose)
    get_reference_genomes(taxa2genomes_file_latest, args.get, args.num_ref, args.output_file_extension,
                          args.output, verbose=args.verbose)
