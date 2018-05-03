#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.01'
__date__ = '27 April 2018'


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
TAXA2PROTEOMES_FILE = "taxa2proteomes.tsv.bz2"


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
    group.add_argument('-a', '--get_all', action='store_true', default=False,
                       help=("If specify -n/--num_ref reference proteomes will be "
                             "downloaded for all taxonomic labels"))
    group.add_argument('-g', '--get_reference_proteomes', type=str,
                       help=('Specify the taxonomic label for which download the set of '
                             'reference proteomes. The label can represents: '
                             'a species ("--get_reference_proteomes s__Escherichia_coli") or a '
                             'a genus ("--get_reference_proteomes g__Escherichia")'))  # DA VERIFICARE!!!

    p.add_argument('-e', '--output_file_extension', type=str, default='.faa.gz',
                   help="Specify path to the extension of the output files")
    p.add_argument('-o', '--output', type=str, required=True,
                   help="Specify path to the output folder where to save the files")
    p.add_argument('-n', '--num_ref', type=int, default=4,
                   help=('Specify how many reference proteomes to download, where -1 '
                         'stands for "all available"'))
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")

    return p.parse_args()


def check_params(args, verbose=False):
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
                info('Downloading "{}" as "{}"\n'.format(url, download_file))

            urlretrieve(url, filename=download_file, reporthook=ReportHook().report)
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" already present!\n'.format(download_file))


def get_reference_proteomes(taxa2proteomes_file, download_url, taxa_label, num_ref,
                            out_file_ext, output, verbose=False):
    download(download_url, taxa2proteomes_file, verbose=verbose)

    core_proteomes = {}
    metadata = None

    for r in bz2.open(taxa2proteomes_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        if taxa_label:
            if r.strip().split('\t')[1].endswith(taxa_label):
                core_proteomes[r.strip().split('\t')[1]] = [tuple(t.split(' ')) for t in r.strip().split('\t')[-1].split(';')[:num_ref]]
            elif taxa_label in r.strip().split('\t')[1]:
                core_proteomes[r.strip().split('\t')[1]] = None
        else:
            core_proteomes[r.strip().split('\t')[1]] = [tuple(t.split(' ')) for t in r.strip().split('\t')[-1].split(';')[:num_ref]]

    if taxa_label:
        if not len(core_proteomes):
            error('no entry found for "{}", please check the taxonomic label provided'
                  .format(taxa_label), exit=True)
        elif len(core_proteomes) > 1:
            error('{} entries found for "{}":\n{}    please check the taxonomic label provided'
                  .format(len(core_proteomes), taxa_label, '    - {}\n'.join(core_proteomes.keys())),
                  exit=True)

    for _, prot_urls in core_proteomes.items():
        for prot, url in prot_urls:
            download(url.format(prot), os.path.join(output, prot + out_file_ext), verbose=verbose)


if __name__ == '__main__':
    args = read_params()
    check_params(args, verbose=args.verbose)
    create_folder(os.path.join(args.output), verbose=args.verbose)
    get_reference_proteomes(TAXA2PROTEOMES_FILE, os.path.join(DOWNLOAD_URL, TAXA2PROTEOMES_FILE),
                            None if args.get_all else args.get_reference_proteomes,
                            args.num_ref, args.output_file_extension, args.output, verbose=args.verbose)
