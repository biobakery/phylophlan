#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.03'
__date__ = '27 April 2018'


import sys
import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap
from urllib.request import urlretrieve
import time


DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
TAXA2CORE_FILE = "taxa2core.tsv.bz2"


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
    group.add_argument('-i', '--input', type=str,
                       help=("Specify the path to either the folder containing the marker "
                             "files or the file of markers, in (multi-)fasta format"))
    group.add_argument('-g', '--get_core_proteins', type=str, default=None,
                       help=('Specify the taxonomic label for which download the set of '
                             'core proteins. The label can represents: '
                             'a species ("--get_core_proteins s__Escherichia_coli") or a '
                             'a genus ("--get_core_proteins g__Escherichia")'))

    p.add_argument('-o', '--output', type=str, required=True, default=None,
                   help="Specify path to the output folder where to save the database")
    p.add_argument('-d', '--db_name', type=str, help="Specify the name of the output database")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help="Specify the extension of the input file(s) specified via -i/--input")
    p.add_argument('-t', '--db_type', default=None, choices=DB_TYPE_CHOICES,
                   help=('Specify the type of the database, where "n" stands for '
                         'nucleotides and "a" for amino acids'))
    p.add_argument('-x', '--output_extension', type=str, default=None,
                   help="Set the database output extension")
    p.add_argument('--overwrite', action='store_true', default=False,
                   help="If specified and the output file exists it will be overwritten")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")

    return p.parse_args()


def check_params(args, verbose=False):
    if not os.path.isdir(args.output):
        error('output param is not a directory', exit=True)

    if args.input:
        if (not os.path.isdir(args.input)) and (not os.path.isfile(args.input)):
            error('input must be either a folder or a file', exit=True)

        if not args.db_name:
            error('-d/--db_name must be specified', exit=True)

        if os.path.isdir(args.input):
            if not args.input_extension:
                error('input is a folder, hence --input_extension must be specified', exit=True)

    if args.get_core_proteins:
        args.input = args.output
        args.input_extension = PROTEOME_EXTENSION
        args.output_extension = PROTEOME_EXTENSION

        if not args.db_name:
            args.db_name = args.get_core_proteins

            if verbose:
                info('Setting db_name to "{}"\n'.format(args.db_name))

    if not args.db_type:
        if not args.output_extension:
            error('either -t (--db_type) or -x (--output_extension) must be specified', exit=True)
        else:
            if not args.output_extension.startswith('.'):
                args.output_extension = '.' + args.output_extension

            if args.output_extension.endswith('.'):
                args.output_extension = args.output_extension[:-1]
    else:
        if not args.output_extension:
            if 'n' == args.db_type:
                args.output_extension = GENOME_EXTENSION
            if 'a' == args.db_type:
                args.output_extension = PROTEOME_EXTENSION

            if verbose:
                info('Setting output extension to "{}"'.format(args.output_extension))
        else:
            error("both -t (--db_type) and -x (--output_extension) were specified, don't "
                  "know which one to pick!", exit=True)


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
                info('Downloading "{}"" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
        except EnvironmentError:
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" already present!\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output database folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info("Output database folder already present\n")


def get_core_proteins(taxa2core_file, download_url, taxa_label, output, output_extension, verbose=False):
    download(download_url, taxa2core_file, verbose=verbose)
    core_proteins = {}
    metadata = None

    for r in bz2.open(taxa2core_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
        elif r.strip().split('\t')[1].endswith(taxa_label):
            core_proteins[r.strip().split('\t')[1]] = [tuple(t.split(' ')) for t in r.strip().split('\t')[-1].split(';')]
        elif taxa_label in r.strip().split('\t')[1]:
            core_proteins[r.strip().split('\t')[1]] = None

    if not len(core_proteins):
        error('no entry found for "{}", please check the taxonomic label provided'
              .format(taxa_label), exit=True)
    elif len(core_proteins) > 1:
        error('{} entries found for "{}":\n{}    please check the taxonomic label provided'
              .format(len(core_proteins), taxa_label, '    - {}\n'.join(core_proteins.keys())),
              exit=True)

    for lbl, prot_urls in core_proteins.items():
        if verbose:
            info('Downloading {} core proteins for {}\n'.format(len(prot_urls), lbl))

        for prot, url in prot_urls:
            download(url.format(prot), os.path.join(output, prot + output_extension), verbose=verbose)


def create_database(db_name, inputt, input_ext, output, overwrite, verbose=False):
    seqs = []

    if os.path.exists(output) and (not overwrite):
        error('output file exists and --overwrite not specified', exit=True)

    if os.path.isdir(inputt):
        for marker in glob.iglob(os.path.join(inputt, '*' + input_ext + '*')):
            seqs += [SeqRecord(record.seq,
                               id='_'.join([db_name.replace('_', '-'), record.id.replace('_', '-'), str(count)]),
                               description='')
                     for count, record in enumerate(SeqIO.parse(bz2.open(marker, 'rt')
                                                                if marker.endswith('.bz2')
                                                                else open(marker),
                                                                "fasta"))]
    else:
        seqs = [SeqRecord(record.seq,
                          id='_'.join([db_name.replace('_', '-'), record.id.replace('_', '-'), str(count)]),
                          description='')
                for count, record in enumerate(SeqIO.parse(bz2.open(inputt, 'rt')
                                                           if inputt.endswith('.bz2')
                                                           else open(inputt), "fasta"))]

    if not seqs:
        error('no sequences found, make sure the input folder/file provided is no empty',
              exit=True)

    if verbose:
        info('Writing output database "{}"\n'.format(output))

    with open(output, 'w') as f:
        SeqIO.write(seqs, f, "fasta")


if __name__ == '__main__':
    args = read_params()
    check_params(args, verbose=args.verbose)
    create_folder(args.output, verbose=args.verbose)

    if args.get_core_proteins:
        get_core_proteins(TAXA2CORE_FILE, os.path.join(TAXA2CORE_FILE, DOWNLOAD_URL),
                          args.get_core_proteins, args.output, args.output_extension, verbose=args.verbose)

    create_database(args.db_name, args.input, args.input_extension,
                    os.path.join(args.output, args.db_name + args.output_extension),
                    args.overwrite, verbose=args.verbose)
