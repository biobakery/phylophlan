#!/usr/bin/env python3


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it), '
              'Mattia Bolzan (mattia.bolzan@unitn.it), '
              'Nicola Segata (nicola.segata@unitn.it)')
__version__ = '0.05'
__date__ = '5 July 2018'


import sys
import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap
from urllib.request import urlretrieve
from urllib.parse import urlencode
from urllib.request import Request
from urllib.request import urlopen
import time


DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
DOWNLOAD_URL = "https://bitbucket.org/nsegata/phylophlan/downloads/"
TAXA2CORE_FILE = "taxa2core_latest.txt"


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
                             'core proteins. The label must represents a species: '
                             '"--get_core_proteins s__Escherichia_coli"'))

    p.add_argument('-o', '--output', type=str, default=None,
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
    if args.input:
        if (not os.path.isdir(args.input)) and (not os.path.isfile(args.input)):
            error('input must be either a folder or a file', exit=True)

        if not args.db_name:
            error('-d/--db_name must be specified', exit=True)

        if os.path.isdir(args.input):
            if not args.input_extension:
                error('input is a folder, hence --input_extension must be specified', exit=True)

    if args.get_core_proteins:
        if args.get_core_proteins[:3] != 's__':  # not sure it's needed, but I don't think we are handling sitautions different than species right now! (Tue 19 Jun 2018)
            error('The taxonomic label provided "{}" does not starts with "s__"'.format(args.get_core_proteins), exit=True)

        if not args.output:
            args.output = args.get_core_proteins

            if verbose:
                info('Setting output folder "{}"'.format(args.output))
        elif not args.ouptut.endswith(args.get_core_proteins):
            args.output = os.path.join(args.output, args.get_core_proteins)

            if verbose:
                info('Setting output folder "{}"'.format(args.output))

        args.input = args.output
        args.input_extension = PROTEOME_EXTENSION
        args.output_extension = PROTEOME_EXTENSION

        if not args.db_name:
            args.db_name = args.get_core_proteins

            if verbose:
                info('Setting db_name to "{}"\n'.format(args.db_name))

    if not args.db_type:
        if not args.output_extension:
            error('either -t/--db_type or -x/--output_extension must be specified', exit=True)
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
            error("both -t/--db_type and -x/--output_extension were specified, don't know which one to use!",
                  exit=True)


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
            error('unable to download "{}"'.format(url))
    elif verbose:
        info('File "{}" present\n'.format(download_file))


def create_folder(output, verbose=False):
    if not os.path.exists(output):
        if verbose:
            info('Creating output folder "{}"\n'.format(output))

        os.mkdir(output, mode=0o775)
    elif verbose:
        info('Output "{}" folder present\n'.format(output))


def get_core_proteins(taxa2core_file, taxa_label, output, output_extension, verbose=False):
    core_proteins = {}
    url = None
    metadata = None

    # taxa2core format
    #
    #   # taxid   taxonomy      UniRefURL       core_list
    #   12345     k__|...|s__   http://..{}..   UP123;UP456;...

    for r in bz2.open(taxa2core_file, 'rt'):
        if r.startswith('#'):
            metadata = r.strip()
            continue

        r_clean = r.strip().split('\t')

        if taxa_label in r_clean[1].split('|'):
            url = r_clean[2]
            core_proteins[r_clean[1]] = r_clean[3].split(';')
        elif taxa_label in r_clean[1]:
            core_proteins[r_clean[1]] = None

    if not len(core_proteins):
        error('no entry found for "{}", please check the taxonomic label provided'.format(taxa_label), exit=True)
    elif len(core_proteins) > 1:
        error('{} entries found for "{}":\n{}    please check the taxonomic label provided'
              .format(len(core_proteins), taxa_label, '    - {}\n'.join(core_proteins.keys())), exit=True)

    retry2download = []

    for lbl, core_prots in core_proteins.items():
        if verbose:
            info('Downloading {} core proteins for {}\n'.format(len(core_prots), lbl))

        for core_prot in core_prots:
            local_prot = os.path.join(output, core_prot + output_extension)
            download(url.format(core_prot), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                retry2download.append(core_prot)

    not_mapped = []

    # try to re-map the ids in case "not mapped" store in not_mapped
    if retry2download:
        if verbose:
            info("Re-trying to download {} core proteins that just failed, please wait as it might take some time\n"
                 .format(len(retry2download)))
        idmapping_url = 'https://www.uniprot.org/uploadlists/'
        contact = "phylophlan@cibiocm.com"
        params = {'from': 'ACC+ID',
                  'to': 'NF90',
                  'format': 'tab',  # or 'list' for only converted clusters
                  'query': ' '.join(retry2download)}
        data = urlencode(params).encode('utf-8')
        request = Request(idmapping_url, data, headers={'User-Agent': 'Python {}'.format(contact)})

        try:
            response = urlopen(request, data)
            uniprotkb2uniref90 = [line.decode().split('\t')[:2] for line in response.readlines()]
        except Exception:
            error('unable convert UniProtKB ID to UniRef90 ID')

        for uniref90_id in (x[1].split('_')[-1] for x in uniprotkb2uniref90[1:]):
            local_prot = os.path.join(output, uniref90_id + output_extension)
            download(url.format(uniref90_id), local_prot, verbose=verbose)

            if not os.path.exists(local_prot):
                not_mapped.append(uniref90_id)

        if (len(uniprotkb2uniref90) - 1) != len(retry2download):
            request_id = uniprotkb2uniref90[0][0].split(':')[1]
            not_mapped_url = 'http://www.uniprot.org/mapping/{}.not'.format(request_id)

            try:
                not_mapped_request = urlopen(not_mapped_url)
                not_mapped_response = [x.decode().strip() for x in not_mapped_request.readlines()]
                not_mapped.extend(not_mapped_response[1:])
            except Exception:
                error('unable fetch not converted IDs')

    if not_mapped:
        nd_out = os.path.join(output, taxa_label + '_core_proteins_not_mapped.txt')

        if verbose:
            info('There are {} core proteins that could not be downloaded, writing thier IDs to "{}"\n'
                 .format(len(not_mapped), nd_out))

        with open(nd_out, 'w') as f:
            f.write('\n'.join(not_mapped))


def create_database(db_name, inputt, input_ext, output, overwrite, verbose=False):
    seqs = []

    if os.path.exists(output) and (not overwrite):
        error('output file exists and --overwrite not specified', exit=True)

    if os.path.isdir(inputt):
        for marker in glob.iglob(os.path.join(inputt, '*' + input_ext + '*')):
            seqs += [SeqRecord(record.seq,
                               id='_'.join([db_name.replace('_', '-'), record.id.replace('_', '-'), str(count)]),
                               description='')
                     for count, record in enumerate(SeqIO.parse(bz2.open(marker, 'rt') if marker.endswith('.bz2')
                                                                else open(marker), "fasta"))]
    else:
        seqs = [SeqRecord(record.seq,
                          id='_'.join([db_name.replace('_', '-'), record.id.replace('_', '-'), str(count)]),
                          description='')
                for count, record in enumerate(SeqIO.parse(bz2.open(inputt, 'rt') if inputt.endswith('.bz2')
                                                           else open(inputt), "fasta"))]

    if not seqs:
        error('no sequences found, make sure the input folder/file provided is not empty', exit=True)

    if verbose:
        info('Writing output database "{}"\n'.format(output))

    with open(output, 'w') as f:
        SeqIO.write(seqs, f, "fasta")


if __name__ == '__main__':
    args = read_params()
    check_params(args, verbose=args.verbose)
    create_folder(args.output, verbose=args.verbose)

    if args.get_core_proteins:
        taxa2core_file_latest = None
        download(os.path.join(DOWNLOAD_URL, TAXA2CORE_FILE), TAXA2CORE_FILE, verbose=args.verbose)

        with open(TAXA2CORE_FILE) as f:
            for r in f:
                if not r.startswith('#'):
                    taxa2core_file_latest = r.strip()
                    break  # file should contains only one line, i.e., the name of the latest taxa2core file

        download(os.path.join(DOWNLOAD_URL, taxa2core_file_latest), taxa2core_file_latest, verbose=args.verbose)
        get_core_proteins(taxa2core_file_latest, args.get_core_proteins, args.output, args.output_extension,
                          verbose=args.verbose)

    create_database(args.db_name, args.input, args.input_extension,
                    os.path.join(args.output, args.db_name + args.output_extension),
                    args.overwrite, verbose=args.verbose)
