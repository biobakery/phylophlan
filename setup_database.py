#!/usr/bin/env python3


import sys
import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse as ap


DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'


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

    p.add_argument('-i', '--input', required=True, type=str,
                   help=("Specify the path to the folder containing the markers file or "
                         "file (multi-fasta) of markers"))
    p.add_argument('-o', '--output', required=False, type=str, default=None,
                   help=("Specify path to the output folder where to save the database, "
                         "if not specified the input path (specified via -i, --input) "
                         "will be used"))
    p.add_argument('-d', '--db_name', required=True, type=str,
                   help="Specify the name of the output database")
    p.add_argument('-e', '--input_extension', type=str, required=False, default=None,
                   help=("Specify the extension of the input file(s) specified via -i, "
                         "--input"))
    p.add_argument('-t', '--db_type', required=False, default=None, choices=DB_TYPE_CHOICES,
                   help=("Specify the type of the database, where 'n' stands for "
                         "nucleotides and 'a' for amino acids"))
    p.add_argument('-x', '--output_extension', type=str, default=None,
                   help="Set the database output extension")
    p.add_argument('--overwrite', action='store_true', default=False,
                   help="If specified and the output file exists it will be overwritten")
    p.add_argument('--verbose', action='store_true', default=False, help="Prints more stuff")

    return p.parse_args()


def check_params(args, verbose=False):
    if (not os.path.isdir(args.input)) and (not os.path.isfile(args.input)):
        error('input must be either a folder or a file', exit=True)

    if not args.output:
        args.output = os.path.dirname(args.input) if os.path.isfile(args.input) else args.input

        if verbose:
            info('Setting output to "{}"\n'.format(args.output))

    if os.path.isdir(args.input):
        if not args.input_extension:
            error('input is a folder, hence --input_extension must be specified',
                  exit=True)

    if not args.db_type:
        if not args.output_extension:
            error('either -t (--db_type) or -x (--output_extension) must be specified',
                  exit=True)
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


def create_folder(output, verbose=False):
    if not os.path.isdir(output):
        if verbose:
            info('Creating output database folder "{}"\n'.format(output))

        os.mkdir(output)
    elif verbose:
        info("Output database folder already present\n")


def create_database(db_name, inputt, input_ext, output, overwrite, verbose=False):
    seqs = []

    if os.path.exists(output) and (not overwrite):
        error('output file exists --overwrite not provided', exit=True)

    if os.path.isdir(inputt):
        for marker in glob.iglob(os.path.join(inputt, '*' + input_ext + '*')):
            seqs += [SeqRecord(record.seq,
                               id='_'.join([db_name, record.id.replace('_', '-'), str(count)]),
                               description='')
                     for count, record in enumerate(SeqIO.parse(bz2.open(marker, 'rt')
                                                                if marker.endswith('.bz2')
                                                                else open(marker),
                                                                "fasta"))]
    else:
        seqs = [SeqRecord(record.seq,
                          id='_'.join([db_name, record.id.replace('_', '-'), str(count)]),
                          description='')
                for count, record in enumerate(SeqIO.parse(bz2.open(inputt, 'rt')
                                                           if inputt.endswith('.bz2')
                                                           else open(inputt), "fasta"))]

    if not seqs:
        error('no sequences found, make sure the input folder/file provided is no empty',
              exit=True)

    if verbose:
        info('Writing output database "{}"'.format(output))

    with open(output) as f:
        SeqIO.write(seqs, f, "fasta")


if __name__ == '__main__':
    args = read_params()
    check_params(args, verbose=args.verbose)
    create_folder(os.path.join(args.db_name, args.output), verbose=args.verbose)
    create_database(args.db_name, args.input, args.input_extension,
                    os.path.join(args.output, args.db_name + args.output_extension),
                    args.overwrite, verbose=args.verbose)
