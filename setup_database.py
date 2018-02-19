#!/usr/bin/env python3

import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
#import sys
#import stat
import argparse as ap
#import configparser as cp

DB_TYPE_CHOICES = ['n', 'a']
GENOME_EXTENSION = '.fna'
PROTEOME_EXTENSION = '.faa'
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
DATABASES_FOLDER = os.path.join(SCRIPT_PATH, 'databases')

def read_params():
    p = ap.ArgumentParser(description="")

    p.add_argument('-i', '--input_folder', required=True,
                   help=("Specify the path to the input folder or file"))

    p.add_argument('-n', '--db_name', required=True,
                   help=("Specify the name of the database"))

    p.add_argument('-d', '--db_type', required=False, choices=DB_TYPE_CHOICES,
                   help=("Specify the type of the database, where 'n' stands for "
                     "nucleotides and 'a' for amino acids"))
    p.add_argument('-e', '--input_extension', type=str, required=True,
                   help=("Specify the extension of input files in the folder)"))

    p.add_argument('--genome_extension', type=str,
                    default=GENOME_EXTENSION,
                    help=("Set the extension for the genomes in your "
                             "inputs, default .fna"))
    p.add_argument('--proteome_extension', type=str,
                    default=PROTEOME_EXTENSION,
                    help=("Set the extension for the proteomes in your "
                             "inputs, default .faa"))
#    p.add_argument('-o', '--output_folder', required=False,
#                   help=("Specify path to the output folder where to save the database, default is input folder") )    

    p.add_argument('--overwrite', action='store_true', default=False,
                   help="If output file exists it will be overwritten")
    p.add_argument('--verbose', action='store_true', default=False,
                   help="Prints more stuff")

    return p.parse_args()

def create_folder(folder, create=False):
    if not os.path.exists(DATABASES_FOLDER):
        os.mkdir(DATABASES_FOLDER)
    else:
        pass
    db_folder = os.path.join(DATABASES_FOLDER, folder)
    if not os.path.isdir(db_folder):
        if create:
            os.mkdir(db_folder)
        else:
            pass
    else:
        print("Database folder is already present")
    return db_folder

def input_files(input_folder, extension):
    inputs = {}
    if os.path.isdir(input_folder):
        inputs = glob.iglob(os.path.join(input_folder, '*' + '.' + extension), recursive=False)
    elif os.path.isfile(input_folder):
        inputs = input_folder
    return inputs

def create_database():
    files = input_files(args.input_folder, args.input_extension)
    database_name = args.db_name
    if args.db_type:
        if 'n' in args.db_type:
            output_file = os.path.join(db_folder, args.db_name + GENOME_EXTENSION) 
        elif 'a' in args.db_type:
            output_file = os.path.join(db_folder, args.db_name + PROTEOME_EXTENSION)
    output_handle = open(output_file, 'w')
    if not os.path.isfile(args.input_folder):
        for marker in files:
            seqid = 0
            seqs = []
            with open(marker) as g:
                for record in SeqIO.parse(g, "fasta"):
    #                seqs.append(SeqRecord(record.seq, id='_'.join([database_name, marker[marker.rfind('/')+1:marker.find('.')], str(seqid)]), description=''))
                    seqs.append(SeqRecord(record.seq, id='_'.join([database_name, record.id.replace('_', '-'), str(seqid)]), description=''))
                    seqid += 1
                SeqIO.write(seqs, output_handle, "fasta")
    else:
        seqs = []
        seqid = 0
        for record in SeqIO.parse(files, "fasta"):
            seqs.append(SeqRecord(record.seq, id='{}_{}_{}'.format(database_name, record.id.replace('_', '-'), seqid), description=''))
            seqid += 1
        SeqIO.write(seqs, output_handle, "fasta")

#    if os.path.isfile(args.input_folder):
#        seqs = []
#        seq_counter = 0
#        for seq_record in SeqIO.parse(files, "fasta"):
#            seqs.append(SeqRecord(seq_record.seq, id='{}_{}_{}'.format(database_name, seq_record.id.replace('_', '-'), seq_counter), description=''))
#            seq_counter += 1
#        SeqIO.write(seqs, output_handle, "fasta")
#    elif os.path.isdir(args.input_folder):
#        for marker in files:
#            seqid = 0
#            seqs = []
#            with open(marker) as g:
#                for record in SeqIO.parse(g, "fasta"):
#    #                seqs.append(SeqRecord(record.seq, id='_'.join([database_name, marker[marker.rfind('/')+1:marker.find('.')], str(seqid)]), description=''))
#                    seqs.append(SeqRecord(record.seq, id='_'.join([database_name, record.id.replace('_', '-'), str(seqid)]), description=''))
#                    seqid += 1
#                SeqIO.write(seqs, output_handle, "fasta")
    output_handle.close()

if __name__ == '__main__':
    args = read_params()
    db_folder = create_folder(args.db_name, True)
    create_database()