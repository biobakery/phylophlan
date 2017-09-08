#!/usr/bin/env python3


import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ############
# # AMPHORA2 #
# ############
# database_name = 'amphora2'

# for marker in glob.iglob('databases/amphora2/*.orig'):
#     print(marker)
#     seqid = 0
#     out = []

#     with bz2.BZ2File(marker, 'rU') as g:
#         for record in SeqIO.parse(g, "fasta"):
#             out.append(SeqRecord(record.seq, id='_'.join([database_name, marker[marker.rfind('/')+1:marker.find('.')], str(seqid)]), description=''))
#             seqid += 1

#     with open(marker[:marker.find('.')]+'.faa', 'w') as f:
#         SeqIO.write(out, f, "fasta")

# #############
# # C. parvum #
# #############
# database_name = 'cparvum'
# seqs = []
# seq_counter = 0

# for seq_record in SeqIO.parse('databases/cparvum/CryptoDB-33_CparvumIowaII_AnnotatedTranscripts.fasta', "fasta"):
#     seqs.append(SeqRecord(seq_record.seq, id='{}_{}_{}'.format(database_name, seq_record.id.replace('_', '-'), seq_counter), description=''))
#     seq_counter += 1

# with open('databases/cparvum/cparvum.fna', 'w') as f:
#     SeqIO.write(seqs, f, "fasta")

##############
# E. rectale #
##############
database_name = 'erectale'
seqs = []
seq_counter = 0

# for seq_record in SeqIO.parse('databases/erectale/core_genes_extracted.fa', "fasta"):
for seq_record in SeqIO.parse('databases/erectale/core_gene_reference_sequences.fa', "fasta"):
    seqs.append(SeqRecord(seq_record.seq, id='{}_{}_{}'.format(database_name, seq_record.id.replace('_', '-'), seq_counter), description=''))
    seq_counter += 1

with open('databases/erectale/erectale.fna', 'w') as f:
    SeqIO.write(seqs, f, "fasta")

