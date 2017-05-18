#!/usr/bin/env python3


import bz2
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


database_name = 'amphora2'

for marker in glob.iglob('databases/amphora2/*.orig'):
    print(marker)
    seqid = 0
    out = []

    # with bz2.open(marker, 'rt') as g:
    with bz2.BZ2File(marker, 'rU') as g:
        for record in SeqIO.parse(g, "fasta"):
            out.append(SeqRecord(record.seq, id='_'.join([database_name, marker[marker.rfind('/')+1:marker.find('.')], str(seqid)]), description=''))
            seqid += 1

    with open(marker[:marker.find('.')]+'.faa', 'w') as f:
        SeqIO.write(out, f, "fasta")
