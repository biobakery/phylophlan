#!/usr/bin/env python


__author__ = 'Nicola Segata (nsegata@hsph.harvard.edu)'
__version__ = '1.10'
__date__ = '17 Sep 2014'


import sys
import os
import tarfile
import argparse as ap
import subprocess as sb
import multiprocessing as mp
import collections
try:
    import cPickle as pickle
except:
    import pickle
import urllib2
from contextlib import closing
from glob import iglob
from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# sys.path.insert(0, 'taxcuration/')
import taxcuration as taxc
import shutil
import time
from bz2 import BZ2File
from tempfile import NamedTemporaryFile
from itertools import chain
import hashlib
try:
    from collections import Counter # works only with python >= 2.7
    collections_counter = True
except:
    collections_counter = False
# from random import randint
# import traceback


download = ""
ppa_fna = "data/ppa.seeds.faa"
ppa_fna_40 = "data/ppa.seeds.40.faa"
ppa_aln = "data/ppafull.aln.faa"
ppa_up2prots = "data/ppafull.up2prots.txt"
ppa_ors2prots = "data/ppafull.orgs2prots.txt"
ppa_tax = "data/ppafull.tax.txt"
ppa_alns = ("data/ppaalns/list.txt", "data/ppaalns/ppa.aln.tar.bz2")
ppa_alns_fol = "data/ppaalns/"
ppa_xml = "data/ppafull.xml"
ppa_udb = "data/ppa.udb"
up2prots = "up2prots.txt"
ors2prots = "orgs2prots.txt"
aln_tot = "aln.fna"
aln_int_tot = "aln.int.fna"
ups2faa_pkl = "ups2faa.pkl"
o_tree = "tree.nwk"
o_inttree = "tree.int.nwk"
p2t_map = "p2t_map.txt"
cleanfaa_fld = "clean_faa/"
old2newprots = "old2newprots.txt"
NOT_ENOUGH_MAPPINGS = "Not enough mappings"
few_maps = "few_maps.txt"
f_protein_ids_map = "protein_ids_map.pkl"


# PEP8 E731
# compressed = lambda x: x+".tar.bz2"
def compressed(x):
    return x+".tar.bz2"


# PEP8 E731
# download_compressed = lambda x: download+os.path.basename(x)+".tar.bz2"
def download_compressed(x):
    return download+os.path.basename(x)+".tar.bz2"


def info(s):
    sys.stdout.write(str(s))
    sys.stdout.flush()


def error(s, init_new_line=False, exit=False, exit_value=1):
    err = '\n' if init_new_line else ''
    err += '[e] '
    sys.stderr.write(err + str(s) + '\n')

    if exit:
        sys.stderr.write('Exiting...\n')

    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def dep_checks():
    for prog in ["FastTree", "usearch", "muscle", "tblastn"]:
        try:
            with open(os.devnull, 'w') as devnull:
                t = sb.Popen([prog], stdout=devnull, stderr=devnull)
                t.wait()
        except OSError:
            t.kill()
            error(prog+" not found or not in system path", exit=True)


def read_params(args):
    p = ap.ArgumentParser(description=
            "NAME AND VERSION:\n"
            "PhyloPhlAn version "+__version__+" ("+__date__+")\n\n"
            "AUTHORS:\n"
            "Nicola Segata (nsegata@hsph.harvard.edu) and Curtis Huttenhower (chuttenh@hsph.harvard.edu)\n\n"
            "DESCRIPTION\n"
            "PhyloPhlAn is a computational pipeline for reconstructing highly accurate and resolved \n"
            "phylogenetic trees based on whole-genome sequence information. The pipeline is scalable \n"
            "to thousands of genomes and uses the most conserved 400 proteins for extracting the \n"
            "phylogenetic signal.\n"
            "PhyloPhlAn also implements taxonomic curation, estimation, and insertion operations.\n\n",
             formatter_class=ap.RawTextHelpFormatter)
    arg = p.add_argument

    arg('inp', metavar='PROJECT NAME', type=str, default=None, nargs='?', help=
        "The basename of the project corresponding to the name of the input data folder inside \n"
        "input/. The input data consist of a collection of multifasta files (extensions allowed are:\n"
        ".faa and .fna, or compressed: .faa.bz2 and .fna.bz2) containing the proteins in each genome.\n"
        "If the project already exists, the already executed steps are not re-ran.\n"
        "The results will be stored in a folder with the project basename in output/\n"
        "Multiple project can be generated and they safetely coexists.")

    arg('-i','--integrate', action='store_true', help=
        "Integrate user genomes into the PhyloPhlAn tree \n")
    arg('-u','--user_tree', action='store_true', help=
        "Build a phylogenetic tree using user genomes only \n")
    arg('-t','--taxonomic_analysis', action='store_true', help=
        "Check taxonomic inconsistencies and refine/correct taxonomic labels\n")
    arg('--tax_test', type=str, default = None, help=
        "nerrors:type:taxl:tmin:tex:name (alpha version, experimental!)\n")

    arg('-c','--clean', action='store_true', help=
        "Clean the final and partial data produced for the specified project.\n"
        " (use --cleanall for removing general installation and database files)\n")
    arg('--cleanall', action='store_true', help=
        "Remove all instalation and database file leaving untouched the initial compressed data \n"
        "that is automatically extracted and formatted at the first pipeline run.\n"
        "Projects are not remove (specify a project and use -c for removing projects).\n")

    arg('--nproc', metavar="N", type=int, default=1, help=
        "The number of CPUs to use for parallelizing the blasting\n"
        "[default 1, i.e. no parallelism]\n")
    # decide which database to use with blast
    arg('--blast_full', action='store_true', default=False, help=
        "If specified, tells blast to use the full dataset of universal proteins\n"
        "[default False, i.e. the small dataset of universal proteins is used]\n")
    # decide if you want perform the faa cleaning of the .faa
    arg('--faa_cleaning', action='store_true', default=False, help=
        "When specified perform a cleaning on the number and the length of proteins, changing also\n"
        "the proteins id such that are unique among all the genomes."
        "If you believe that your genomes have unique ids you can skip this step.")
    # protein filtering params
    arg('--min_num_proteins', type=int, default=100, help=
        "This parameter is used when the --faa_cleaning is specified. When performing the cleaning,\n"
        "genomes with less than this number of proteins will be discarded.\n"
        "[default minimum number of proteins is 100]\n")
    arg('--min_len_protein', type=int, default=50, help=
     "This parameter is used when the --faa_cleaning is specified. When performing the cleaning,\n"
     "proteins shorter that this value will be discarded.\n"
     "[default minimum length for each protein is 50]\n")
    # minimum number of unversal proteins mapped
    arg('--min_uprots', type=int, default=0, help=
        "This parameter is used to filter both usearch and tblastn mapping phases.\n"
        "Genomes that mapp less than this number of universal proteins, will be discarded\n"
        "[default minimum number of universal proteins mapped is 0, i.e., no genomes will be omitted]\n")
    # user defined data folder
    arg('--c_dat', type=str, default=None, help=
        'Custom path to the folder where to store the intermediate files\n')
    # user defined input folder
    arg('--c_in', type=str, default=None, help=
        'Custom path to the folder containing the folder with the input data\n')
    # user defined output folder
    arg('--c_out', type=str, default=None, help=
        'Custom path to the output folder where to save the results\n')
    # user defined universal protenis folder
    arg('--c_up', type=str, default=None, help=
        'Path to the file (fasta) or folder containing the *.faa files to use as universal proteins. Default None.\n')

    arg('-v', '--version', action='store_true', help=
        "Prints the current PhyloPhlAn version and exit\n")

    return vars(p.parse_args())


def init():
    cf_fna, cf_udb = None, None

    if download:
        for f in [ppa_fna, ppa_aln, ppa_xml, ppa_up2prots, ppa_ors2prots]:
            u = urllib2.urlopen(download_compressed(f))

            if not os.path.exists(f):
                info("Downloading and extracting "+f+"... ")

                with closing(tarfile.open(fileobj=StringIO(u.read()))) as inp:
                    inp.extractall(path=os.path.dirname(f))

                info("Done!\n")
    else:
        for f in [ppa_fna, ppa_fna_40, ppa_tax, ppa_aln, ppa_xml, ppa_up2prots, ppa_ors2prots]:
            if not os.path.exists(f):
                info("Extracting \""+f+"\"... ")

                with closing(tarfile.open(compressed(f))) as inp:
                    inp.extractall(path=os.path.dirname(f))

                info("Done!\n")

        for t, f in [ppa_alns]:
            if not os.path.exists(t):
                info("Extracting \""+f+"\"... ")

                with closing(tarfile.open(f)) as inp:
                    inp.extractall(path=os.path.dirname(f))

                info("Done!\n")

    if not os.path.exists(ppa_udb):
        info("Generating \""+ppa_udb+"\" (usearch indexed DB)...")
        sb.call(["usearch", "-quiet", "-makeudb_ublast", ppa_fna, "-output", ppa_udb]) # usearch8.0.1623
        info("Done!\n")

    if cf_up:
        if os.path.isfile(cf_up): # it is the fasta file containing all the universal proteins
            cf_fna = cf_up
            cf_udb = cf_up[:cf_up.rfind('.')+1]+'udb'
        elif os.path.isdir(cf_up): # it is a folder with a fasta file for each universal protein
            if (cf_up == './') or (cf_up == '../'):
                cf_up_filename = 'cf_up'
            else:
                cf_up_filename = [i for i in cf_up.split('/') if i][-1]

            cf_fna = cf_up+cf_up_filename+'.fna'
            cf_udb = cf_up+cf_up_filename+'.udb'

            if not os.path.isfile(cf_fna):
                with open(cf_fna, 'w') as f:
                    for i in iglob(cf_up+'*'):
                        with open(i) as g:
                            f.write(g.read())
            else:
                info('File "'+cf_fna+'" already present in "'+cf_up+'" folder.\n')
        else: # what's that??
            error("Not sure about the format of your custom set of universal proteins. I'll use the default ones!\n")
            return None, None

        if cf_fna and cf_udb:
            if not os.path.isfile(cf_udb):
                info("Generating custom \""+cf_udb+"\" (usearch indexed DB)...")
                sb.call(["usearch", "-quiet", "-makeudb_ublast", cf_fna, "-output", cf_udb])
                info("Done!\n")
            else:
                info('usearch database "'+cf_udb+'" already present in "'+cf_up+'" folder.\n')

    return cf_fna, cf_udb


def clean_all(loc_dat):
    not_rm = []

    for f in chain(iglob(loc_dat+'*.txt'), iglob(loc_dat+'*.faa'), iglob(loc_dat+'*.xml'),
                   iglob(loc_dat+'*.udb'), iglob(loc_dat+'*.wdb'), iglob(os.path.dirname(ppa_alns[1])+'/*.txt'),
                   iglob(os.path.dirname(ppa_alns[1])+'/*.score'), iglob(os.path.dirname(ppa_alns[1])+'/*.aln')):
        if os.path.isfile(f):
            os.remove(f)
        else:
            not_rm.append(f)

    if not_rm:
        msg = "The following file"
        msg += "s were" if len(not_rm) > 1 else " was"
        msg += " not removed: "
        info(msg + ', '.join(not_rm) + "\n")


def clean_project():
    if os.path.exists(dat_fol):
        shutil.rmtree(dat_fol)

    if os.path.exists(out_fol):
        shutil.rmtree(out_fol)


def get_inputs(proj, params):
    if not os.path.isdir(inp_fol):
        error("No "+proj+" folder found in '"+inp_fol+"'", init_new_line=True, exit=True)

    files = list(os.listdir(inp_fol))
    faa_in = []
    fna_in = []

    for f in files:
        if 'faa' in f.split('.'):
            cut = -2 if f.endswith('.bz2') else -1
            faa_in.append('.'.join(f.split('.')[:cut]))

        if 'fna' in f.split('.'):
            cut = -2 if f.endswith('.bz2') else -1
            fna_in.append('.'.join(f.split('.')[:cut]))

    if not (len(faa_in) + len(fna_in)):
        error("No '.faa' or '.fna' input files found in '"+inp_fol+"'", init_new_line=True, exit=True)

    txt_in = [l for l in files if os.path.splitext(l)[1] == ".txt"]
    if len(txt_in) > 1:
        error("More than one '.txt' input files found in input/\n"
              "[No more than one txt file (the taxonomy) allowed", init_new_line=True, exit=True)

    tax_in = [l for l in files if os.path.splitext(l)[1] == ".tax"]
    if len(tax_in) > 1:
        error("More than one '.tax' input files found in input/\n"
              "[No more than one tax file (the taxonomy for taxonomic analysis) allowed", init_new_line=True, exit=True)

    mdt_in = [l for l in files if os.path.splitext(l)[1] == ".metadata"]
    if len(tax_in) > 1:
        error("More than one '.metadata' input files found in input/\n"
              "[No more than one metadata file allowed", init_new_line=True, exit=True)

    if not os.path.isdir(dat_fol):
        os.mkdir(dat_fol)

    # clean the faa files by filter on the number of proteins and their minimum length
    if faa_in and params['faa_cleaning']:
        t0 = time.time()
        faa_in = clean_faa(proj, faa_in, params['min_num_proteins'], params['min_len_protein']) # , verbose=True)
        info("\nCleaning faa inputs took "+str(int(time.time()-t0))+" s\n")

    fna_in = uniq_filenames(faa_in, fna_in) # make uniq filename for fna inputs
    if not os.path.exists(dat_fol+ors2prots):
        with open(dat_fol+ors2prots, 'w') as outf:
            for f in faa_in:
                fld = inp_fol if not params['faa_cleaning'] else dat_cln_fol

                if os.path.isfile(fld+f+'.faa'):
                    prots = sorted([l.split()[0][1:] for l in open(fld+f+'.faa') if l.startswith('>')])
                else:
                    with BZ2File(fld+f+'.faa.bz2', 'rU') as inp:
                        prots = sorted([l.split()[0][1:] for l in inp if l.startswith('>')])
                outf.write('\t'.join([f] + prots) + '\n')

            for f, fo in fna_in:
                if os.path.isfile(inp_fol+f+'.fna'):
                    prots = sorted([l.split()[0][1:] for l in open(inp_fol+f+'.fna') if l.startswith('>')])
                else:
                    with BZ2File(inp_fol+f+'.fna.bz2', 'rU') as inp:
                        prots = sorted([l.split()[0][1:] for l in inp if l.startswith('>')])
                outf.write('\t'.join([fo] + prots) + '\n')

    return faa_in, fna_in, txt_in[0] if txt_in else None, tax_in[0] if tax_in else None, mdt_in[0] if mdt_in else None


def check_inp(inps):
    init_dig = [l for l in inps if l[0].isdigit()]

    if init_dig:
        error("The following genome IDs start with a digit:\n    "+"\n    ".join(init_dig) +
              "\nPlease rename accordingly the corresponding input file names", init_new_line=True, exit=True)

    ppa_ids = [l[1:] for l in open(ppa_aln) if l.startswith('>')]
    init_dup = [l for l in inps if l in ppa_ids]

    if init_dup:
        error("The following genome IDs are already in PhyloPhlAn:\n    "+"\n    ".join(init_dig),
              init_new_line=True, exit=True)


def clean_faa(proj, faa_in, min_num_proteins, min_len_protein, verbose=False):
    if not os.path.isdir(dat_cln_fol): # create the directory if it does not exists
        os.mkdir(dat_cln_fol)
    elif os.path.isfile(dat_fol+f_protein_ids_map): # load the protein_ids_map mapping and return it
        with open(dat_fol+f_protein_ids_map) as f:
            return pickle.load(f).keys()

    protein_ids_map = dict() # will contains the old protein id and the new protein id assigned

    for f in faa_in:
        records = []
        out = []
        seqid = 0
        msg = ""
        old_id = dict()

        if os.path.isfile(inp_fol+f+'.faa'):
            with open(inp_fol+f+'.faa', 'rU') as h_in:
                records = list(SeqIO.parse(h_in, "fasta"))
        elif os.path.isfile(inp_fol+f+'.faa.bz2'):
            with BZ2File(inp_fol+f+'.faa.bz2', 'rU') as h_in:
                records = list(SeqIO.parse(h_in, "fasta"))
        else:
            info(' '.join(["[clean_faa]", "File:", f, "not found (.faa or .faa.bz2)"]) + "\n")
            continue

        if len(records) >= min_num_proteins: # check the number of proteins of the genome
            for protein in records:
                if len(protein.seq) >= min_len_protein: # check the length of the protein
                    pid = f+'_'+str(seqid)
                    out.append(SeqRecord(protein.seq, id=pid, description=''))
                    old_id[pid] = protein.id
                    seqid += 1
        else:
            msg = "not enough proteins!"

        # if some proteins has been removed, checked that the good ones kept are at least min_num_proteins
        if len(out) >= min_num_proteins:
            if len(records) - len(out):
                msg = str(len(records)-len(out)) + " too short proteins removed!"

            # save the old and new mapping
            if f in protein_ids_map:
                app = protein_ids_map[f].append()
                for t in ((old_id[s.id], s.id) for s in out):
                    app(t)
            else:
                protein_ids_map[f] = [(old_id[s.id], s.id) for s in out]

            # write the new .faa
            with open(dat_cln_fol+f+'.faa', 'w') as h_out:
                SeqIO.write(out, h_out, "fasta")
        else:
            msg = "not enough good proteins!"

        if msg and verbose:
            info(' '.join(["[clean_faa]", f, msg]) + "\n")

    # write the mapping of the old and new proteins id of all genomes
    with open(dat_fol+old2newprots, 'w') as f:
        for k, v in protein_ids_map.iteritems():
            f.write('\t'.join([k] + ['('+a+','+b+')' for a,b in v]))

    # serialize the protein_ids_map mapping
    with open(dat_fol+f_protein_ids_map, 'w') as f:
        pickle.dump(protein_ids_map, f, pickle.HIGHEST_PROTOCOL)

    return protein_ids_map.keys()


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def exe_usearch(x):

    def screen_usearch_wdb(inpf):
        tab = ((l[0].split(' ')[0], l[1].split('_')[1], float(l[-1])) for l in (ll.strip().split('\t') for ll in open(inpf)))
        f2p, f2b, p2b = {}, {}, {}

        for f,p,b in tab:
            if f in f2p and b < f2b[f]:
                continue
            if p in p2b and b < p2b[p]:
                continue
            if p in p2b:
                for ff in f2p.keys():
                    if p in f2p[ff]:
                        del f2p[ff]
                        del f2b[ff]

            f2p[f], f2b[f], p2b[p] = p, b, b

        with open(inpf, "w") as outf:
            # there should be at least min_uprots of universal proteins mapped
            if len(f2p) >= min_uprots:
                for k, v in f2p.iteritems():
                    outf.write("\t".join([str(k), str(v)]) + "\n")
            else:
                outf.write("\t".join([inpf[inpf.rfind('/')+1:], NOT_ENOUGH_MAPPINGS]) + "\n")

    #
    if not terminating.is_set():
        try:
            info("Starting "+x[7][x[7].rfind('/')+1:]+"\n")
            tmp_faa = None
            t = None

            if not os.path.isfile(x[3]):
                tmp_faa = NamedTemporaryFile(suffix='.faa', prefix=x[3][x[3].rfind('/')+1:x[3].find('.')],
                    dir=x[3][:x[3].rfind('/')+1])
                with open(tmp_faa.name, 'w') as h_out:
                    with BZ2File(x[3]+'.bz2', 'rU') as h_in:
                        records = list(SeqIO.parse(h_in, "fasta"))

                    SeqIO.write(records, h_out, "fasta")

                cmd = x[:3] + [tmp_faa.name] + x[4:]
            else:
                cmd = x

            t = sb.Popen(cmd)
            t.wait()
            shutil.copy2(x[7], x[7]+'.bkp') # copy the results of usearch
            screen_usearch_wdb(x[7])

            if tmp_faa:
                tmp_faa.close()

            info(x[7][x[7].rfind('/')+1:]+" generated!\n")
        except:
            if t:
                t.kill()

            if tmp_faa:
                tmp_faa.close()

            terminating.set()
            error(' '.join(x))
            return
    else:
        terminating.set()

    return


def faa2ppafaa(inps, nproc, proj, faa_cleaning):
    loc_inp = dat_cln_fol if faa_cleaning else inp_fol

    terminating = mp.Event()
    pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)
    mmap = [(loc_inp+i+'.faa', dat_fol+i+'.b6o') for i in inps if not os.path.exists(dat_fol+i+'.b6o')]

    if not os.path.isdir(dat_fol): # create the data folder if does not exists
        os.mkdir(dat_fol)

    if not mmap:
        info("All usearch runs already performed!\n")
    else:
        info("Looking for PhyloPhlAn proteins in input faa files\n")
        us_cmd = [["usearch", "-quiet", "-ublast", i, "-db", cf_udb if cf_udb else ppa_udb, "-blast6out", o, "-threads", "1",
            "-evalue", "1e-10", "-maxaccepts", "8", "-maxrejects", "32"] for i, o in mmap] # usearch8.0.1517

        try:
            pool.map(exe_usearch, us_cmd)
            pool.close()
        except:
            pool.terminate()
            pool.join()
            error("Quitting faa2ppafaa()")
            return

    if mmap:
        pool.join()
        info("All usearch runs performed!\n")

    if not os.path.exists(dat_fol+up2prots):
        up2p = collections.defaultdict(list)
        too_few_maps = set()

        for i in inps:
            for p, up in (l.strip().split('\t') for l in open(dat_fol+i+'.b6o').readlines()):
                if not up.startswith(NOT_ENOUGH_MAPPINGS):
                    up2p[up].append(p)
                else:
                    too_few_maps |= set([i])

        # write the skipped genomes
        if too_few_maps:
            with open(dat_fol+few_maps, 'w') as fout:
                for e in too_few_maps:
                    fout.write(e + '\n')

        # write the mapping between universal proteins and proteins mapped
        if up2p:
            with open(dat_fol+up2prots, 'w') as outf:
                for k, v in up2p.iteritems():
                    outf.write("\t".join([k]+v) + "\n")


def blastx_exe(x):
    if not terminating.is_set():
        try:
            info("Starting "+x[6][x[6].rfind('/')+1:]+"\n")
            tmp_fna = None
            t = None

            if not os.path.isfile(x[2]):
                tmp_fna = NamedTemporaryFile(suffix='.fna', prefix=x[2][x[2].rfind('/')+1:x[2].find('.')],
                    dir=x[2][:x[2].rfind('/')+1])
                with open(tmp_fna.name, 'w') as h_out:
                    with BZ2File(x[2]+'.bz2', 'rU') as h_in:
                        records = list(SeqIO.parse(h_in, "fasta"))

                    SeqIO.write(records, h_out, "fasta")

                cmd = x[:2] + [tmp_fna.name] + x[3:]
            else:
                cmd = x

            with open(os.devnull, 'w') as devnull:
                t = sb.Popen(cmd, stderr=devnull) # tblastn quiet homemade!
                t.wait()

            if tmp_fna: tmp_fna.close()
            info(x[6][x[6].rfind('/')+1:]+" generated!\n")
        except:
            if t: t.kill()
            if tmp_fna: tmp_fna.close()
            terminating.set()
            error(' '.join(x))
            return
    else:
        terminating.set()

    return


def blast(inps, nproc, proj, blast_full=False):
    loc_dat = dat_fol+"tblastn/"
    terminating = mp.Event()
    pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)
    mmap = [(inp_fol+i+'.fna', loc_dat+io+'.b6o') for i, io in inps.iteritems() if (not os.path.exists(loc_dat+io+'.b6o'))]

    if not os.path.isdir(loc_dat): # create the tmp directory if does not exists
        os.mkdir(loc_dat)

    if not mmap:
        info('All tblastn runs already performed!\n')
    else:
        info('Looking for PhyloPhlAn proteins in input fna files\n')

        if cf_up: # custom set of universal proteins
            dataset = cf_fna
            evalue = '1e-40'
        elif blast_full: # full dataset
            dataset = ppa_fna
            evalue = '1e-60'
        else: # small dataset
            dataset = ppa_fna_40
            evalue = '1e-40'

        us_cmd = [['tblastn', '-subject', i, '-query', dataset, '-out', o, '-outfmt', '6', '-evalue', evalue]
                  for i, o in mmap] # tblastn2.2.28+

        try:
            pool.map(blastx_exe, us_cmd)
            pool.close()
        except:
            pool.terminate()
            pool.join()
            error("Quitting blast()")
            return

        pool.join()
        info('All tblastn runs performed!\n')


def gens2prots(inps, proj, faa_cleaning):
    loc_inp = dat_cln_fol if faa_cleaning else inp_fol

    if not os.path.isdir(dat_fol): # create the tmp directory if does not exists
        os.mkdir(dat_fol)

    if os.path.exists(dat_fol+ups2faa_pkl):
        return

    ups2prots = dict([(ll[0], set(ll[1:])) for ll in (l.strip().split('\t') for l in open(dat_fol+up2prots))])
    prots2ups = {}

    for k, v in ups2prots.items():
        for vv in v:
            if vv in prots2ups:
                error(str(vv) + " already in dict!")
            prots2ups[vv] = k

    ups2faa = collections.defaultdict(set)
    e = None

    for i in inps:
        inp = None

        if os.path.isfile(loc_inp+i+'.faa'):
            inp = open(loc_inp+i+'.faa', 'rU')
        else:
            inp = BZ2File(loc_inp+i+'.faa.bz2', 'rU')

        for l in inp:
            if l.startswith(">"):
                if e in prots2ups:
                    ups2faa[prots2ups[e]].add(s+"\n")

                e, s = l.strip().split()[0][1:], ""

            s += l

        if e in prots2ups:
            ups2faa[prots2ups[e]].add(s+"\n")

        if inp:
            inp.close()

    for k, v in ups2faa.items():
        with open(dat_fol+k+'.faa', 'w') as outf:
            for vv in v:
                outf.write("".join(vv))

    with open(dat_fol+ups2faa_pkl, 'w') as outf:
        pickle.dump(ups2faa, outf, pickle.HIGHEST_PROTOCOL)


def screen(stat, cols, sf=None, unknown_fraction=0.5, n=10):
    lena, nsc = len(cols[0]), []

    if sf and os.path.exists(sf):
        with open(sf) as sfinp:
            scores = [(ll[:12], ll[12:]) for ll in [l for l in sfinp.readlines()]]
    else:
        scores = [("1.0", c) for c in cols]

    for sc, st in zip(scores, stat):
        ssc, lsc = sc

        try:
            float(ssc)
        except:
            continue

        if (len(st) == 1) or \
           (len(st) == 2 and ("-" in st or "X" in st)) or \
           (len(st) == 3 and "X" in st and "-" in st):
            continue

        nn = 0
        for k in ["-","X"]:
            if k in st:
                nn += st[k]

        if nn > lena*unknown_fraction:
            continue

        nsc.append((float(ssc), lsc.strip()))

    nsc = sorted(nsc, key=lambda x: x[0], reverse=True)
    ret = [v for _, v in nsc][-n:]
    # ret = [v for _, v in nsc]
    # return ret if ret else [list(['-']*lena+['\n'])] # I'm not sure about the '\n', I believe that SeqIO.write() will complain!
    return ret if ret else [list(['-']*lena)]


def aln_subsample(inp_f, out_f, scores, unknown_fraction, namn):
    with open(inp_f) as inp:
        fnas = list(SeqIO.parse(inp, "fasta"))

    ids = [f.id for f in fnas]
    cols = zip(*[f.seq for f in fnas])
    col_alf = [set(x) for x in cols]
    col_sta = [dict([(a, seq.count(a)) for a in s]) for s, seq in zip(col_alf, cols)]
    col_sta_ok = screen(col_sta, cols, sf=scores, unknown_fraction=unknown_fraction, n=namn)
    raws = zip(*col_sta_ok)
    recs = [SeqRecord(seq=Seq("".join(r)), id=nid, description=nid) for r, nid in zip(raws, ids)]

    with open(out_f,"w") as out:
        SeqIO.write(recs, out, "fasta")


def exe_muscle(x):
    if not terminating.is_set():
        try:
            if len(x) < 14:
                info("Running muscle on "+x[3] + "\n")
            else:
                info("Running muscle on "+x[4] + " and "+x[6]+"\n")

            t = None
            with open(os.devnull, 'w') as devnull:
                t = sb.Popen(x[:-2], stderr=devnull) # quiet mode
                t.wait()

            if x[9]:
                pn = 20
            else:
                pn = max(int(max(int((400.0-int(x[-1][1:]))*30.0/400.0), 1)**2 /30.0), 3)

            if len(x) < 14:
                aln_subsample(x[5], x[-2], x[7], 0.1, pn)
                info(x[-2] + " generated (from "+x[3]+")!\n")
            else:
                aln_subsample(x[8], x[-2], x[10], 0.1, pn)
                info(x[-2] + " generated (from "+x[4]+" and "+x[6]+")!\n")
        except:
            if t: t.kill()
            terminating.set()
            error(' '.join(x))
            return
    else:
        terminating.set()

    return


def faa2aln(nproc, proj, integrate=False):
    ppa_fol = cf_up if cf_up else ppa_alns_fol

    with open(dat_fol+ups2faa_pkl) as inp:
        prots = pickle.load(inp)

    if cf_up:
        up_list = (i[i.rfind('/')+1:i.rfind('.')] for i in iglob(cf_up+'*')) # assumes the folder contains all and only the needed files
    else:
        up_list = ('p{num:04d}'.format(num=aa) for aa in range(400))

    mmap = [(dat_fol+i+".faa", dat_fol+i+".aln", dat_fol+i+".sc",
             ppa_fol+i+".aln", ppa_fol+i+".aln.score",
             dat_fol+i+".int.aln", dat_fol+i+".int.sc",
             dat_fol+i+".sub.aln", dat_fol+i+".int.sub.aln",
             cf_up, i in prots, i) for i in up_list]

    if mmap:
        if integrate:
            for _, _, _, _, _, _, _, _, _, _, _, i in mmap:
                with open(dat_fol+i+".faa", 'a') as f:
                    with open(ppa_fol+i+".faa", 'r') as g:
                        f.write(g.read())

        us_cmd = [["muscle", "-quiet",
                   "-in", i,
                   "-out", o,
                   "-scorefile", s,
                   "-maxiters", "2",
                   so, pn] for i,o,s,o2,s2,oi,si,so,soi,up,pres,pn in mmap if not os.path.exists(o) and pres]

        if us_cmd:
            info("Looking for PhyloPhlAn proteins to align\n")
            info(str(len(us_cmd))+" alignments to be performed\n")
            terminating = mp.Event()
            pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

            try:
                pool.map(exe_muscle, us_cmd)
                pool.close()
            except:
                pool.terminate()
                pool.join()
                error("Quitting faa2aln()")
                return

            pool.join()
            info("All alignments performed!\n")

    # if integrate:
    #     us_cmd = [ ["muscle","-quiet",
    #                 "-profile",
    #                 "-in1",o,
    #                 "-in2",o2,
    #                 "-out",oi,
    #                 "-scorefile",si,
    #                 "-maxiters","2",
    #                 soi,pn] for i,o,s,o2,s2,oi,si,so,soi,pres,pn in mmap
    #                     if not os.path.exists(soi) and pres]
    #     for i,o,s,o2,s2,oi,si,so,soi,pres,pn in mmap: #paralellizable
    #         if not os.path.exists(soi) and not pres:
    #             pnn = max( int( max( int((400.0-int(pn[1:]))*30.0/400.0),1)**2 / 30.0), 3)
    #             aln_subsample( o2, soi, s2, 0.1, pnn )
    #     if us_cmd:
    #         info("Merging alignments from user genomes with all the PhyloPhlAn alignments\n")
    #         info(str(len(us_cmd))+" alignments to be merged\n")
    #         terminating = mp.Event()
    #         pool = mp.Pool(initializer=initt, initargs=(terminating, ), processes=nproc)

    #         try:
    #             pool.map(exe_muscle, us_cmd)
    #             pool.close()
    #         except:
    #             pool.terminate()
    #             pool.join()
    #             error("Quitting faa2aln() - integrate")
    #             return

    #         pool.join()
    #         info("All alignments already merged with PhyloPhlAn alignments!\n")

    if os.path.exists(dat_fol+up2prots):
        info("All alignments already computed!\n")
        return

    # up2p = collections.defaultdict( list )
    # for i in inps:
    #     for p,up in (l.strip().split('\t')
    #             for l in open(dat_fol+i+".b6o").readlines()):
    #         up2p[up].append( p )
    # if integrate:
    #     with open(ppa_up2prots) as inpf:
    #         for l in (ll.strip().split('\t') for ll in inpf):
    #             up2p[l[0]] += l[1:]

    # with open(dat_fol+up2prots,"w") as outf:
    #     for k,v in up2p.items():
    #         outf.write( "\t".join([k]+v) + "\n" )


def aln_merge(proj, integrate):
    loc_dat = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot

    if os.path.exists(loc_dat):
        info("All alignments already merged!\n")
        return

    up2p = dict([(l[0], set(l[1:])) for l in (ll.strip().split('\t') for ll in open(dat_fol+up2prots))])

    if integrate:
        for l in (ll.strip().split('\t') for ll in open(ppa_up2prots)):
            if l[0] in up2p:
                up2p[l[0]] |= set(l[1:])
            else:
                up2p[l[0]] = set(l[1:])

    # check if there are proteins id duplicated and if yes warn the user and exit
    if collections_counter:
        print "TEST", [k for k, v in Counter(sum([list(x) for x in up2p.values()], [])).iteritems() if v > 1]

    all_prots = set.union(*up2p.values()) # all proteins id
    genomes_to_skip = [r.strip() for r in open(dat_fol+few_maps, 'r')] if os.path.isfile(dat_fol+few_maps) else []
    t2p = dict([(l[0], set(l[1:])) for l in (ll.strip().split('\t') for ll in open(dat_fol+ors2prots)) if l[0] not in genomes_to_skip])

    if integrate:
        for l in (ll.strip().split('\t') for ll in open(ppa_ors2prots)):
            t2p[l[0]] = set(l[1:])

    # map proteins id to genomes id
    p2t = dict()
    if os.path.exists(dat_fol+p2t_map):
        p2t = dict([(l[0], l[1]) for l in (row.strip().split('\t') for row in open(dat_fol+p2t_map, 'r'))])

    for t, p in t2p.iteritems():
        for pp in p:
            if pp in all_prots: # if the protein is present in the current alignments
                p2t[pp] = t  # map proteins to their genome

    all_g = set(t2p.keys()) # all genomes id
    aln = dict([(t, "") for t in t2p.keys()]) # dictionary that will contains all alignments
    up = dict([(t, 0) for t in t2p.keys()]) # dictionary that will contains the number of universal proteins mapped

    for f in sorted((dat_fol+ff for ff in os.listdir(dat_fol) if ff.endswith(".sub.aln"))):
        g2aln = dict()
        lenal = 0 # keep the lenght of the current alignments

        for seq in SeqIO.parse(f, "fasta"):
            if not lenal:
                lenal = len(seq.seq)

            g2aln[p2t[seq.id]] = seq

        for g in all_g: # for all genomes
            if g in g2aln: # if there is alignment
                up[g] += 1
                aln[g] += g2aln[g]
            else: # otherwise add gaps
                aln[g] += SeqRecord(Seq('-'*lenal), id=str(g))

    out_faas = []
    for k, v in aln.iteritems():
        v.id = k
        v.description = ""
        out_faas.append(v)

    SeqIO.write(out_faas, loc_dat, "fasta")
    info("All alignments merged into "+loc_dat+"!\n")

    if not os.path.isdir(out_fol):
        os.mkdir(out_fol) # create output directory if does not exists

    # write statistics file
    info("Writing "+out_fol+"aln_stats.tsv ")
    with open(out_fol+'aln_stats.tsv', 'w') as f:
        f.write('\t'.join(['id', 'tot_length', 'aln_length', 'tot_gaps', 'up_mapped'])+'\n')

        for k, v in aln.iteritems():
            tot_len = len(v)
            tot_gaps = str(v.seq).count('-')
            f.write('\t'.join([str(k), str(tot_len), str(tot_len-tot_gaps), str(tot_gaps), str(up[k])])+'\n')

    info('Done!\n')


def fasttree(proj, integrate):
    # dat_fol = cf_data if cf_data else "data/"
    # dat_fol += proj+"/"
    # out_fol = cf_output if cf_output else "output/"
    # out_fol += proj+"/"
    loc_out = out_fol+proj+"."+ (o_inttree if integrate else o_tree)

    if not os.path.isdir(out_fol):
        os.mkdir(out_fol)

    aln_in = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot

    if os.path.exists(loc_out):
        info("Final tree already built ("+loc_out+")!\n")
        return

    info("Start building the tree with FastTree\n")
    cmd = ["FastTree", "-quiet",
           # "-fastest", "-noml", "-gamma", "-pseudo",
           "-bionj", "-slownni", "-mlacc", "2", "-spr", "4"]
    t = None
    with open(os.devnull, 'w') as devnull:
        t = sb.Popen(cmd, stdin=open(aln_in), stdout=open(loc_out, 'w'), stderr=devnull) # quiet mode
        t.wait()

    # convert aln_in from fasta to phylip format
    # cmd = [ "raxmlHPC-PTHREADS-SSE3",
    #         "-T","#threads",
    #         "-s","alignment in phylip format",
    #         "-n",loc_out,
    #         "-m","PROTCATWAG",
    #         "-p","1982"
    #         ]
    # sb.call( cmd )
    info("Tree built! The output newick file is in "+loc_out+"\n")


def circlader(proj, integrate, tax=None):
    # inp_fol = cf_input if cf_input else "input/"
    # inp_fol += proj+"/"
    # dat_fol = cf_data if cf_data else "data/"
    # dat_fol += proj+"/"
    # out_fol = cf_output if cf_output else "output/"
    # out_fol += proj+"/"
    # out_img = "output/"+proj+"/imgs/"

    tree = out_fol+proj+"."+ (o_inttree if integrate else o_tree)
    xtree_rr = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"xml"
    xtree_an = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
    pngtree = out_fol+proj+".tree"+ (".int" if integrate else "")+".png"
    annotation_f = dat_fol+proj+(".int" if integrate else "")+".annot.txt"
    archaea_f = dat_fol+proj+(".int" if integrate else "")+".archaea.txt"
    tax_d = [l.strip().split('\t') for l in open(inp_fol+tax)] if tax else []
    # o2t = dict([(o,t.split('.')) for t,o in tax_d])
    t2o = dict(tax_d)

    if os.path.exists(pngtree):
        info("Output tree images already generated!\n")
        return

    # info("Generating tree output images ...")
    c_circ = "circlader/circlader.py"
    c_circ_an = "circlader/circlader_annotate.py"
    c_reroot = "./taxcuration/tree_reroot.py"
    names = [l.strip().split('\t')[0] for l in open(dat_fol+ors2prots)]

    if integrate:
        return

        int_names = [l.strip().split('\t')[0] for l in open(ppa_ors2prots)]
        mt4tax_d = [l.strip().split('\t')[::-1] for l in open(ppa_tax)]
        # mt4o2t = dict([(o,t.split('.')) for t,o in mt4tax_d])
        mt4t2o = dict(mt4tax_d)
        mt4archaea = [o for t,o in mt4t2o.items() if t.split(".")[0] == "d__Archaea"]

        with open(archaea_f, "w") as out:
            out.write("\n".join(mt4archaea) +"\n")

        sb.call([c_reroot, "-s", "lca", "-f", archaea_f, tree, xtree_rr])
        annotation_f = dat_fol+proj+".annot.txt"

        with open(annotation_f, "w") as outf:
            for v1,v2 in [('clade_fill_color','r')]:
                outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n")

            for v1,v2 in [('clade_size','7.0'), ]:
                outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names+int_names]) + "\n")

            for v1,v2 in [('clade_size','2.0')]:
                outf.write("\t".join(["*",v1,v2]) + "\n")

            for v1,v2 in [('sub_branches_angle_reduction','0.0'), ('clade_edge_line_width','0.2'), ('branch_line_width','0.45'), ]:
                outf.write("\t".join([v1,v2]) + "\n")
    else:
        n = 1 if len(names) < 4 else int(round(min(len(names)*0.05,100)))
        archaea = [o for t,o in t2o.items() if t.split(".")[0] == "d__Archaea"]

        if t2o and archaea:
            with open(archaea_f, "w") as out:
                out.write("\n".join(archaea) +"\n")

            sb.call([c_reroot, "-s", "lca", "-f", archaea_f, tree, xtree_rr])
        elif n == 1 and len(names) < 8:
            sb.call([c_reroot, "-s","longest_edge", tree, xtree_rr])
        elif n == 1:
            sb.call([c_reroot, "-s", "longest_internal_edge", tree, xtree_rr])
        else:
            sb.call([c_reroot, "-s", "longest_internal_edge_n", "-n", str(n), tree, xtree_rr])

        return

        with open(annotation_f, "w") as outf:
            for v1,v2 in [('clade_fill_color','r'), ('clade_size','7.0'), ]:
                outf.write("\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n")

            for v1,v2 in [('clade_size','2.0')]:
                outf.write("\t".join(["*",v1,v2]) + "\n")

            for v1,v2 in [('sub_branches_angle_reduction','0.0'), ('clade_edge_line_width','0.2'), ('branch_line_width','0.45'), ]:
                outf.write("\t".join([v1,v2]) + "\n")

    sb.call([c_circ_an, "--annot", annotation_f, xtree_rr, xtree_an])
    sb.call([c_circ, "--dpi", "300", xtree_an, pngtree])
    info("Done!\n")


def tax_curation(proj, tax, integrate=False, mtdt=None, inps=None):
    # inp_fol = cf_input if cf_input else "input/"
    # inp_fol += proj+"/"
    # out_fol = cf_output if cf_output else "output/"
    # out_fol += proj+"/"
    taxin = ppa_tax if integrate else inp_fol+tax
    taxout = None if integrate else out_fol+"new."+tax
    taxdiffs = None if integrate else out_fol+"diff."+tax
    taxboth = None if integrate else out_fol+"comp."+tax
    intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
    operation = "imputation" if integrate else "curation"

    if not os.path.exists(intree):
        intree = out_fol+proj+".tree."+ ("int." if integrate else "")+"nwk"

    if tax:
        info("Starting taxonomic curation: detecting suspicious taxonomic labels ")
        taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, to_integrate=inp_fol+tax)
        taxCleaner.remove_tlabels()
        info("Done!\n")
        info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy ")
    else:
        taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, inps=inps)
        info("Trying to impute taxonomic labels for taxa newly integrated into the tree ")

    taxCleaner.infer_tlabels()
    info("Done!\n")

    info("Writing taxonomic "+operation+" outputs ")
    tid2img = taxCleaner.write_output(proj, images=False)
    info("Done!\n")
    info("Writing final pdf reports ")
    if integrate:
        pass
    elif mtdt:
        taxCleaner.write_report(out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images=True, typ="refined")
        taxCleaner.write_report(out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images=True, typ="corrected")

    if not integrate:
        taxCleaner.write_new_taxonomy(taxout, taxdiffs, taxboth)

    info("Done!\n")


def tax_imputation(proj, tax, integrate=False, mtdt=None, inps=None):
    # inp_fol = cf_input if cf_input else "input/"
    # inp_fol += proj+"/"
    # out_fol = cf_output if cf_output else "output/"
    # out_fol += proj+"/"
    taxin = ppa_tax if integrate else inp_fol+tax
    # taxout = None if integrate else out_fol+"new."+tax
    # taxdiffs = None if integrate else out_fol+"diff."+tax
    # taxboth = None if integrate else out_fol+"comp."+tax
    intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
    operation = "imputation" if integrate else "curation"

    if not os.path.exists(intree):
        intree = out_fol+proj+".tree."+ ("int." if integrate else "")+"nwk"

    info("Trying to impute taxonomic labels for taxa newly integrated into the tree ")
    taxCleaner = taxc.TaxCleaner(taxin, intree, integrate=integrate, inps=inps)
    taxCleaner.infer_tlabels()
    info("Done!\n")
    info("Writing taxonomic "+operation+" outputs ")
    taxCleaner.write_output(proj, images=False, imp_only=True)
    info("Done!\n")


def tax_curation_test(proj, tax, nerrorrs=10, error_type='wrong', taxl='s', tmin=3, tex=None, name='err.txt', integrate=False, descr=None):
    err_tax = out_fol+name
    intree = out_fol+proj+".tree.reroot."+ ("int." if integrate else "")+"annot.xml"
    taxCleaner = taxc.TaxCleaner(inp_fol+tax, intree)
    taxCleaner.introduce_errors(nerrorrs, err_tax, error_type, taxl, tmin, tex)
    taxCleaner = taxc.TaxCleaner(err_tax, intree)
    info("Starting taxonomic curation: detecting suspicious taxonomic labels ")
    taxCleaner.remove_tlabels()
    info("Done!\n")
    info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy ")
    taxCleaner.infer_tlabels()
    info("Done!\n")
    info("Writing taxonomic curation outputs and evaluations ")
    taxCleaner.write_output("log.txt", proj, outname=name) # , images=True )
    taxCleaner.evaluate(proj, inp_fol+tax, err_tax, name, descr=descr)
    info("Done!\n")


def notinclusters(xx, yy, clusters):
    reverse = True if xx > yy else False

    for x, y in clusters:
        if reverse and (x > y):
            if ((x <= xx) and (xx <= y)) and ((x <= yy) and (yy <= y)):
                return False
        elif not reverse and (x < y):
            if ((y <= xx) and (xx <= x)) and ((y <= yy) and (yy <= x)):
                return False

    return True


def longest_not_overlapping(points):
    spoints = sorted(points, key=lambda x: abs(x[1]-x[0]), reverse=True)
    lno = [spoints[0]]

    for a, b in spoints[1:]:
        to_add = True

        for x, y in lno:
            if not (((a < x) and (b < x) and (a < y) and (b < y)) or
                    ((a > x) and (b > x) and (a > y) and (b > y))):
                to_add = False
                break

        if to_add: lno.append((a, b))

    return lno


def find_clusters(reprs, points, f1, f2, reverse=True):
    clusters = set()
    changes = True

    while changes:
        changes = False
        tmp_clusters = set()

        for x, y in set(reprs).difference(clusters):
            xyrev = True if x > y else False

            for a, b in points:

                if reverse and xyrev:
                    a, b = b, a

                if ((x != a) and (y != b)) and notinclusters(a, b, tmp_clusters):
                    if ((a <= x) and (b >= x) and (b <= y)) or ((a >= x) and (a <= y) and (b >= y)):
                        tmp_clusters = tmp_clusters.union(set([(f1(x, a), f2(y, b))]))

            if tmp_clusters:
                break
            else:
                if (x, y) not in clusters:
                    clusters = clusters.union(set([(x, y)]))

        if tmp_clusters:
            changes = True
            to_remove = set([(x, y)])

            for a, b in reprs:
                for x, y in tmp_clusters:
                    if (a >= x and a <= y and b >= x and b <= y) and ((a, b) not in to_remove):
                        to_remove = to_remove.union(set([(a, b)]))

            reprs = list(set(reprs).difference(to_remove))
            reprs.append((f1([a for a, _ in tmp_clusters]), f2([b for _, b in tmp_clusters])))
            tmp_clusters = set()

    return clusters


def fake_proteome(proj, fna_out, faa_cleaning):
    loc_inp = dat_cln_fol if faa_cleaning else inp_fol
    b6o_files = iglob(dat_fol+'tblastn/*.b6o')
    fake_proteomes = []
    contigs = {}
    done = True
    p2t = dict()
    fin_dict = dict([(fo, fi) for fi, fo in fna_out.iteritems()])

    # check if there are already present the result files
    for _, f in fna_out.iteritems():
        if not os.path.isfile(loc_inp+f+'.faa') and done:
            done = False

    if done:
        return None

    if not os.path.isdir(loc_inp):
        os.mkdir(loc_inp) # create the directory if it does not exists

    info('Reading tblastn output files (.b6o)')

    for f in b6o_files:
        key = f[f.rfind('/')+1:f.rfind('.')]

        with open(f) as hf:
            for r in hf:
                if key not in fake_proteomes:
                    fake_proteomes.append(key)

                tmp_lst = r.strip().split('\t')
                dic = {} if key not in contigs else contigs[key]

                tmp_tmp_lst = [tmp_lst[0]]
                tmp_tmp_lst.append(float(tmp_lst[2]))
                tmp_tmp_lst += [int(e) for e in tmp_lst[3:10]]
                tmp_tmp_lst += [float(e) for e in tmp_lst[10:]]

                if tmp_lst[1] in dic:
                    dic[tmp_lst[1]].append(tmp_tmp_lst)
                else:
                    dic[tmp_lst[1]] = [tmp_tmp_lst]

                contigs[key] = dic

    info(': '+str(len(contigs))+' files loaded!\n')

    for f, cc in contigs.items():
        fin = fin_dict[f]
        clusters = {}

        for c, d in cc.items():
            clusters[c] = set()
            yerev = []
            norev = []

            for dd in d:
                if dd[7] > dd[8]:
                    yerev.append((dd[7], dd[8]))
                else:
                    norev.append((dd[7], dd[8]))

            if norev:
                lno = longest_not_overlapping(norev) # find longest and not-overlapping segments
                clusters[c] |= find_clusters(lno, norev, min, max, reverse=False) # search for clusters

            if yerev:
                lno = longest_not_overlapping(yerev) # find longest and not-overlapping segments
                clusters[c] |= find_clusters(lno, yerev, max, min, reverse=True) # search for clusters

        # open input file and extract and translate the clusters
        proteome = []
        ff = None

        if os.path.isfile(inp_fol+fin+'.fna'):
            ff = open(inp_fol+fin+'.fna')
            info('Reading input file: '+inp_fol+fin+'.fna\n')
        elif os.path.isfile(inp_fol+fin+'.fna.bz2'):
            ff = BZ2File(inp_fol+fin+'.fna.bz2')
            info('Starting '+inp_fol+fin+'.fna.bz2\n')
        else:
            info('File not found: '+inp_fol+fin+'.fna(.bz2)\n')

        if ff:
            for record in SeqIO.parse(ff, 'fasta'):
                for contig, cpoints in clusters.items():
                    if record.id in contig:
                        for s, e in cpoints:
                            reverse = False

                            if s > e:
                                s, e = e, s
                                reverse = True

                            sequence1 = Seq(str(record.seq)[s:e])

                            if s-1 < 0:
                                sequence2 = Seq('N' + str(record.seq)[s:e])
                            else:
                                sequence2 = Seq(str(record.seq)[s-1:e])

                            if s-2 < 0:
                                if s-1 < 0:
                                    sequence3 = Seq('NN' + str(record.seq)[s:e])
                                else:
                                    sequence3 = Seq('N' + str(record.seq)[s:e+1])
                            else:
                                sequence3 = Seq(str(record.seq)[s-2:e])

                            # if sequence no div by 3, add Ns
                            while (len(sequence1) % 3) != 0:
                                sequence1 += Seq('N')

                            while (len(sequence2) % 3) != 0:
                                sequence2 += Seq('N')

                            while (len(sequence3) % 3) != 0:
                                sequence3 += Seq('N')

                            if reverse:
                                sequence1 = sequence1.reverse_complement()
                                sequence2 = sequence2.reverse_complement()
                                sequence3 = sequence3.reverse_complement()

                            rev = ':c' if reverse else ':' # reverse or not

                            seqid1 = f+'_'+record.id+rev+str(s)+'-'+str(e)
                            seqid2 = f+'_'+record.id+rev+str(s-1)+'-'+str(e)
                            seqid3 = f+'_'+record.id+rev+str(s-2)+'-'+str(e)

                            aminoacids1 = Seq.translate(sequence1)
                            aminoacids2 = Seq.translate(sequence2)
                            aminoacids3 = Seq.translate(sequence3)

                            proteome.append(SeqRecord(aminoacids1, id=seqid1, description=''))
                            proteome.append(SeqRecord(aminoacids2, id=seqid2, description=''))
                            proteome.append(SeqRecord(aminoacids3, id=seqid3, description=''))

                            # save map from seqid to genome
                            p2t[seqid1] = f
                            p2t[seqid2] = f
                            p2t[seqid3] = f

            ff.close()

            # write output file
            with open(loc_inp+f+'.faa', 'w') as ff:
                SeqIO.write(proteome, ff, 'fasta')
                info(loc_inp+f+'.faa generated!\n')

    # write the partial mapping of proteins into genomes
    if p2t:
        with open(dat_fol+p2t_map, 'w') as f:
            for k, v in p2t.iteritems():
                f.write(str(k)+'\t'+str(v)+'\n')

    return fake_proteomes


def uniq_filenames(faa_in, fna_in):
    new_fna_in = []

    for fna in fna_in:
        if fna in faa_in:
            hashh = hashlib.sha1()
            hashh.update(fna)
            new_fna_in.append((fna, fna+'_'+str(hashh.hexdigest()[:4])))
        else:
            new_fna_in.append((fna, fna))

    return new_fna_in


if __name__ == '__main__':
    pars = read_params(sys.argv)

    if ('v' in pars and pars['v']) or ('version' in pars and pars['version']):
        info("PhyloPhlAn version "+__version__+" ("+__date__+")\n")
        sys.exit(0)

    global cf_data
    cf_data = pars['c_dat']
    if cf_data:
        cf_data += '/' if not cf_data.endswith('/') else ''

    if pars['cleanall']:
        if ('taxonomic_analysis' in pars and pars['taxonomic_analysis']) or \
           ('user_tree' in pars and pars['user_tree']) or ('integrate' in pars and pars['integrate']):
            error("--cleanall is in conflict with -t, -u, and -i", exit=True)
        else:
            info("Cleaning \"data/\" and \"data/ppaalns/\" folders... ")
            clean_all(cf_data if cf_data else 'data/')
            info("Done!\n")

        sys.exit(0)

    projn = pars['inp']
    if not projn:
        error("Project name not provided.", exit=True)

    global cf_input
    cf_input = pars['c_in']
    if cf_input:
        cf_input += '/' if not cf_input.endswith('/') else ''

    global cf_output
    cf_output = pars['c_out']
    if cf_output:
        cf_output += '/' if not cf_output.endswith('/') else ''

    global dat_fol
    global dat_cln_fol
    dat_fol = cf_data if cf_data else 'data/'
    dat_fol += projn+'/'
    dat_cln_fol = dat_fol+cleanfaa_fld

    global inp_fol
    global inp_cln_fol
    inp_fol = cf_input if cf_input else 'input/'
    inp_fol += projn+'/'
    inp_cln_fol = inp_fol+cleanfaa_fld

    global out_fol
    out_fol = cf_output if cf_output else 'output/'
    out_fol += projn+'/'

    if pars['clean']:
        if ('taxonomic_analysis' in pars and pars['taxonomic_analysis']) or \
           ('user_tree' in pars and pars['user_tree']) or \
           ('integrate' in pars and pars['integrate']):
            error("-c/--clean is in conflict with -t, -u, and -i", exit=True)
        else:
            info("Cleaning project \""+projn+"\"... ")
            clean_project()
            info("Done!\n")

        sys.exit(0)

    global cf_up
    cf_up = pars['c_up']
    if cf_up:
        cf_up += '/' if not cf_up.endswith('/') else ''

    dep_checks()

    global cf_fna
    global cf_udb
    cf_fna, cf_udb = init()

    info("Loading and checking input files: ")
    t0 = time.time()
    faa_in, fna_in, tax, rtax, mtdt = get_inputs(projn, pars)
    check_inp(faa_in + [f for _, f in fna_in])
    info(str(len(faa_in)+len(fna_in))+" input files loaded in "+str(int(time.time()-t0))+" s!\n")

    if not tax and pars['taxonomic_analysis'] and not pars['integrate']:
        error("No taxonomy file found for the taxonomic analysis", exit=True)

    global min_uprots
    min_uprots = pars['min_uprots']
    t0 = time.time()

    if fna_in:
        try:
            blast(dict(fna_in), pars['nproc'], projn, pars['blast_full'])
        except:
            error("Quitting PhyloPhlAn [blast error]", exit=True)

        fake = fake_proteome(projn, dict(fna_in), pars['faa_cleaning'])

        if fake:
            faa_in += fake

    if faa_in:
        try:
            faa2ppafaa(faa_in, pars['nproc'], projn, pars['faa_cleaning'])
            gens2prots(faa_in, projn, pars['faa_cleaning'])
        except:
            error("Quitting PhyloPhlAn [usearch error]", exit=True)

    t1 = time.time()
    info("Mapping finished in "+str(int(t1-t0))+" secs.\n")

    try:
        faa2aln(pars['nproc'], projn, pars['integrate'])
    except:
        error("Quitting PhyloPhlAn [muscle error]", exit=True)

    t2 = time.time()
    info("Aligning finished in "+str(int(t2-t1))+" secs ("+str(int(t2-t0))+" total time).\n")
    aln_merge(projn, pars['integrate'])
    fasttree(projn, pars['integrate'])
    t4 = time.time()
    info("Tree building finished in "+str(int(t4-t2))+" secs ("+str(int(t4-t0))+" total time).\n")
    circlader(projn, pars['integrate'], tax)

    # if pars['integrate'] and tax and not pars['taxonomic_analysis']:
    #     tax_imputation( projn, tax, mtdt = mtdt, integrate = pars['integrate'], inps = inps )

    if pars['integrate'] and pars['taxonomic_analysis']:
        tax_imputation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)
        sys.exit()

    if not pars['taxonomic_analysis']:
        sys.exit()

    if pars['tax_test']:
        nerrorrs, error_type, taxl, tmin, tex, name = pars['tax_test'].split(":")
        tax_curation_test(projn, tax,
                          nerrorrs=int(nerrorrs),
                          error_type=error_type,
                          taxl=taxl,
                          tmin=int(tmin),
                          tex=int(tex) if tex else None,
                          name=name,
                          descr=pars['tax_test'])
    else:
        tax_curation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)
