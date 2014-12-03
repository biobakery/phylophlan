#!/usr/bin/env python


__author__  = 'Nicola Segata (nsegata@hsph.harvard.edu)'
__version__ = '1.10'
__date__    = '17 Sep 2014'


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
from glob import glob
from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
sys.path.insert(0,'taxcuration/')
import taxcuration as taxc
import shutil
import time


download = ""
ppa_fna = "data/ppa.seeds.faa"
ppa_fna_40 = "data/ppa.seeds.40.faa"
ppa_aln = "data/ppafull.aln.faa"
ppa_up2prots = "data/ppafull.up2prots.txt"
ppa_ors2prots = "data/ppafull.orgs2prots.txt"
ppa_tax = "data/ppafull.tax.txt"
ppa_alns = ("data/ppaalns/list.txt","data/ppaalns/ppa.aln.tar.bz2")
ppa_alns_fol = "data/ppaalns/"
ppa_xml = "data/ppafull.xml"
ppa_wdb = "data/ppa.wdb"
up2prots = "up2prots.txt"
ors2prots = "orgs2prots.txt"
aln_tot = "aln.fna"
aln_int_tot = "aln.int.fna"
ups2faa_pkl = "ups2faa.pkl"
o_tree = "tree.nwk"
o_inttree = "tree.int.nwk"
p2t_map = "p2t_map.txt"


compressed = lambda x: x+".tar.bz2"
download_compressed = lambda x: download+os.path.basename(x)+".tar.bz2"


def info(s):
    sys.stdout.write(str(s))
    sys.stdout.flush()


def error(s):
    sys.stderr.write('[E] ' + str(s) + '\n')
    sys.stderr.write('Exiting... \n')
    sys.exit(1)


def dep_checks():
    for prog in ["FastTree", "usearch", "muscle", "tblastn"]:
        try:
            with open(os.devnull, 'w') as devnull:
                sb.call([prog], stdout=devnull, stderr=devnull)
        except OSError:
            error(prog+" not found or not in system path")


def read_params(args):
    p = ap.ArgumentParser( description=
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
             formatter_class=ap.RawTextHelpFormatter )
    arg = p.add_argument

    arg( 'inp', metavar='PROJECT NAME', type=str, default=None, nargs='?', help=
        "The basename of the project corresponding to the name of the input data folder inside \n"
        "input/. The input data consist of a collection of multifasta files (extension .faa)\n"
        "containing the proteins in each genome. \n"
        "If the project already exists, the already executed steps are not re-ran.\n"
        "The results will be stored in a folder with the project basename in output/\n"
        "Multiple project can be generated and they safetely coexists." )

    arg( '-i','--integrate', action='store_true', help=
         "Integrate user genomes into the PhyloPhlAn tree \n")
    arg( '-u','--user_tree', action='store_true', help=
         "Build a phylogenetic tree using user genomes only \n" )
    arg( '-t','--taxonomic_analysis', action='store_true', help=
         "Check taxonomic inconsistencies and refine/correct taxonomic labels\n")
    arg( '--tax_test', type=str, default = None, help=
         "nerrors:type:taxl:tmin:tex:name (alpha version, experimental!)\n" )

    arg( '-c','--clean', action='store_true', help=
         "Clean the final and partial data produced for the specified project.\n"
         " (use --cleanall for removing general installation and database files)\n")
    arg( '--cleanall', action='store_true', help=
            "Remove all instalation and database file leaving untouched the initial compressed data \n"
            "that is automatically extracted and formatted at the first pipeline run.\n"
            "Projects are not remove (specify a project and use -c for removing projects).\n")

    arg( '--nproc', metavar="N", type=int, default=1, help =
         "The number of CPUs to use for parallelizing the blasting\n"
         "[default 1, i.e. no parallelism]\n" )

    arg( '--blast_full', action='store_true', default=False, help =
         "If specified, tells blast to use the full dataset of universal proteins\n"
         "[default False, i.e. the small dataset of universal proteins is used]\n" )

    arg( '-v','--version', action='store_true', help=
         "Prints the current PhyloPhlAn version and exit\n" )

    return vars(p.parse_args())


def init():
    if download:
        for f in [ppa_fna,ppa_aln,ppa_xml,ppa_up2prots,ppa_ors2prots]:
            u = urllib2.urlopen( download_compressed(f)  )
            if not os.path.exists( f ):
                info("Downloading and extracting "+f+" ... ")
                with closing(tarfile.open( fileobj = StringIO(u.read()) )) as inp:
                    inp.extractall(path=os.path.dirname(f))
                info("Done!\n")
    else:
        for f in [ppa_fna, ppa_fna_40, ppa_tax, ppa_aln, ppa_xml, ppa_up2prots, ppa_ors2prots]:
            if not os.path.exists( f ):
                info("Extracting "+f+" ... ")
                with closing(tarfile.open( compressed(f)  )) as inp:
                    inp.extractall(path=os.path.dirname(f))
                info("Done!\n")

        for t,f in [ppa_alns]:
            if not os.path.exists( t ):
                info("Extracting "+f+" ... ")
                with closing(tarfile.open( f )) as inp:
                    inp.extractall(path=os.path.dirname(f))
                info("Done!\n")

    if not os.path.exists( ppa_wdb ):
        info("Generating "+ppa_wdb+" (usearch indexed DB)... ")
        sb.call( ["usearch","-quiet",
                  "--makewdb",ppa_fna,
                  "--output",ppa_wdb])
        info("Done!\n")


def clean_all():
    files = glob('data/*.txt')
    files += glob('data/*.faa')
    files += glob('data/*.xml')
    files += glob('data/*.wdb')
    files += glob(os.path.dirname(ppa_alns[1]) + '/*.txt')
    files += glob(os.path.dirname(ppa_alns[1]) + '/*.score')
    files += glob(os.path.dirname(ppa_alns[1]) + '/*.aln')

    for f in files:
        if f not in []:
            sb.call(['rm', '-f', f])


def clean_project( proj ):
    sb.call(["rm","-rf","data/"+proj])
    sb.call(["rm","-rf","output/"+proj])


def get_inputs(proj):
    inp_fol = "input/"+proj+"/"
    dat_fol = "data/"+proj+"/"

    if not os.path.isdir(inp_fol):
        error( "No "+proj+" folder found in input/" )

    files = list(os.listdir(inp_fol))
    faa_in = [os.path.splitext(l)[0] for l in files if os.path.splitext(l)[1] == '.faa']
    fna_in = [os.path.splitext(l)[0] for l in files
              if (os.path.splitext(l)[1] == '.fna') and (os.path.splitext(l)[0] not in faa_in)]

    if len(faa_in) + len(fna_in) < 1:
        error("No '.faa' or '.fna' input files found in " + str(inp_fol))

    txt_in = [l for l in files if os.path.splitext(l)[1] == ".txt"]
    if len( txt_in ) > 1:
        error( "More than one '.txt' input files found in input/\n"
              "[No more than one txt file (the taxonomy) allowed" )
    tax_in = [l for l in files if os.path.splitext(l)[1] == ".tax"]
    if len( tax_in ) > 1:
        error( "More than one '.tax' input files found in input/\n"
              "[No more than one txt file (the taxonomy for taxonomic analysis) allowed" )
    mdt_in = [l for l in files if os.path.splitext(l)[1] == ".metadata"]
    if len( tax_in ) > 1:
        error( "More than one '.metadata' input files found in input/\n"
              "[No more than one txt file (the taxonomy for taxonomic analysis) allowed" )

    if not os.path.isdir(dat_fol):
        os.mkdir( dat_fol )

    if not os.path.exists(dat_fol + ors2prots):
        with open(dat_fol + ors2prots, 'w') as outf:
            for f in faa_in:
                prots = sorted([l.split()[0][1:] for l in open(inp_fol + f + '.faa') if l.startswith('>')])
                outf.write('\t'.join([f] + prots) + '\n')

            for f in fna_in:
                prots = sorted([l.split()[0][1:] for l in open(inp_fol + f + '.fna') if l.startswith('>')])
                outf.write('\t'.join([f] + prots) + '\n')

    return faa_in, fna_in, txt_in[0] if txt_in else None, tax_in[0] if tax_in else None, mdt_in[0] if mdt_in else None


def check_inp(inps):
    init_dig = [l for l in inps if l[0].isdigit()]

    if init_dig:
        error("The following genome IDs start with a digit\n" + "\n".join(init_dig) +
             "\nPlease rename accordingly the corresponding input file names")

    ppa_ids = [l[1:] for l in open(ppa_aln) if l.startswith('>')]
    init_dup = [l for l in inps if l in ppa_ids]

    if init_dup:
        error("The following genome IDs are already in PhyloPhlAn\n" + "\n".join(init_dig))


def exe_usearch(x):

    def screen_usearch_wdb( inpf ):
        tab = ((l[0].split(' ')[0],l[1].split('_')[1],float(l[-1]))
                for l in (ll.strip().split('\t')
                    for ll in open(inpf) ) )
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
        with open(inpf,"w") as outf:
            for k,v in f2p.items():
                outf.write( "\t".join([str(k),str(v)]) +"\n" )

    try:
        info( "Starting "+x[5] + "...\n" )
        sb.call( x )
        screen_usearch_wdb( x[5] )
        info( x[5] + " generated!\n" )
    except OSError:
        error( "OSError: fatal error running usearch." )
        return
    except ValueError:
        error( "ValueError: fatal error running usearch." )
        return
    except KeyboardInterrupt:
        error( "KeyboardInterrupt: usearch process interrupted." )
        return


def faa2ppafaa( inps, nproc, proj ):
    inp_fol = "input/"+proj+"/"
    dat_fol = "data/"+proj+"/usearch/"
    pool = mp.Pool( nproc )
    mmap = [(inp_fol+i+'.faa', dat_fol+i+'.b6o') for i in inps if not os.path.exists(dat_fol+i+'.b6o')]

    if not os.path.isdir(dat_fol): os.mkdir(dat_fol) # create the tmp directory if does not exists

    if not mmap:
        info("All usearch runs already performed!\n")
    else:
        info("Looking for PhyloPhlAn proteins in input faa files\n")
        us_cmd = [ ["usearch","-quiet",
                    "-wdb",ppa_wdb,
                    "-blast6out",o,
                    "-query",i,
                    "-evalue","1e-40"] for i,o in mmap ]
        pool.map_async( exe_usearch, us_cmd )
        pool.close()
        pool.join()
        info("All usearch runs performed!\n")

    if os.path.exists(dat_fol+up2prots):
        return
    up2p = collections.defaultdict( list )
    if not os.path.exists(dat_fol+up2prots):
        for i in inps:
            for p,up in (l.strip().split('\t')
                    for l in open(dat_fol+i+'.b6o').readlines()):
                up2p[up].append(p)
        with open(dat_fol+up2prots, 'w') as outf:
            for k,v in up2p.items():
                outf.write( "\t".join([k]+v) + "\n" )


def blastx_exe(x):
    try:
        info("Starting " + x[6] + "...\n")
        sb.call(x)
        info(x[6] + " generated!\n")
    except OSError:
        error("OSError: fatal error running tblastn.")
        return
    except ValueError:
        error("ValueError: fatal error running tblastn.")
        return
    except KeyboardInterrupt:
        error("KeyboardInterrupt: tblastn process interrupted.")
        return


def blast(inps, nproc, proj, blast_full=False):
    inp_fol = 'input/' + proj + '/'
    dat_fol = 'data/' + proj + '/tblastn/'
    pool = mp.Pool(nproc)
    mmap = [(inp_fol+i+'.fna', dat_fol+i+'.b6o') for i in inps if not os.path.exists(dat_fol+i+'.b6o')]

    if not os.path.isdir(dat_fol): os.mkdir(dat_fol) # create the tmp directory if does not exists

    if not mmap:
        info('All tblastn runs already performed!\n')
    else:
        info('Looking for PhyloPhlAn proteins in input fna files\n')

        if blast_full: # full dataset
            dataset = ppa_fna
        else: # small dataset
            dataset = ppa_fna_40

        us_cmd = [['tblastn', '-subject', i, '-query', dataset, '-out', o, '-outfmt', '6', '-evalue', '1e-40']
                  for i, o in mmap]

        pool.map_async(blastx_exe, us_cmd)
        pool.close()
        pool.join()
        info('All tblastn runs performed!\n')

    if not os.path.exists(dat_fol+up2prots):
        dic = collections.defaultdict(list)
        uniq = set()

        for i in inps:
            dicc = collections.defaultdict(list)

            # for each universal protein found tuples
            for lst in (l.strip().split('\t') for l in open(dat_fol+i+'.b6o').readlines()):
                upid = lst[0].split('_')[1]
                contig_id = str(lst[1])
                sstart = int(lst[8])
                send = int(lst[9])
                bitscore = float(lst[-1])

                if (contig_id, sstart, send) not in uniq:
                    uniq.add((contig_id, sstart, send))

                    if upid not in dicc:
                        dicc[upid] = lst
                    elif bitscore > float(dicc[upid][-1]):
                        dicc[upid] = lst

            dic[i] = dicc

        updic = collections.defaultdict(list)
        p2t = dict()

        # open each file
        for k, v in dic.iteritems():
            with open(inp_fol+k+'.fna', 'rU') as f:
                # parse the fasta (nucleotides) file
                for record in SeqIO.parse(f, 'fasta'):
                    # look for the contigs id
                    for kk, vv in v.iteritems():
                        if str(vv[1]) in record.id:
                            sstart = int(vv[8])
                            send = int(vv[9])
                            reverse = False

                            if sstart > send:
                                sstart, send = send, sstart
                                reverse = True

                            sequence = Seq(str(record.seq)[sstart-1:send])
                            if reverse: sequence = sequence.reverse_complement()
                            rev = ':c' if reverse else ':' # reverse or not
                            seqid = record.id+rev+str(sstart)+'-'+str(send)
                            aminoacids = Seq.translate(sequence)
                            updic[kk].append(SeqRecord(aminoacids, id=seqid, description=''))
                            p2t[seqid] = k # save map from seqid to genome

        # write the partial mapping of proteins into genomes
        with open('data/'+proj+'/'+p2t_map, 'w') as f:
            for k, v in p2t.iteritems():
                f.write(str(k) + '\t' + str(v) + '\n')

        # write mapping file dat_fol+up2prots
        with open(dat_fol+up2prots, 'w') as f:
            for k in sorted(updic):
                f.write('\t'.join([k] + [sr.id for sr in updic[k]]) + '\n')

        # write a fasta file for each up
        for k, v in updic.iteritems():
            with open(dat_fol+k+'.faa', 'w') as f:
                SeqIO.write(v, f, 'fasta')

        # write pickle file
        with open(dat_fol+ups2faa_pkl, 'w') as f:
            pickle.dump([k for k, _ in updic.iteritems()], f, pickle.HIGHEST_PROTOCOL)


def gens2prots(inps, proj ):
    inp_fol = "input/"+proj+"/"
    dat_fol = "data/"+proj+"/usearch/"

    if not os.path.isdir(dat_fol): os.mkdir(dat_fol) # create the tmp directory if does not exists

    if os.path.exists( dat_fol+ups2faa_pkl ):
        return

    ups2prots = dict([ (ll[0],set([i for i in ll[1:]])) for ll in
                            (l.strip().split('\t')
                                for l in open(dat_fol+up2prots))])

    prots2ups = {}
    for k,v in ups2prots.items():
        for vv in v:
            if vv in prots2ups:
                error(str(vv) + " already in dict!")
            prots2ups[vv] = k
    ups2faa = collections.defaultdict( set )
    e = None
    for i in inps:
        with open(inp_fol+i+".faa", 'U') as inp:
            for l in inp:
                if l.startswith(">"):
                    if e in prots2ups:
                        ups2faa[prots2ups[e]].add(s+"\n")
                    e, s = l.strip().split()[0][1:], ""
                s += l
            if e in prots2ups:
                ups2faa[prots2ups[e]].add(s+"\n")

    for k,v in ups2faa.items():
        with open(dat_fol+k+'.faa', 'w') as outf:
            for vv in v:
                outf.write( "".join(vv) )
    with open(dat_fol+ups2faa_pkl, 'w') as outf:
        pickle.dump(ups2faa, outf, pickle.HIGHEST_PROTOCOL)


def screen( stat, cols, sf = None, unknown_fraction = 0.5, n = 10 ):
    lena, nsc = len(cols[0]), []
    if sf and os.path.exists(sf):
        with open(sf) as sfinp:
            scores = [(float(ll[:12]),ll[12:]) for ll in [l for l in sfinp.readlines()]]
    else:
        scores = [(1.0,c) for c in cols]
    for sc,st in zip(scores,stat):
        if len(st) == 1:
            continue
        if "-" in st and len(st) == 2:
            continue
        if "X" in st and len(st) == 2:
            continue
        if "X" in st and "-" in st and len(st) == 3:
            continue
        nn = 0
        for k in ["-","X"]:
            if k in st:
                nn += st[k]
        if nn > lena*unknown_fraction:
            continue
        nsc.append( sc )
    nsc = sorted(nsc,key=lambda x:x[0],reverse = True)
    ret = [v for s,v in nsc][-n:]
    return ret if ret else [list(['-']*lena+['\n'])]


def aln_subsample( inp_f, out_f, scores, unknown_fraction, namn ):
    with open(inp_f) as inp:
        fnas = list(SeqIO.parse(inp, "fasta"))
    ids = [f.id for f in fnas]
    cols = zip(*[f.seq for f in fnas])
    # l = len(cols[0])
    col_alf = [set(x) for x in cols]
    col_sta = [dict([(a,seq.count(a)) for a in s]) for s,seq in zip(col_alf,cols)]
    col_sta_ok = screen(col_sta,cols, sf = scores,
                        unknown_fraction = unknown_fraction, n = namn )
    raws = zip(*col_sta_ok)
    recs = [SeqRecord(seq=Seq("".join(r)),id=nid,description=nid)
                for r,nid in zip(raws,ids)]
    with open(out_f,"w") as out:
        SeqIO.write( recs, out, "fasta")


def exe_muscle(x):
    try:
        if len(x) < 13:
            info( "Running muscle on "+x[3] + "...\n" )
        else:
            info( "Running muscle on "+x[4] + " and "+x[6]+"...\n" )
        sb.call( x[:-2] )
        pn = max( int( max( int((400.0-int(x[-1][1:]))*30.0/400.0),1)**2 / 30.0), 3)
        if len(x) < 13:
            aln_subsample( x[5], x[-2], x[7], 0.1, pn )
            info( x[-2] + " generated (from "+x[3]+")!\n" )
        else:
            aln_subsample( x[8], x[-2], x[10], 0.1, pn )
            info( x[-2] + " generated (from "+x[4]+" and "+x[6]+")!\n" )

    except OSError, e:
        error( "OSError: fatal error running muscle." )
        raise e
    except ValueError, e:
        error( "ValueError: fatal error running muscle." )
        raise e
    except KeyboardInterrupt, e:
        error( "KeyboardInterrupt: usearch process muscle." )
        raise e
    except Exception, e:
        error( e )
        raise e


def faa2aln( nproc, proj, integrate = False ):
    dat_fol = "data/"+proj+"/"
    with open( dat_fol + ups2faa_pkl ) as inp:
        prots = pickle.load(inp)
    mmap = [(   dat_fol+i+".faa",dat_fol+i+".aln",dat_fol+i+".sc",
                ppa_alns_fol+i+".aln",ppa_alns_fol+i+".aln.score",
                dat_fol+i+".int.aln",dat_fol+i+".int.sc",
                dat_fol+i+".sub.aln",dat_fol+i+".int.sub.aln",i in prots,
                i)
                for i in ('p{num:04d}'.format(num=aa) for aa in range(400))
                        ]
    if mmap:
        us_cmd = [ ["muscle","-quiet",
                    "-in",i,
                    "-out",o,
                    "-scorefile",s,
                    "-maxiters","2",
                    so,pn] for i,o,s,o2,s2,oi,si,so,soi,pres,pn in mmap
                        if not os.path.exists(o) and pres]
        if us_cmd:
            info("Looking for PhyloPhlAn proteins to align\n")
            info(str(len(us_cmd))+" alignments to be performed\n")
            pool = mp.Pool( nproc )
            pool.map_async( exe_muscle, us_cmd )
            pool.close()
            pool.join()
            info("All alignments performed!\n")
    if integrate:
        us_cmd = [ ["muscle","-quiet",
                    "-profile",
                    "-in1",o,
                    "-in2",o2,
                    "-out",oi,
                    "-scorefile",si,
                    "-maxiters","2",
                    soi,pn] for i,o,s,o2,s2,oi,si,so,soi,pres,pn in mmap
                        if not os.path.exists(soi) and pres]
        for i,o,s,o2,s2,oi,si,so,soi,pres,pn in mmap: #paralellizable
            if not os.path.exists(soi) and not pres:
                pnn = max( int( max( int((400.0-int(pn[1:]))*30.0/400.0),1)**2 / 30.0), 3)
                aln_subsample( o2, soi, s2, 0.1, pnn )
        if us_cmd:
            info("Merging alignments from user genomes with all the PhyloPhlAn alignments\n")
            info(str(len(us_cmd))+" alignments to be merged\n")
            pool = mp.Pool( nproc )
            pool.map_async( exe_muscle, us_cmd )
            pool.close()
            pool.join()
            info("All alignments already merged with PhyloPhlAn alignments!\n")

    if os.path.exists( dat_fol+up2prots ):
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
    dat_fol = "data/"+proj+"/"
    outf = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot

    if os.path.exists( outf ):
        info("All alignments already merged!\n")
        return

    up2p = dict([(l[0],set(l[1:])) for l in
                    (ll.strip().split('\t') for ll in
                        open(dat_fol+up2prots))])

    all_prots = set.union(*up2p.values()) # all proteins id

    if integrate:
        for l in (ll.strip().split('\t') for ll in open(ppa_up2prots)):
            if l[0] in up2p:
                up2p[l[0]] |= set(l[1:])
            else:
                up2p[l[0]] = set(l[1:])

    t2p = dict([(l[0],set(l[1:])) for l in
                (ll.strip().split('\t') for ll in open(dat_fol+ors2prots))])
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
    aln = dict([(t,"") for t in t2p.keys()]) # dictionary that will contains all alignments

    for f in sorted((dat_fol+ff for ff in os.listdir(dat_fol) if (ff.endswith(".int.sub.aln") and integrate) or (not integrate and ff.endswith(".sub.aln") and not ff.endswith(".int.sub.aln")))):
        g2aln = dict()
        lenal = 0 # keep the lenght of the current alignments

        for seq in SeqIO.parse(f, "fasta"):
            if not lenal: lenal = len(seq.seq)
            g2aln[p2t[seq.id]] = seq

        for g in all_g: # for all genomes
            if g in g2aln: # if there is alignment
                aln[g] += g2aln[g]
            else: # otherwise add gaps
                aln[g] += SeqRecord(Seq('-' * lenal), id=str(g))

###########################################################
    # faas = []
    # faas = SeqIO.to_dict( faas )
    # aln = dict([(t,"") for t in t2p.keys()])
    # salnk = sorted(aln.keys())

    # for up,pts in sorted(up2p.items(),key=lambda x:x[0]):
    #     toadd, ll = [], -1
    #     for k in salnk:
    #         ppp = []
    #         for x in list(t2p[k]):
    #             for y in list(pts):
    #                 if str(x) in str(y):
    #                     ppp += [y]
    #         if ppp:
    #             if ll < 0:
    #                 ll = len( faas[ppp[0]] )
    #             aln[k] += faas[ppp[0]]
    #         else:
    #             toadd.append(k)
    #     for k in toadd:
    #         ss = SeqRecord( Seq("".join( '-' * ll )) )
    #         ss.id = k
    #         aln[k] += ss
###########################################################

    out_faas = []
    for k, v in aln.iteritems():
        v.id = k
        v.description = ""
        out_faas.append(v)

    SeqIO.write(out_faas, outf, "fasta")
    info("All alignments merged into "+outf+"!\n")


def fasttree( proj, integrate ):
    out_fol = "output/"+proj+"/"
    dat_fol = "data/"+proj+"/"
    outt = out_fol+proj+"."+ ( o_inttree if integrate else o_tree)
    if not os.path.isdir(out_fol):
        os.mkdir( out_fol )
    aln_in = dat_fol+aln_int_tot if integrate else dat_fol+aln_tot
    if os.path.exists( outt ):
        info("Final tree already built ("+outt+")!\n")
        return
    info("Start building the tree with FastTree ... \n")
    cmd = [ "FastTree", "-quiet",
            #"-fastest","-noml"
            #"-gamma",
            "-bionj","-slownni",
            "-mlacc","2", # "-pseudo",
            "-spr","4"
            ]
    sb.call( cmd,
             stdin =  open(aln_in),
             stdout = open(outt,"w") )

    # convert aln_in from fasta to phylip format
    #cmd = [ "raxmlHPC-PTHREADS-SSE3",
    #        "-T","#threads",
    #        "-s","alignment in phylip format",
    #        "-n",outt,
    #        "-m","PROTCATWAG",
    #        "-p","1982"
    #        ]
    #sb.call( cmd )
    info("Tree built! The output newick file is in "+outt+"\n")


def circlader( proj, integrate, tax = None ):
    inp_fol = "input/"+proj+"/"
    dat_fol = "data/"+proj+"/"
    out_fol = "output/"+proj+"/"
    # out_img = "output/"+proj+"/imgs/"

    tree = out_fol+proj+"."+ ( o_inttree if integrate else o_tree)
    xtree_rr = out_fol+proj+".tree.reroot."+ ( "int." if integrate else "")+"xml"
    xtree_an = out_fol+proj+".tree.reroot."+ ( "int." if integrate else "")+"annot.xml"
    pngtree = out_fol+proj+".tree"+ ( ".int" if integrate else "")+".png"
    annotation_f = dat_fol+proj+( ".int" if integrate else "")+".annot.txt"
    archaea_f = dat_fol+proj+( ".int" if integrate else "")+".archaea.txt"

    tax_d = [l.strip().split('\t') for l in open(inp_fol+tax)] if tax else []
    # o2t = dict([(o,t.split('.')) for t,o in tax_d])
    t2o = dict(tax_d)

    if os.path.exists( pngtree ):
        info("Output tree images already generated!\n")
        return

    #info("Generating tree output images ...")
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
        with open( archaea_f, "w" ) as out:
            out.write( "\n".join( mt4archaea ) +"\n" )
        sb.call( [  c_reroot, "-s", "lca", "-f", archaea_f,
                    tree, xtree_rr ] )

        annotation_f = dat_fol+proj+".annot.txt"
        with open( annotation_f, "w") as outf:
            for v1,v2 in [('clade_fill_color','r')]:
                outf.write( "\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n" )
            for v1,v2 in [  ('clade_size','7.0'),
                         ]:
                outf.write( "\n".join(["\t".join([nn,v1,v2]) for nn in names+int_names]) + "\n" )
            for v1,v2 in [('clade_size','2.0')]:
                outf.write( "\t".join(["*",v1,v2]) + "\n" )
            for v1,v2 in [  ('sub_branches_angle_reduction','0.0'),
                            ('clade_edge_line_width','0.2'),
                            ('branch_line_width','0.45'),
                         ]:
                outf.write( "\t".join([v1,v2]) + "\n" )

    else:
        n = 1 if len(names) <  4 else int(round(min(len(names)*0.05,100)))
        archaea = [o for t,o in t2o.items() if t.split(".")[0] == "d__Archaea"]
        if t2o and archaea:
            with open( archaea_f, "w" ) as out:
                out.write( "\n".join( archaea ) +"\n" )
            sb.call( [  c_reroot, "-s", "lca", "-f", archaea_f,
                        tree, xtree_rr ] )
        elif n == 1 and len(names) < 8:
            sb.call( [ c_reroot, "-s","longest_edge", tree, xtree_rr ] )
        elif n == 1:
            sb.call( [  c_reroot, "-s", "longest_internal_edge",
                        tree, xtree_rr ] )
        else:
            sb.call( [  c_reroot, "-s", "longest_internal_edge_n", "-n",
                        str(n), tree, xtree_rr ] )

        return
        with open( annotation_f, "w") as outf:
            for v1,v2 in [  ('clade_fill_color','r'),
                            ('clade_size','7.0'),
                         ]:
                outf.write( "\n".join(["\t".join([nn,v1,v2]) for nn in names]) + "\n" )
            for v1,v2 in [('clade_size','2.0')]:
                outf.write( "\t".join(["*",v1,v2]) + "\n" )
            for v1,v2 in [  ('sub_branches_angle_reduction','0.0'),
                            ('clade_edge_line_width','0.2'),
                            ('branch_line_width','0.45'),
                         ]:
                outf.write( "\t".join([v1,v2]) + "\n" )

    sb.call( [ c_circ_an, "--annot", annotation_f, xtree_rr, xtree_an ] )
    sb.call( [ c_circ, "--dpi", "300", xtree_an, pngtree ] )
    info( "Done!\n" )


def tax_curation( proj, tax, integrate = False, mtdt = None, inps = None ):
    inp_fol = "input/"+proj+"/"
    out_fol = "output/"+proj+"/"
    taxin = ppa_tax if integrate else inp_fol+tax
    taxout = None if integrate else out_fol+"new."+tax
    taxdiffs = None if integrate else out_fol+"diff."+tax
    taxboth = None if integrate else out_fol+"comp."+tax
    intree = out_fol+proj+".tree.reroot."+ ( "int." if integrate else "")+"annot.xml"
    if not os.path.exists(intree):
        intree = out_fol+proj+".tree."+ ( "int." if integrate else "")+"nwk"
    operation = "imputation" if integrate else "curation"

    if tax:
        info("Starting taxonomic curation: detecting suspicious taxonomic labels ... ")
        taxCleaner = taxc.TaxCleaner( taxin, intree, integrate = integrate, to_integrate = inp_fol+tax )
        taxCleaner.remove_tlabels()
        info("Done!\n")
        info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy ... ")
    else:
        taxCleaner = taxc.TaxCleaner( taxin, intree, integrate = integrate, inps = inps  )
        info("Trying to impute taxonomic labels for taxa newly integrated into the tree... ")

    taxCleaner.infer_tlabels()
    info("Done!\n")

    info("Writing taxonomic "+operation+" outputs ... ")
    tid2img = taxCleaner.write_output( proj, images = False )
    info("Done!\n")
    info("Writing final pdf reports ... ")
    if integrate:
        pass
    elif mtdt:
        taxCleaner.write_report( out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images = True, typ = "refined" )
        taxCleaner.write_report( out_fol+"/"+"tax"+operation+"_report", tid2img, inp_fol+mtdt, images = True, typ = "corrected" )

    if not integrate:
        taxCleaner.write_new_taxonomy( taxout, taxdiffs, taxboth )

    info("Done!\n")


def tax_imputation( proj, tax, integrate = False, mtdt = None, inps = None ):
    inp_fol = "input/"+proj+"/"
    out_fol = "output/"+proj+"/"
    taxin = ppa_tax if integrate else inp_fol+tax
    # taxout = None if integrate else out_fol+"new."+tax
    # taxdiffs = None if integrate else out_fol+"diff."+tax
    # taxboth = None if integrate else out_fol+"comp."+tax
    intree = out_fol+proj+".tree.reroot."+ ( "int." if integrate else "")+"annot.xml"
    if not os.path.exists(intree):
        intree = out_fol+proj+".tree."+ ( "int." if integrate else "")+"nwk"
    operation = "imputation" if integrate else "curation"

    info("Trying to impute taxonomic labels for taxa newly integrated into the tree... ")
    taxCleaner = taxc.TaxCleaner( taxin, intree, integrate = integrate, inps = inps  )

    taxCleaner.infer_tlabels()
    info("Done!\n")

    info("Writing taxonomic "+operation+" outputs ... ")
    taxCleaner.write_output( proj, images = False, imp_only = True )
    info("Done!\n")


def tax_curation_test( proj, tax,
                       nerrorrs = 10, error_type = 'wrong',
                       taxl = 's', tmin = 3, tex = None,
                       name = 'err.txt',
                       integrate = False, descr = None ):
    inp_fol = "input/"+proj+"/"
    out_fol = "output/"+proj+"/"
    err_tax = out_fol+name
    intree = out_fol+proj+".tree.reroot."+ ( "int." if integrate else "")+"annot.xml"


    taxCleaner = taxc.TaxCleaner( inp_fol+tax, intree )
    taxCleaner.introduce_errors( nerrorrs, err_tax, error_type, taxl, tmin, tex )
    taxCleaner = taxc.TaxCleaner( err_tax, intree )

    info("Starting taxonomic curation: detecting suspicious taxonomic labels ... ")
    taxCleaner.remove_tlabels()
    info("Done!\n")
    info("Trying to impute taxonomic labels for taxa with suspicious or incomplete taxonomy ... ")
    taxCleaner.infer_tlabels()
    info("Done!\n")

    info("Writing taxonomic curation outputs and evaluations... ")
    taxCleaner.write_output( "log.txt", proj, outname = name ) # , images = True )
    taxCleaner.evaluate( proj, inp_fol+tax, err_tax, name, descr = descr )
    info("Done!\n")


def merge_usearch_blast(inps, proj):
    dat_fol = 'data/'+proj+'/'
    usearch_fol = 'data/'+proj+'/usearch/'
    tblastn_fol = 'data/'+proj+'/tblastn/'
    usearch_files = []
    tblastn_files = []

    if os.path.isdir(usearch_fol):
        usearch_files = [os.path.normcase(f) for f in os.listdir(usearch_fol)]

    if os.path.isdir(tblastn_fol):
        tblastn_files = [os.path.normcase(f) for f in os.listdir(tblastn_fol)]

    # check if I have already copied files in dat_fol folder
    dat_fol_num = len(os.listdir(dat_fol))
    if usearch_files and tblastn_files:
        dat_fol_num -= 4 # bad!
    elif (usearch_files and not tblastn_files) or (tblastn_files and not usearch_files):
        dat_fol_num -= 2 # bad!

    if dat_fol_num:
        return

    if usearch_files and not tblastn_files: # just usearch ran -> NO MERGE, just copy
        for f in usearch_files:
            shutil.copy2(usearch_fol+f, dat_fol+f)
    elif tblastn_files and not usearch_files: # just tblastn ran -> NO MERGE, just copy
        for f in tblastn_files:
            shutil.copy2(tblastn_fol+f, dat_fol+f)
    elif usearch_files and tblastn_files: # usearch and tblastn ran -> MERGE!
        # copy .b6o files
        for f in [ff for ff in usearch_files if '.b6o' in ff]:
            shutil.copy2(usearch_fol+f, dat_fol+f)

        for f in [ff for ff in tblastn_files if '.b6o' in ff]:
            shutil.copy2(tblastn_fol+f, dat_fol+f)

        # merge .faa files
        sett = list(set([f for f in usearch_files if '.faa' in f]) & set([f for f in tblastn_files if '.faa' in f]))

        for f in sett:
            records = []
            with open(usearch_fol+f, 'rU') as ff:
                records = list(SeqIO.parse(ff, 'fasta'))

            with open(tblastn_fol+f, 'rU') as ff:
                records += list(SeqIO.parse(ff, 'fasta'))

            with open(dat_fol+f, 'w') as ff:
                SeqIO.write(records, ff, 'fasta')

        for f in usearch_files:
            if f not in sett:
                shutil.copy2(usearch_fol+f, dat_fol+f)

        for f in tblastn_files:
            if f not in sett:
                shutil.copy2(tblastn_fol+f, dat_fol+f)

        # merge up2prots files
        up2prot = collections.defaultdict(list)
        with open(usearch_fol+up2prots, 'r') as f:
            for row in f:
                up = row.split('\t')[0].strip()
                prot = [s.strip() for s in row.split('\t')[1:]]
                up2prot[up] = prot

        with open(tblastn_fol+up2prots, 'r') as f:
            for row in f:
                up = row.split('\t')[0].strip()
                prot = [s.strip() for s in row.split('\t')[1:]]

                if up in up2prot:
                    up2prot[up] += prot
                else:
                    up2prot[up] = prot

        with open(dat_fol+up2prots, 'w') as f:
            for k, v in up2prot.iteritems():
                f.write('\t'.join([k] + v) + '\n')

        # merge ups2faa_pkl files
        with open(usearch_fol+ups2faa_pkl, 'r') as f:
            pckl = set(pickle.load(f))

        with open(tblastn_fol+ups2faa_pkl, 'r') as f:
            pckl |= set(pickle.load(f))

        with open(dat_fol+ups2faa_pkl, 'w') as f:
            pickle.dump(list(pckl), f, pickle.HIGHEST_PROTOCOL)
    else: # something very bad happened!
        error("Merge phase, no files found!")


if __name__ == '__main__':
    pars = read_params( sys.argv )
    projn = pars['inp']

    if pars['cleanall']:
        if ('taxonomic_analysis' in pars and pars['taxonomic_analysis']) or (
                'user_tree' in pars and pars['user_tree']) or (
                        'integrate' in pars and pars['integrate']):
            error("--cleanall is in conflict with -t, -u, and -i")
        else:
            clean_all()
        sys.exit(0)

    if pars['clean']:
        if ('taxonomic_analysis' in pars and pars['taxonomic_analysis']) or (
                'user_tree' in pars and pars['user_tree']) or (
                        'integrate' in pars and pars['integrate']):
            error("--cleanall is in conflict with -t, -u, and -i")
        else:
            clean_project(projn)
        sys.exit(0)

    if 'v' in pars and pars['v']:
        sys.stdout.write("PhyloPhlAn version "+__version__+" ("+__date__+")\n")
        sys.exit(0)
    if projn == None:
        error("Project name not provided.")

    dep_checks()
    init()
    faa_in, fna_in, tax, rtax, mtdt = get_inputs(projn)

    if not tax and pars['taxonomic_analysis'] and not pars['integrate']:
        error("No taxonomy file found for the taxonomic analysis")

    check_inp(faa_in+fna_in)
    t0 = time.time()

    if (pars['nproc'] > 1) and faa_in and fna_in:
        cores = float(pars['nproc'])
        tasks = len(faa_in) + len(fna_in)
        usearch_cpus = int(round(len(faa_in) * (cores/tasks))) + 1
        blast_cpus = int(round(len(fna_in) * (cores/tasks))) + 1

     # balancing cpu load
        while (usearch_cpus + blast_cpus) > cores:
            if usearch_cpus > 1:
                usearch_cpus -= 1
            else:
                blast_cpus -= 1

        if (blast_cpus > len(fna_in)) and (usearch_cpus < len(faa_in)):
            surplus = blast_cpus - len(fna_in)
            blast_cpus -= surplus
            usearch_cpus += surplus

        if (usearch_cpus > len(faa_in)) and (blast_cpus < len(fna_in)):
            surplus = usearch_cpus - len(faa_in)
            usearch_cpus -= surplus
            blast_cpus += surplus
     # balancing cpu load

        usearch_thr = mp.Process(target=faa2ppafaa, args=(faa_in, usearch_cpus, projn))
        usearch_thr.start()

        blast_thr = mp.Process(target=blast, args=(fna_in, blast_cpus, projn, pars['blast_full']))
        blast_thr.start()

        usearch_thr.join()
        blast_thr.join()

        gens2prots(faa_in, projn)
    else:
        if fna_in:
            blast(fna_in, pars['nproc'], projn, pars['blast_full'])

        if faa_in:
            faa2ppafaa(faa_in, pars['nproc'], projn)
            gens2prots(faa_in, projn)

    # merging phase for .faa and .fna files
    merge_usearch_blast(faa_in+fna_in, projn)

    t1 = time.time()
    sys.stdout.write("Mapping finished in "+str(int(t1-t0))+" secs.\n\n" )

    faa2aln( pars['nproc'], projn, pars['integrate'] )

    t2 = time.time()
    sys.stdout.write("Aligning finished in "+str(int(t2-t1))+" secs ("+str(int(t2-t0))+" total time).\n\n" )

    aln_merge(projn, pars['integrate'])
    fasttree(projn,pars['integrate'])

    t4 = time.time()
    sys.stdout.write("Tree building finished in "+str(int(t4-t2))+" secs ("+str(int(t4-t0))+" total time).\n\n" )

    circlader(projn, pars['integrate'], tax)

    #if pars['integrate'] and tax and not pars['taxonomic_analysis']:
    #    tax_imputation( projn, tax, mtdt = mtdt, integrate = pars['integrate'], inps = inps )
    if pars['integrate'] and pars['taxonomic_analysis']:
        tax_imputation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)
        sys.exit()

    if not pars['taxonomic_analysis']:
        sys.exit()

    if pars['tax_test']:
        nerrorrs, error_type, taxl, tmin, tex, name = pars['tax_test'].split(":")
        tax_curation_test( projn, tax,
                           nerrorrs = int(nerrorrs),
                           error_type = error_type,
                           taxl = taxl,
                           tmin = int(tmin),
                           tex = int(tex) if tex else None,
                           name = name,
                           descr = pars['tax_test'] )
    else:
        tax_curation(projn, tax, mtdt=mtdt, integrate=pars['integrate'], inps=faa_in+fna_in)
