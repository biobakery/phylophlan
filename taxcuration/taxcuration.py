import sys
import os
import copy 
import subprocess as sb 
import collections
import tempfile as tf
import random as rnd
from scipy import stats

# circlader_lib to be updated!!!!
import pyphlan as circ
import pyphlan as ppa

col_list = ["#0000FF","#800080","#FF0000","#800000","#FFFF00","#808000", 
            "#00FF00","#008000","#00FFFF","#008080","#C0C0C0","#000080",
            "#FF00FF","#808080","#FFC0CB","#FFA500","#FF6347","#F0E68C",
            "#9400D3","#7B68EE","#9ACD32","#3CB371","#40E0D0","#4682B4",
            "#00BFFF","#8B4513","#BC8F8F"]

def get_tmp_file(closed = False):
    tmpf = tf.NamedTemporaryFile(mode='w',delete=False)
    if closed:
        tmpf.close()
        return tmpf.name
    return tmpf

def annotate_nodes( l, nodes, fout ):
    if nodes is None:
        fout.write( "\n".join(["\t".join([v1,v2]) for v1,v2 in l]) + "\n" )
        return
    for v1,v2 in l:
        fout.write( "\n".join(["\t".join([nn,v1,v2]) for nn in nodes]) + "\n" )

class Taxon:

    sTaxLevs = 'dpcofgst'
    cTaxLevSep = '.'
    unkwn = '?'
    
    error_type = ['unclassified_tlab','spp','redundant_clade','ltcs_out']

    def __init__( self, tid = None, sTax = None ):
        self.names = dict([(c,"") for c in self.sTaxLevs])
        self.newnames = {} 
        self.torem = {}
        self.refined = {}
        self.corrected = {}
        self.new_sp = None
        if tid and sTax:
            self.load_taxon( tid, sTax)

    def load_taxon( self, tid, sTax ):
        self.tid = tid
        for v in sTax.split(self.cTaxLevSep):
            self.names[v[0]] = v
        self.analyze()

    def fta( self, tmax = None, kwn_only = False, newnames = False, viz = False ):
        nnames = self.newnames if newnames else self.names
        ret = []
        if tmax:
            tmax_i = self.sTaxLevs.find( tmax )
        for i,stmp in enumerate(self.sTaxLevs): 
            if tmax and i > tmax_i:
                continue
            if nnames[stmp] == self.unkwn:
                if kwn_only:
                    return None
                ret.append( stmp+"__"+self.unkwn )
            else:
                ret.append( nnames[stmp] )
        ret = self.cTaxLevSep.join( ret )
        if viz:
            ret = ret.replace("."," ").replace("__",":")
        return ret

    def remove_all( self ):
        for t in self.sTaxLevs[:-1]:
            self.newnames[t] = self.unkwn
        self.newnames['t'] = self.names['t']

    def copy2new( self, from_taxon, tmax = 's' ):
        self.newnames = {}
        copy = True
        for t in self.sTaxLevs[:-1]:
            if copy:
                self.newnames[t] = from_taxon.names[t]
            else:
                self.newnames[t] = self.unkwn 
            if t == tmax:
                copy = False 
        self.newnames['t'] = self.names['t']

    def get_changes( self ):
        ret, last = [], ""
        for t in self.sTaxLevs[:-1]:
            if self.names[t] == self.newnames[t] and not ret:
                last = t
                continue
            if not ret:
                ret.append( self.names[last] )
            ret.append( self.names[t] if self.names[t] != self.unkwn else t+"__"+self.unkwn )
        return self.cTaxLevSep.join( ret )

    def is_unchanged( self ):
        for t in self.sTaxLevs[:-1]:
            if self.names[t] != self.newnames[t]:
                return False
        return True

    def is_removed( self ):
        if len(self.newnames) == 0:
            return False
        for t in self.sTaxLevs[:-1]:
            if self.names[t] != self.newnames[t]:
                if self.newnames[t] == self.unkwn:
                    return True
        return False

    def is_refinement( self ):
        if len(self.newnames) == 0:
            return None
        for t in self.sTaxLevs[:-1]:
            if self.names[t] != self.newnames[t]:
                if self.names[t] == self.unkwn:
                    return True
                else:
                    return False
        return False

    def is_corrected( self ):
        for t in self.sTaxLevs[:-1]:
            if self.names[t] != self.newnames[t]:
                if self.names[t] != self.unkwn and self.newnames[t] != self.unkwn:
                    return True
                else:
                    return False
        return False


    def analyze( self ):
        uncl = ['unclassified','__bacterium','__proteobacterium','__cf_','__sp_','__?']
        self.uncl = 0
        self.glev = False
        for s in self.sTaxLevs:
            if any([unclc in self.names[s] for unclc in uncl]):
                self.names[s] = self.unkwn
                if self.uncl == 0 and s == 's':
                    self.glev = True
                self.uncl += 1
                self.torem['unclassified_tlab'] = {'val':2}


class TaxCleaner:

    def __init__( self, fTax, nTree, integrate = False, to_integrate = None, inps = None ):
        self.nTree = ppa.PpaTree(nTree)
        self.taxa = {}
        self.var = {}
        self.clades2terms = None 
        self.dists = None
        self.totbrlen = None 
        self.to_skip = set() 
        #self.cTree = circ.CircTree( nTree ) 
        with open(fTax) as inpf:
            for tid,tax in (s.strip().split('\t') for s in inpf):
                self.taxa[tid] = Taxon( tid = tid, sTax = tax )
        if integrate:
            self.to_skip = set( self.taxa.keys()  )
            if to_integrate:
                with open(to_integrate) as inpf:
                    for tid,tax in (s.strip().split('\t') for s in inpf):
                        self.taxa[tid] = Taxon( tid = tid, sTax = tax )
                        self.taxa[tid].torem = { 'unknown': {"val":2} } 
            if inps:
                for tid in inps:
                    self.taxa[tid] = Taxon( tid = tid, sTax = ".".join([ttlev+"__?" for ttlev in Taxon.sTaxLevs]) )
                    self.taxa[tid].torem = { 'unknown': {"val":2} }

    def introduce_errors( self, nerrors, outfile, error_type, taxl, tmin = 3, tex = None ):
        c2t = collections.defaultdict(list)
        f2t = {}
        for tid,t in self.taxa.items():
            tfta = t.fta( tmax = taxl, kwn_only = True )
            if tfta and not tfta.count("?") and not 'unclassified_tlab' in t.torem:
                c2t[ tfta ].append( tid )
        if tex:
            cc2t = dict([(k,v) for k,v in c2t.items() if len(v) == tex ])
        else:
            cc2t = dict([(k,v) for k,v in c2t.items() if len(v) >= tmin ])
       
        if error_type == 'unrefine':
            sel_cl = rnd.sample( cc2t.keys(), nerrors )
            for i,c in enumerate(sel_cl):
                rtid = rnd.choice(cc2t[c])
                fr = self.taxa[rtid].fta()
                self.taxa[rtid].names['s'] = Taxon.unkwn
                self.taxa[rtid].names['t'] = "t__synth_err_"+str(i)
                f2t[fr] = self.taxa[rtid].fta()

        if error_type == 'wrong':
            remm = {}
            while len(remm) < nerrors:
                sel_cl = rnd.choice( cc2t.keys() )
                sel_cl_po = Taxon.cTaxLevSep.join(sel_cl.split(Taxon.cTaxLevSep)[:-1])
                pot = [l for l in cc2t.keys() if l and l.count(sel_cl_po) and sel_cl != l] 
                if pot:
                    remm[sel_cl] = pot

            for i,c in enumerate(remm):
                rtid = rnd.choice(cc2t[c])
                toc = rnd.choice( remm[c] )
                fr = self.taxa[rtid].fta()
                to = self.taxa[ rnd.choice(cc2t[toc]) ]
                for k,v in self.taxa[rtid].names.items()[:-1]:
                    self.taxa[rtid].names[k] = to.names[k]
                self.taxa[rtid].names['t'] = "t__synth_err_"+str(i)
                f2t[fr] = self.taxa[rtid].fta()

        with open( outfile, "w") as outf:
            for k,v in self.taxa.items():
                outf.write("\t".join([k,v.fta()]) + "\n" )
        
        with open( outfile+".gt", "w") as outf:
            for k,v in f2t.items():
                outf.write("\t".join([k,v])+"\n")

    def lcca( self, t, c2t = None, parent_steps = 15, tlev = 's' ):
        # get the path from the root to leaf t
        node_path = self.nTree.tree.get_path(t.tid)
        if node_path:
            node_path = list(node_path)
        if not node_path or len(node_path) < 2:
            return None,None,None
        # for all nodes in the path starting from the highest ones....
        for ii,p in enumerate(node_path[-parent_steps:]):
            terms = c2t[p] if c2t else list(p.get_terminals())
            
            # extract al terminals different from t with their tax label 
            descn = [(term,val) for term,val in 
                        ((tt,self.taxa[tt.name].fta( tmax = tlev, kwn_only = True ))
                            for tt in terms if tt.name != t.tid and tt.name in self.taxa)
                        if val]
            
            if not descn or len(descn) < 2:
                return None,None,None

            terms,descn = zip(*descn)
            # if there's only one label in the subtree ...
            if len(set(descn)) == 1 and descn[0] != t.fta( tmax = tlev, kwn_only = True ):
                for c in p.clades:
                    ttterm = c2t[c] if c2t else list(c.get_terminals())
                    if set([at.name for at in terms]) <= set([tt.name for tt in ttterm]):
                        #return None,None,None
                        continue
                return p,terms,descn[0]
        return None,None,None

    def comp_variability( self ):
        if self.dists is None:
            self.dists = ppa.dist_matrix( self.nTree.tree ) 
       
        self.var = {}
        for tlev in Taxon.sTaxLevs[:-1]:
            tc2t = self.nTree.get_c2t()
            c2t = collections.defaultdict(list)
            for tid,taxon in self.taxa.items():
                sTaxon = taxon.fta(tmax = tlev, kwn_only = True)
                c2t[sTaxon].append( tid )
            dists = []
            for c,t in c2t.items():
                if len(t) < 2:
                    continue
                for t1 in t:
                    dists.append( min([self.dists[t1][t2] for t2 in t if t1 != t2]) )
            dists = sorted(dists)
            self.var[tlev] = {}
            for i in range(100):
                self.var[tlev][i] = dists[int(len(dists)*float(i)/100.0)] if dists else None

    def dist_rank( self, taxon, terms ):
        if self.dists is None:
            all_dists = [   min([self.nTree.tree.distance(l1.name,l2.name)
                                for l1 in terms
                                    if l1.name != l2.name and l1.name != taxon.tid])
                                for l2 in terms
                                    if l2.name != taxon.tid]
            mind = min([self.nTree.tree.distance(l.name,taxon.tid) for l in terms if l.name != taxon.tid])
        else:
            all_dists = [   min([self.dists[l1.name][l2.name]
                                for l1 in terms
                                    if l1.name != l2.name and l1.name != taxon.tid])
                                for l2 in terms
                                    if l2.name != taxon.tid]
            mind = min([self.dists[l.name][taxon.tid] for l in terms if l.name != taxon.tid])
        sdists = sorted(all_dists+[mind])
        r = float(sdists.index(mind))/float(len(sdists))
        return r

    def onevsall_dist( self, taxon, tlev = 's' ):
        terms = self.nTree.tree.get_terminals()
        if self.dists is None:
            pass
        else:
            all_dists = [(self.dists[taxon.tid][t.name],t.name) for t in terms if t.name != taxon.tid and len(self.taxa[t.name].torem) == 0]
        all_dists = sorted( all_dists, key = lambda x:x[0] )
        ret = [all_dists[0]]
        cl_lab = self.taxa[all_dists[0][1]].fta( tmax = tlev, kwn_only = True )
        if cl_lab is None:
            return None
        for d,t in all_dists[1:]:
            if self.taxa[t].fta( tmax = tlev, kwn_only = True ) != cl_lab:
                break
            ret.append((d,t))
        return ret


    def save_circ_img( self, taxon, cortype, outdir, tlev = 'g' ):
        c_circlader = "src/circlader/circlader.py"
        c_annotate = "src/circlader/circlader_annotate.py" 

        tree_tf = get_tmp_file() #
        tbl = 0.0
        z,a,b,d = '', 'g', 'f', ' no sp'
        if tlev == 'f':
            z,a,b,d = 'g','f','o', ' no genus'
        if tlev == 'o':
            z,a,b,d = 'f','o','c', ' no family'
        if tlev == 'c':
            z,a,b,d = 'o','c','p', ' no class'
        g,f = taxon.newnames[a],taxon.newnames[b]
        
        
        
        #if g == '?':
        #    return None
        #return outdir+"/"+taxon.tid+tlev+".png"
        
    
    
        taxa = set([t for t,v in self.taxa.items() if v.names[a] == g and v.names[b] == f])
        taxa.add( taxon.tid )
        ntaxa = 0
        tlind = 6 
        while ntaxa < 7:
            if ntaxa:
                tlind -= 1
                tl = dict([(tt,taxon.newnames[tt]) for tt in Taxon.sTaxLevs[:tlind]])
                for t,v in self.taxa.items():
                    if all([v.names[tt] == tl[tt] for tt in Taxon.sTaxLevs[:tlind]]):
                        taxa.add( t )
            ctree = copy.deepcopy( self.nTree )
            ctree.subtree( strategy = 'lca', fn = list(taxa) ) 
            ntaxa = len(ctree.tree.get_terminals())
            tbl = ctree.tree.total_branch_length() 
            ctree.export( tree_tf.name )
            tree_tf.close()

        ann_tf = get_tmp_file()
        title = g+" (tot. branch length "+str(round(tbl,3))+")"
        annotate_nodes( [  ('sub_branches_angle_reduction','0.0'),
                           ('clade_edge_line_width','0.2'), 
                           ('annotation_wings_width',".4"),
                           ('height',"5.5"),
                           ('start_rotation','0'),
                           ('annotation_legend_font_size','12'),
                           ('class_label_font_size','12'),
                           ('title',title), 
                           ],
                        None, ann_tf )
        annotate_nodes( [  ('clade_size','10.0'),
                           ('clade_marker_style',"o"),
                        ],
                        ["*"], ann_tf )
        if tlev == 'g':
            species = dict([(t,self.taxa[t].names['g'][3:]+(" "+self.taxa[t].names['s'][3:] if self.taxa[t].names['s'][3:] else " spp")) for t in taxa])
            species[taxon.tid] = taxon.newnames['g'][3:]+(" "+taxon.newnames['s'][3:] if taxon.newnames['s'][3:] else " spp")
            gspp = taxon.newnames['g'][3:]+" spp"
        else:
            species = dict([(t,self.taxa[t].names[z][3:] if self.taxa[t].names[a][3:] else self.taxa[t].names[a][3:]+d) for t in taxa])
            species[taxon.tid] = taxon.newnames[z][3:] if taxon.newnames[b][3:] else taxon.newnames[a][3:]+d
            gspp = taxon.newnames[z][3:]+d
        
        for i,l in enumerate(sorted(list(set(species.values())),key=lambda x:species.values().count(x),reverse=True)):
            if not l.strip():
                continue
            annotate_nodes( [  ('label',l),
                               ('clade_fill_color','k' if l.endswith(gspp) else col_list[i%len(col_list)]),
                               ('clade_size','30.0' if l.endswith(gspp) else '50.0'),
                               ('annotation_color','k' if l.endswith(gspp) else col_list[i%len(col_list)]),
                              #('annotation_ext_offset',"2.5"),
                            ],
                            [l], ann_tf )
        for k,v in species.items():
            if not v.strip():
                continue
            if k == taxon.tid:
                fr = "from "+taxon.get_changes()
                annotate_nodes( [('class',v),
                             ('annotation_label',taxon.names['t']+": "+fr),
                             ('clade_marker_style',"*"),
                             #('clade_fill_color','b'),
                             ('clade_size','220.0'),
                             ('clade_edge_line_width','1.4'),
                             ('annotation_label_rotation','90'),
                             ],
                             [taxon.tid], ann_tf )
                continue
            annotate_nodes( [('class',v),
                             ('annotation_label',self.taxa[k].names['t']),
                             ('annotation_label_rotation','90'),
                             ],
                             [k], ann_tf )

        ann_tf.close()
        tree_tmp_tf = get_tmp_file(closed = True)
        sb.call( [ c_annotate, "--annot", ann_tf.name, tree_tf.name, tree_tmp_tf ] )
        outtree = outdir+"/"+taxon.tid+tlev+".png"
        sb.call([c_circlader,"--dpi","100",tree_tmp_tf,outtree]) 
        
        #os.remove( ann_tf.name  ) 
        os.remove( tree_tmp_tf  ) 
        os.remove( tree_tf.name  ) 
        if g == '?':
            return None
        return outtree

    def get_lists(self):
        refined =       {0:[],1:[],2:[]}
        corrected =     {0:[],1:[],2:[]}
        removed =       {0:[],1:[],2:[]}
        incompleted =   {0:[],1:[],2:[]}
        newspecies =    []
        
        for tid,taxon in self.taxa.items():
            ref,cor,rem = taxon.refined, taxon.corrected, taxon.torem
            if cor:
                val = max([v['val'] for v in cor.values()])
                """
                if taxon.is_removed():
                    removed[val].append( [taxon.tid,taxon.fta(),taxon.fta(newnames=True),str(val)] )
                else:
                """
                corrected[val].append( [taxon.tid,taxon.fta(),taxon.fta(newnames=True),str(val)] )
            elif ref:
                val = max([v['val'] for v in ref.values()])
                #sup = 0.0
                #if 'ltcs' in ref:
                #    sup += 10.0 + (1 - ref['ltcs']['r']) * ref['ltcs']['n']
                #elif 'nn' in ref:
                #    sup += (1 - ref['nn']['nnd']) * ref['nn']['n']
                    
                refined[val].append( [taxon.tid,taxon.fta(),taxon.fta(newnames=True),str(val)] )
            elif rem:
                val = max([v['val'] for v in rem.values()])
                if taxon.is_removed():
                    removed[val].append( [taxon.tid,taxon.fta(),taxon.fta(newnames=True),str(val)] )
                else:
                    incompleted[val].append( [taxon.tid,taxon.fta(),"",str(val)] )
            if taxon.new_sp:
                newspecies.append( [taxon.tid,taxon.fta( tmax = 'g' )+'.s__'+taxon.new_sp+"."+taxon.names['t']] )

        return (refined,corrected,removed,incompleted,newspecies)

    def write_report( self, report_pdf, tid2img, metadata = None, images = False, typ = "refined" ):

        if metadata:
            data = [l.strip().split('\t') for l in open(metadata)]
            keys = [k.strip() for k in data[0][1:]]

            t2d = {}
            for l in data[1:]:
                t2d[l[0]] = dict([(k,v) for k,v in zip(keys,l[1:])])


        def report( typ ):
            corfn = report_pdf+"_"+typ+'_taxa.tex'

            corf = open(corfn,"w")
            corf.write( "\\documentclass[8pt]{article}\n" )
            corf.write( "\\usepackage[margin=0.75in]{geometry}\n" )
            corf.write( "\\usepackage{graphicx}\n" )
            corf.write( "\\title{List of "+typ+" taxa with high confidence}\n" )
            corf.write( "\\author{\\\\ \\\\Document automatically generated by mToL400 pipeline v. 1.0 \\\\ \\\\ \\\\ "
                        "Nicola Segata (nsegata@hsph.harvard.edu) \\\\ Curtis Huttenhower (chuttenh@hsph.harvard.edu)}\n" )
            corf.write( "\\begin{document}\n" )
            corf.write( "\\maketitle" )
            corf.write( "\\begin{abstract}" )
            if typ == 'refined':
                corf.write( "This document contains suggested taxonomic refinements for public genomes (IMG JGI version 3.5,"
                            "downloaded February 2012) lacking identifiers at one or more taxonomic levels (from phylum to species) "
                            "and for which the automated mToL400 pipeline could confidently assign a name." )
            elif typ == 'corrected':
                corf.write( "This document contains suggested taxonomic corrections for public genomes (IMG JGI version 3.5,"
                            "downloaded February 2012) for which the automated mToL400 could confidently detected and corrected "
                            "taxonomic labels." )

            corf.write( "\n \\\\ \n"
                        "The relevant information associated with each case is reported (e.g. sequencing center and data of submission). "
                        "In addition, for genomes with known genus-level taxonomic labels, "
                        "we have provided a phylogenetic tree rooted into the least common ancestor of all available sequenced strains "
                        "from that particular genus (as defined by IMG taxonomy). "
                        "The star indicatess the "
                        "genome under investigation, and its label corresponds to the new taxonomic assignment. "
                        "The cases are reported in decreasing order of confidence."
                        "\n \\\\ \\\\ \n"
                        "Please do not hesitate to contact us for any comment or suggestion.\n" )
            corf.write( "\\end{abstract}\n" )
            corf.write( "\\clearpage\n" ) 

            if typ == 'refined':
                ordtid = [(tid,taxon) for tid,taxon in self.taxa.items() if taxon.refined and max([v['val'] for v in taxon.refined.values()]) == 2]
                def fsort( x ):
                    ref = x[1].refined
                    sup = 0.0 
                    if 'ltcs' in ref: 
                        sup += 10.0 + (1 - ref['ltcs']['r']) * ref['ltcs']['n']
                    elif 'nn' in ref:
                        sup += (1 - ref['nn']['nnd']) * ref['nn']['n']
                    return sup

                ordtid = sorted( ordtid, key = fsort, reverse = True ) 
            elif typ == 'corrected':
                ordtid = [(tid,taxon) for tid,taxon in self.taxa.items() if taxon.corrected and max([v['val'] for v in taxon.corrected.values()]) == 2]
               
                def fsort( x ):
                    ref = x[1].corrected
                    sup = 0.0 
                    if 'ltcs' in ref: 
                        sup += 10.0 + (1 - ref['ltcs']['r']) * ref['ltcs']['n']
                    elif 'nn' in ref:
                        sup += (1 - ref['nn']['nnd']) * ref['nn']['n']
                    return sup

                ordtid = sorted( ordtid, key = fsort, reverse = True )
            else:
                return

            for tid,taxon in ordtid:
                pr = None
                if typ == 'refined':
                    pr = taxon.refined
                if typ == 'corrected':
                    pr = taxon.corrected
                if typ == 'removed':
                    pr = taxon.torem
                if pr:
                    val = max([v['val'] for v in pr.values()])
                    if int(val) != 2:
                        continue
                    m = t2d[tid] if metadata else None
                    out = ["\\section{Taxon "+tid+(": "+m['Genome Name']+"}" if m else "}")]
                    if m:
                        ids = ", ".join([s+":"+m[s] for s in ['NCBI Taxon ID','RefSeq Project ID','GenBank Project ID']])
                        out.append( "\\paragraph{Other taxon IDs:} "+ids ) 
                    out.append( "\\paragraph{Original taxonomy:} "+taxon.fta(viz=True) )
                    out.append( "\\paragraph{"+typ+" taxonomy:} "+taxon.fta(newnames=True,viz=True) )
                    if m:
                        out.append( "\\paragraph{Add Date:} "+m['Add Date']+" (IMG release "+m["IMG Release"]+")" ) 
                        out.append( "\\paragraph{Sequencing Center:} "+m['Sequencing Center']+"\n" ) 
                        out.append( "\\paragraph{Sequencing Status:} "+m['Status']+"\n" ) 
                    if tid in tid2img and tid2img[tid]:
                        out.append( "\\begin{center}\n" )
                        out.append( "\\includegraphics[width=130mm]{"+tid2img[tid].replace("\_","_")+"}" )
                        out.append( "\\end{center}\n" )
                    out.append( "\\clearpage" ) 

                    corf.write( "\n".join( [o.replace("_","\_") for o in out] ) +"\n" )
    
            corf.write( "\\end{document}\n" )
            corf.close()
            sb.call( ["pdflatex","-output-directory",os.path.dirname(corfn),corfn] )

        report( typ )

    def write_new_taxonomy( self, outf, outdiffs, outboth ):
        outt = {}
        diffs = {}
        both = {}
        for tid,taxon in self.taxa.items():
            both[tid] = [taxon.fta()]
            ref,cor,rem = taxon.refined, taxon.corrected, taxon.torem
            if cor and max([v['val'] for v in cor.values()]) == 2:
                outt[tid] = taxon.fta(newnames=True)
                both[tid].append( outt[tid] )
                diffs[tid] = [taxon.fta(),outt[tid]]
            elif ref and max([v['val'] for v in ref.values()]) == 2:
                outt[tid] = taxon.fta(newnames=True)
                both[tid].append( outt[tid] )
                diffs[tid] = [taxon.fta(),outt[tid]]
            elif rem and max([v['val'] for v in rem.values()]) == 2 and taxon.is_removed():
                outt[tid] = taxon.fta(newnames=True)
                both[tid].append( outt[tid] )
                diffs[tid] = [taxon.fta(),outt[tid]]
            else:
                outt[tid] = taxon.fta()
                both[tid].append( outt[tid] )


        with open( outf, "w" ) as out:
            for k,v in sorted( outt.items(), key=lambda x:x[0]):
                out.write( "\t".join( [str(k),v] ) + "\n" )
        with open( outdiffs, "w" ) as out:
            out.write( "\t".join( ["genome_ID","original","new"] ) + "\n" )
            for k,v in sorted( diffs.items(), key=lambda x:x[0]):
                out.write( "\t".join( [str(k),v[0],v[1]] ) + "\n" )
        with open( outboth, "w" ) as out:
            out.write( "\t".join( ["genome_ID","original","new"] ) + "\n" )
            for k,v in sorted( both.items(), key=lambda x:x[0]):
                out.write( "\t".join( [str(k),v[0],v[1]] ) + "\n" )


    def write_output( self, proj, outname = None, images = False ):
        out_fol = "output/"+proj+"/"
        out_imgs = {}

        if images:
            for f in ['refined','corrected','removed','incompleted']:
                out_imgs[f] = {}
                if not os.path.exists(  out_fol+f+"/" ):
                    os.mkdir( out_fol+f+"/")
                for c in [0,1,2]:
                    out_imgs[f][c] = out_fol+f+"/"+str(c)+"/"
                    if not os.path.exists( out_imgs[f][c] ):
                        os.mkdir( out_imgs[f][c] )


        nspp, narem, ntorem = 0, 0, 0
        nuncl_tot, nuncl = 0, dict([(t,0) for t in Taxon.sTaxLevs[:-1]])
        nred, nltcsout = 0,0
        tid2img = {}
        for tid,taxon in self.taxa.items():
            if tid in self.to_skip:
                continue

            if taxon.refined and images:
                if int(v['val']) == 2: 
                    tid2img[tid] = self.save_circ_img( taxon, 'refined', out_imgs['refined'][v['val']] )
            if taxon.corrected and images:
                if int(v['val']) == 2: 
                    tid2img[tid] = self.save_circ_img( taxon, 'corrected', out_imgs['corrected'][v['val']] )


            if taxon.glev: nspp += 1       
            f = False
            for t in Taxon.sTaxLevs[:-1]:
                if taxon.names[t] == Taxon.unkwn:
                    nuncl[t] += 1
                    f = True
            if f: nuncl_tot += 1
            if taxon.torem and ('ltcs_out' in taxon.torem or 'redundant_clade' in taxon.torem):
                narem += 1
                if 'ltcs_out' in taxon.torem:
                    nltcsout += 1
                if 'redundant_clade' in taxon.torem:
                    nred += 1
            if taxon.torem: ntorem += 1
        
        refined,corrected,removed,incompleted,newspecies = self.get_lists()

        for conf in [0,1,2]:
            exten = ".txt" if outname is None else "_"+outname
            if conf in refined:
                out = "\n".join(["\t".join(r) for r in refined[conf] if r[0] not in self.to_skip and "d__?" not in r[1]])
                if out:
                    with open(out_fol+"refined_conf_"+str(conf)+exten,"w") as outf:
                        outf.write(out+"\n")
            if conf in refined:
                out = "\n".join(["\t".join(r) for r in refined[conf] if r[0] not in self.to_skip and "d__?" in r[1]])
                if out:
                    with open(out_fol+"imputed_conf_"+str(conf)+exten,"w") as outf:
                        outf.write( out +"\n")
            if conf in corrected:
                out = "\n".join(["\t".join(r) for r in corrected[conf] if r[0] not in self.to_skip])
                if out:
                    with open(out_fol+"corrected_conf_"+str(conf)+exten,"w") as outf:
                        outf.write( out +"\n") 
            if conf in removed:
                out = "\n".join(["\t".join(r) for r in removed[conf] if r[0] not in self.to_skip])
                if out:
                    with open(out_fol+"removed_conf_"+str(conf)+exten,"w") as outf:
                        outf.write(out +"\n") 
            if conf in incompleted:
                out = "\n".join(["\t".join(r) for r in incompleted[conf] if r[0] not in self.to_skip])
                if out:
                    with open(out_fol+"incomplete_conf_"+str(conf)+exten,"w") as outf:
                        outf.write(out +"\n") 
        if newspecies:
            with open(out_fol+"subgenera.txt","w") as outf:
                outf.write("\n".join(["\t".join(r) for r in newspecies if r[0] not in self.to_skip]) +"\n")
        return tid2img

    def evaluate( self, proj, true_tax, wrong_tax, pred_name, descr = None ):
        out_fol = "output/"+proj+"/"
        with open( wrong_tax ) as inpf:
            d_errs =  dict([l.strip().split('\t') for l in inpf if l.count("t__synth_err_")])
        with open( true_tax ) as inpf:
            d_cors = dict([l.strip().split('\t') for l in inpf if l.strip().split('\t')[0] in d_errs])
        
        types = ["refined","corrected","removed","incomplete"] 

        ress = []

        for conf in [0,1,2]:
            detected = 0
            out = {}
            nres = {"refined":{'ok':0,'ko':0}, "corrected":{'ok':0,'ko':0},"removed":0,"incomplete":0}

            for typ in types:
                with open(  out_fol+typ+"_conf_"+str(conf)+"_"+pred_name ) as inpf:
                    for tid,ftax,ttax,cconf in (l.strip().split('\t') for l in inpf if l.count("t__synth_err_")):
                        detected += 1
                        ok = ttax.split(".t__")[0] == d_cors[tid].split(".t__")[0] 
                        out[tid] = [str(int(ok)),typ,str(conf),d_cors[tid],ftax,ttax]
                        if typ in ["refined","corrected"]:
                            if ok: nres[typ]['ok'] += 1
                            else: nres[typ]['ko'] += 1
                        else:
                            nres[typ] += 1

            ress.append( [   descr,conf,detected,
                             nres['refined']['ok'],nres['refined']['ko'],
                             nres['corrected']['ok'],nres['corrected']['ko'],
                             nres['removed'],nres['incomplete'] ] )

        with open( out_fol+"evaluate_"+pred_name, "w" ) as outf:
            outf.write( "\n".join( ["\t".join([k]+v) for k,v in out.items()]) +"\n"  )
            outf.write( "\n"+"\t".join(['descr','conf','detected',
                                        'refined_ok','refined_ko',
                                        'corrected_ok','corrected_ko']+types[-2:]) + "\n")
            for rr in ress:
                outf.write( "\t".join([str(s) for s in rr]) + "\n") 

    def infer_tlabels( self, synth_err = False ):
        if self.clades2terms is None:
            self.clades2terms = ppa.clades2terms( self.nTree.tree ) 
        
        self.totbrlen = self.nTree.tree.total_branch_length()
        if self.dists is None:
            self.dists = ppa.dist_matrix( self.nTree.tree )
        if not len(self.var) or self.var is None:
            self.comp_variability()
       
        # cicle over all taxa and focus on the one that are marked to be suspicious
        for tid,taxon in self.taxa.items():
            if tid in self.to_skip:
                continue
            if taxon.torem:
                if synth_err and 'synth_err' not in taxon.names['t']:
                    continue
                # currently 
                conf = max(v['val'] for v in taxon.torem.values())
                if conf < 2:
                    continue
                
                # for taxonomic labels from phylum to species
                for tlev in Taxon.sTaxLevs[1:-1]:
                    # find the largest otherwise fully consistent subtree comprising the taxon of interest
                    lcca,lcca_terms,off_c = self.lcca( taxon, c2t = self.clades2terms, tlev = tlev, parent_steps = 0 )
                    if lcca:
                        terms = [self.taxa[lt.name].fta(tmax=tlev) for lt in lcca_terms]
                        r = self.dist_rank( taxon, lcca_terms )

                        lconf = 0
                        if r < 0.95: lconf = 1
                        if r < 0.9: lconf = 2
                        
                        if lconf < max( [-1] + 
                                        [tv['val'] for tv in taxon.refined.values()] +
                                        [tv['val'] for tv in taxon.corrected.values()] ):
                            continue

                        taxon.copy2new( self.taxa[lcca_terms[0].name], tmax = tlev )
                        
                        if taxon.is_refinement():
                            taxon.refined['ltcs'] = {'val':lconf,'r':r,'n':len(lcca_terms)}
                        elif taxon.is_corrected():
                            taxon.corrected['ltcs'] = {'val':lconf,'r':r,'n':len(lcca_terms)}
                    else:
                        nns = self.onevsall_dist( taxon, tlev = tlev )
                        if nns is None:
                            continue
                        nnd = nns[0][0]

                        if nnd >= self.var[tlev][90]: # self.totbrlen * 1e-2:
                            continue
                        
                        lconf = 0
                        if nnd < self.var[tlev][75]:
                            lconf = 1
                        if nnd < self.var[tlev][50] or nnd < self.totbrlen * 1e-5:
                            lconf = 2
                       
                        if lconf < max( [-1] + 
                                        [tv['val'] for tv in taxon.refined.values()] +
                                        [tv['val'] for tv in taxon.corrected.values()] ):
                            continue
                        taxon.copy2new( self.taxa[nns[0][1]], tmax = tlev )

                        if taxon.is_refinement():
                            taxon.refined['nn'] = {'val':lconf,'nnd':nnd,'n':len(nns)}
                        elif taxon.is_corrected():
                            taxon.corrected['nn'] = {'val':lconf,'nnd':nnd,'n':len(nns)} 

        new_species = False
        if new_species:
            clusters = collections.defaultdict(list)
            for tid,taxon in ((a,b) for a,b in self.taxa.items() if b.torem):
                if 's' in self.taxa[tid].newnames and self.taxa[tid].newnames['s'] != '?':
                    continue
                found = False
                for cid,tids in clusters.items():
                    for t in tids:
                        if t in self.to_skip:
                            continue
                        dist = self.dists[tid][t]
                        if dist > self.var['s'][90]: continue
                        knwn_clos = [v for k,v in self.taxa.items() if k != t and self.dists[t][k] < dist]
                        knwn_clos = [v for v in knwn_clos if v.names['s'] == '?' and ('s' not in v.newnames or v.newnames['s'] == '?')]
                        if len( knwn_clos ) > 0:
                            continue
                        clusters[cid].append( tid )
                        found = True
                        break
                    if found:
                        break
                if not found:
                    clusters[tid].append(tid)

            cln = 0
            for k,v in clusters.items():
                for t in v:
                    self.taxa[t].new_sp = 'ns_{:03d}'.format(cln)
                cln += 1
    
    def remove_tlabels( self, synth_err = False ):
        # precompute all vs all distances and clade to terminal mappings
        if self.dists is None:
            self.dists = ppa.dist_matrix( self.nTree.tree )
        if self.clades2terms is None:
            self.clades2terms = ppa.clades2terms( self.nTree.tree ) 
        
        # Check whether lower level labels are strongly overrefined
        for tid,taxon in self.taxa.items():
            if taxon.names['s'] == Taxon.unkwn:
                continue
            # check at species level only 
            if not synth_err or (synth_err and 'synth_err' in taxon.names['t']): # taxon.names['s'] != '?':

                #if "cholera" in taxon.names['s']:
                #    print taxon.names

                # identify the largest monophyletic subtree around the taxon exluding the taxon itself
                lcca,lcca_terms,off_c = self.lcca( taxon, c2t = self.clades2terms,  parent_steps = 15 )                 
                if not lcca:
                    #if "cholera" in taxon.names['s']:
                    #    print "not lcca"
                    continue
                

                # if such monphyletic tree exists then compute the closest distances statistics
                all_dists = [   min([self.dists[l1.name][l2.name] 
                                    for l1 in lcca_terms 
                                        if l1.name != l2.name and l1.name != taxon.tid]) 
                                for l2 in lcca_terms 
                                    if l2.name != taxon.tid]
                mind = min([self.dists[l.name][taxon.tid] for l in lcca_terms if l.name != taxon.tid])
                
                sdists = sorted(all_dists+[mind])
                
                # r is the "rank" of the distance to the closest taxon wrt the closest distances in the clade
                r = float(sdists.index(mind))/float(len(sdists))
                nlcca = len(lcca_terms)
                
                #if "cholera" in taxon.names['s']:
                #    print all_dists, mind
                #    print sdists
                #    print "r",r
                #    print "nlcca",nlcca

                if r <= 0.6 and nlcca > 7:
                    taxon.torem['redundant_clade'] = {'val':2,'r':r,'nlcca':nlcca}
                elif r <= 0.7 and nlcca > 5:
                    taxon.torem['redundant_clade'] = {'val':1,'r':r,'nlcca':nlcca} 
                elif r <= 0.8 and nlcca > 2:
                    taxon.torem['redundant_clade'] = {'val':0,'r':r,'nlcca':nlcca}
                else:
                    continue
                taxon.remove_all()
        
        ids2clades = dict([(t.name,t.name) for t in self.nTree.tree.get_terminals()])
        tc2t = self.nTree.get_c2t()
      
        
        self.d_ltcs = {}
        self.d_lca = {}
        self.d_ltcs_fra = {}


        # for all tax levels from species to phyla
        for tlev in Taxon.sTaxLevs[1:-1][::-1]:
            c2t = collections.defaultdict(list)
            for tid,taxon in self.taxa.items():
                sTaxon = taxon.fta(tmax = tlev, kwn_only = True) if tlev == 's' else taxon.fta(tmax = tlev )
                c2t[sTaxon].append( tid )
            
            if tlev not in self.d_ltcs:
                self.d_ltcs[tlev] = {}
                self.d_lca[tlev] = {}
                self.d_ltcs_fra[tlev] = {}

            # Check the consistency of the labels in each clade 
            for c,terms in c2t.items(): 
                # If the clade has only two taxa, no way to fix it!
                lterms = len(terms)
                if lterms < 3 or c is None:
                    continue
                # Find the LCA and LTCS for the taxa in the clade
                lca = self.nTree.lca( terms, terminals2clades = ids2clades )
                ltcs = self.nTree.ltcs( terms, tc2t = tc2t, terminals2clades = ids2clades,
                                        lca_precomputed = lca )
                lca_terms = set(self.clades2terms[lca])
                ltcs_terms = set(self.clades2terms[ltcs])
                
                ltcs_fra = float(len(ltcs_terms))/float(lterms)
                

                #if ltcs_fra < 0.7:
                #    continue
                if len(ltcs_terms) < 4:
                    continue
            
                self.d_ltcs[tlev][c], self.d_lca[tlev][c], self.d_ltcs_fra[tlev][c] = ltcs, lca, ltcs_fra
                
                # Find the taxa outside the LTCS  
                out_terms = set([t for t in lca_terms if t.name in terms]) - ltcs_terms
                outs = [c]
                if len(out_terms):
                    # find the radius of the LTCS
                    #rad = max(ltcs.depths().values()) 
                    rad = stats.scoreatpercentile(list(ltcs.depths().values()), 75)
                    
                    rtmp = [(t,self.nTree.tree.distance(ltcs,t)/rad if rad != 0.0 else 1000.0) for t in out_terms]
                    #rtmp = [(t,dists[ltcs.name][t.name]/rad) for t in out_terms]
                    for r,d in rtmp:

                        if d < 1.0:
                            continue

                        dltcs = sorted(
                                [(o,self.dists[o.name][r.name]) for o in ltcs_terms 
                                    if o.name != r.name and o.name not in self.to_skip], #  or o not in lca_terms],
                             key=lambda x:x[1]   )
                        dlca = sorted(
                                [(o,self.dists[o.name][r.name]) for o in lca_terms 
                                    if o.name != r.name and o.name not in self.to_skip],
                                key=lambda x:x[1])

                        if not dltcs or not dlca:
                            continue

                        # if the closest taxon with the same label is not in the
                        # LCA then remove the name
                        cltcslca_r = dltcs[0][1]/dlca[0][1] if dlca[0][1] != 0 else 10.0
                        
                        wl = self.taxa[dlca[0][0].name].fta(tmax = tlev, kwn_only = True ) if tlev == 's' \
                                else self.taxa[dlca[0][0].name].fta(tmax = tlev )
                        rl = self.taxa[r.name].fta(tmax = tlev,kwn_only = True ) if tlev == 's' \
                                else self.taxa[r.name].fta(tmax = tlev )
                        
                        
                        if ( wl != rl or ltcs_fra > 0.9 )  and cltcslca_r > 1.0:
                            taxon = self.taxa[r.name]
                            if ltcs_fra > 0.8 and d > 2.0 and cltcslca_r > 2.0:
                                taxon.torem['ltcs_out'] = {'val':2,'off_cl':c,'ltcs_fra':ltcs_fra,'d':d,'cltcslca_r':cltcslca_r}
                                taxon.remove_all()
                            elif ltcs_fra > 0.7 and d > 1.25 and cltcslca_r > 1.5:
                                if 'ltcs_out' not in taxon.torem or taxon.torem['ltcs_out']['val'] == 0:
                                    taxon.torem['ltcs_out'] = {'val':1,'off_cl':c,'ltcs_fra':ltcs_fra,'d':d,'cltcslca_r':cltcslca_r}
                                    taxon.remove_all()
                            else:
                                taxon.torem['ltcs_out'] = {'val':0,'off_cl':c,'ltcs_fra':ltcs_fra,'d':d,'cltcslca_r':cltcslca_r}
                                taxon.remove_all()
        
