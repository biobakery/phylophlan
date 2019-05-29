 #!/usr/bin/env python3

__author__ = ('Francesco Asnicar (f.asnicar@unitn.it)',
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it)')
__version__ = '0.01'
__date__ = '29 May 2019'


import argparse as ap
import os
import sys
import time
import datetime
import pandas as pd
import itertools
from Bio import Phylo
import tempfile

NEWICK = 'newick'
TREE_TYPES = ['newick','nexus','phyloxml','cdao','nexml']
PHYLOGENETIC_THR = 0.05
MUT_PERCENTAGE_THR = 0.05 

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
    p = ap.ArgumentParser(description="",
                          epilog=(""),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str, required=True,
                   help='Specify the file of the tree generated from phylophlan.py')
    p.add_argument('-m', '--mutation_rates', type=str, required=True,
                   help='Specify the file of the mutation rates generated from phylophlan.py')
    p.add_argument('--p_threshold', type=float, default=PHYLOGENETIC_THR, 
                   help='Specify the phylogenetic threshold you want to test on the tree.'
                   'The phylogenetic distance between any node from the same subtree will be less than this threshold' )
    p.add_argument('--m_threshold', type=float, default=MUT_PERCENTAGE_THR,
                   help='Specify the percentage of mutation rate you want to test on th tree.'
                   'The mutation rate between any node from the same subtree will be less than this threshold')
    p.add_argument('--tree_format', choices=TREE_TYPES, default=NEWICK,
                   help='Specify the format of the input tree.')
    p.add_argument('-o','--output', type=str, default=None,
                   help='Specify the name of the output file, if not specified default is stdout')
    p.add_argument('--overwrite', action='store_true', default=False,
                   help='If specified, will overwrite the output file if it exists' )
    p.add_argument('-s','--separator', type=str, default='\t',
                   help='Specify the separator you want in the output')
    p.add_argument('--verbose', action='store_true', default=False,
                   help='Write more stuff' )
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_strain_finder.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan_strain_finder.py version and exit")
  

    return p.parse_args()


def check_params(args, verbose=False):
    if verbose:
        info('Checking for parameters...\n')
            
    if (not os.path.isfile(args.input)):
        error('input file {} does not exist'.format(args.input), exit=True)
    
    if (not os.path.isfile(args.mutation_rates)):
       error('mutation_rates file {} does not exist'.format(args.map), exit=True)
      
    if (not args.overwrite) and args.output and os.path.isfile(args.output) :
        args.output = args.output+str(datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        
    if not args.output:
        args.output = sys.stdout

    if args.p_threshold < 0:
        error('p_threshold should be a positive number')
    
    if args.m_threshold < 0:
        error('m_threshold should be a positive number')

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    if len(node_path)>1:
        return node_path[-2]
    else:
        return child_clade

def check_thr(p, l, tree, md, p_thr, m_thr, verbose=False):
    if p==l:
        if verbose:
            info('Root reached,\n'
                'return {} as root of the subtree\n'.format(l))
        return l

    if any((tree.distance(other_children,l_children))/tree.total_branch_length()>=p_thr for other_children in p.get_terminals() for l_children in l.get_terminals()):
        if verbose:
            info('Not every leaf under {} respects the phylogenetic threshold,\n'
                 'return {} as root of the subtree\n'.format(p,l))
        return l 

    else:
        sons = [s.name for s in p.get_terminals()]
        tup = list(itertools.combinations(sons,2))
        tup = [(sorted(t)) for t in tup]
        if any(float(md[(t[0],t[1])]) > m_thr for t in tup):
            if verbose:
                info('Not every leaf under {} respects the mutation_rates threshold,\n'
                'return {} as root of the subtree\n'.format(p,l))
            return l
        else:            
            return check_thr(get_parent(tree,p), p, tree, md, p_thr, m_thr, verbose)


def phylophlan_strain_finder():
    args = read_params()
    check_params(args, args.verbose)
    tree = Phylo.read(args.input, args.tree_format)

    mut_rates = pd.read_csv(args.mutation_rates, sep='\t', dtype={'ids':str})
    mut_rates.set_index('ids', inplace=True)

    if args.verbose:
        info('Reading mutation_rates table...\n')
    
    mydict = dict([((i,c),mut_rates.at[i,c]) for i in mut_rates.index for c in mut_rates.columns if i<c])

    tested_valid = []

    for l in tree.get_terminals():
        if any(l in x.get_terminals() for x in tested_valid):
            continue
        r = (check_thr(get_parent(tree, l), l, tree, mydict, args.p_threshold, args.m_threshold, args.verbose)) 
        tested_valid.append(r)
    

    if args.verbose:
        info('Creating output...\n')

    f = args.output

    if isinstance(args.output, str):
        f = open(args.output, 'w')

    print('#phylogenetic_threshold{}{}'.format(args.separator, args.p_threshold), file=f)
    print('#mutation_rate_threshold{}{}' .format(args.separator, args.m_threshold), file=f)
    print('#total_branch_length{}{}'.format(args.separator,tree.total_branch_length()), file=f)
    print(args.separator.join(['#subtree', 'min_dist', 'mean_dist', 'max_dist', 'min_mut', 
                               'mean_mut', 'max_mut', 'distances', 'mutation_rates']), file=f)
    

    for test in tested_valid: 
        subtree = Phylo.Newick.Tree(test).format(NEWICK).strip()

        sons = [s for s in test.get_terminals()]

        if len(sons) == 1:
            continue                        

        tup = list(itertools.combinations(sons,2))
        m_rate = []
        distances= []
        m_min = d_min = tree.total_branch_length()
        m_max = m_mean = d_max = d_mean = count =  0
                
        for t in tup:                    
            if t[0].name < t[1].name:
                m = float(mydict[(t[0].name,t[1].name)]) 
            else:
                m = float(mydict[(t[1].name,t[0].name)])
            d = tree.distance(t[0],t[1])
            distances.append(t[0].name+','+t[1].name+':'+str(d))
            m_rate.append(t[0].name+','+t[1].name+':'+ str(m))
            m_min = min(m_min,m)
            m_max = max(m_max, m)
            m_mean = m_mean + m
            d_min = min(d_min, d)
            d_max = max(d_max, d)
            d_mean = d_mean + d
            count = count +1
        d_mean = float(d_mean/count)
        m_mean = float(m_mean/count)

        print((args.separator).join([subtree, str(d_min), str(d_mean), str(d_max),  
              str(m_min), str(m_mean), str(m_max),'|'.join(distances), '|'.join(m_rate)]), file=f)

    if isinstance(args.output, str):
        if not f.closed:
            f.close()


if __name__ == '__main__':    
    t0 = time.time()
    phylophlan_strain_finder()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)