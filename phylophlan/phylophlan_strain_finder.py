#!/usr/bin/env python


__author__ = ('Francesco Asnicar (f.asnicar@unitn.it), '
              'Claudia Mengoni (claudia.mengoni@studenti.unitn.it)')
__version__ = '3.0.11'
__date__ = '19 May 2021'


import argparse as ap
import os
import sys
import time
import datetime
import pandas as pd
import itertools
from Bio import Phylo


if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn {} requires Python 3, your current Python version is {}.{}.{}"
                    .format(__version__, sys.version_info[0], sys.version_info[1], sys.version_info[2]))

NEWICK = 'newick'
TREE_TYPES = ['newick', 'nexus', 'phyloxml', 'cdao', 'nexml']
PHYLOGENETIC_THR = 0.05
MUT_PERCENTAGE_THR = 0.05
OUTPUT_EXTENSIONS = {';': '.txt',
                     ',': '.csv',
                     '\t': '.tsv'}


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
    p = ap.ArgumentParser(description=("The phylophlan_strain_finder.py script analyzes the phylogeny and the mutation rates table "
                                       "generated from the phylophlan.py script and returns sub-trees representing the same strain, "
                                       "according to both a phylogenetic threshold (computed on the normalized pairwise phylogenetic "
                                       "distances) and a mutation rate threshold (computed on the aligned sequences of the markers used "
                                       "in the phylogenetic analysis)"),
                          formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input', type=str, required=True,
                   help='Specify the file of the phylogenetic tree as generated from phylophlan.py')
    p.add_argument('-m', '--mutation_rates', type=str, required=True,
                   help='Specify the file of the mutation rates as generated from phylophlan.py')
    p.add_argument('--phylo_thr', type=float, default=PHYLOGENETIC_THR,
                   help='Maximum phylogenetic distance threshold for every pair of nodes in the same subtree (inclusive)')
    p.add_argument('--mutrate_thr', type=float, default=MUT_PERCENTAGE_THR,
                   help='Maximum mutation rate ratio for every pair of nodes in the same subtree (inclusive)')
    p.add_argument('--tree_format', choices=TREE_TYPES, default=NEWICK, help='Specify the format of the input tree')
    p.add_argument('-o', '--output', type=str, default=None, help='Specify the output filename, if not specified will be stdout')
    p.add_argument('--overwrite', action='store_true', default=False, help='Overwrite the output file if exists')
    p.add_argument('-s', '--separator', type=str, default='\t', choices=OUTPUT_EXTENSIONS.keys(),
                   help='Specify the separator to use in the output')
    p.add_argument('--citation', action='version',
                   version=('Asnicar, F., Thomas, A.M., Beghini, F. et al. '
                            'Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using PhyloPhlAn 3.0. '
                            'Nat Commun 11, 2500 (2020). '
                            'https://doi.org/10.1038/s41467-020-16366-7'),
                   help="Show citation")
    p.add_argument('--verbose', action='store_true', default=False, help='Write more stuff')
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_strain_finder.py version {} ({})'.format(__version__, __date__),
                   help="Prints the current phylophlan_strain_finder.py version and exit")

    return p.parse_args()


def check_params(args, verbose=False):
    if (not os.path.isfile(args.input)):
        error('input file {} does not exist'.format(args.input), exit=True)

    if (not os.path.isfile(args.mutation_rates)):
        error('mutation_rates file {} does not exist'.format(args.map), exit=True)

    if args.output is None:
        args.output = sys.stdout
    else:
        if(os.path.dirname(args.output)):
            if not os.path.exists(os.path.dirname(args.output)):
                error('output path does not exists: "{}"'.format(args.output), exit=True)

        if not os.path.basename(args.output):
            old = args.output
            out_ext = OUTPUT_EXTENSIONS[args.separator] if args.separator in OUTPUT_EXTENSIONS else '.txt'
            args.output = os.path.join(args.output, 'phylophlan_strain_finder.{}'.format(out_ext))
            info('No output filename specified "{}", writing output to "{}"'.format(old, args.output))

        if (not args.overwrite) and os.path.isfile(args.output):
            old = args.output
            args.output = '{}_{}{}'.format(os.path.splitext(args.output)[0],
                                           datetime.datetime.today().strftime('%Y%m%d%H%M%S'),
                                           os.path.splitext(args.output)[1])
            info('Output file "{}" exists, new output filename is "{}"\n'.format(old, args.output))

    if args.phylo_thr < 0.0:
        error('phylo_thr should be a positive number', exit=True)

    if args.mutrate_thr < 0.0:
        error('mutrate_thr should be a positive number', exit=True)

    if verbose:
        info('Arguments: {}\n'.format(vars(args)))


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)

    if len(node_path) > 1:
        return node_path[-2]
    else:
        return child_clade


def check_thr(p, l, tree, md, p_thr, m_thr, verbose=False):
    if p == l:
        if verbose:
            info('Root reached, return {} as root of the subtree\n'.format(l))

        return l

    if any(((tree.distance(other_children, l_children)) / tree.total_branch_length()) > p_thr
           for other_children in p.get_terminals()
           for l_children in l.get_terminals()):
        if verbose:
            info('Not every leaf under {} respects the phylogenetic threshold,\n'
                 'return {} as root of the subtree\n'.format(p, l))

        return l
    else:
        sons = [s.name for s in p.get_terminals()]
        tup = list(itertools.combinations(sons, 2))
        tup = [(sorted(t)) for t in tup]
        if any(md[(t[0], t[1])] > m_thr for t in tup):
            if verbose:
                info('Not every leaf under "{}" respects the mutation_rates threshold, return "{}"" as root of the subtree\n'.format(p, l))
            return l
        else:
            return check_thr(get_parent(tree, p), p, tree, md, p_thr, m_thr, verbose=verbose)


def phylophlan_strain_finder():
    args = read_params()

    if args.verbose:
        info('phylophlan_strain_finder.py version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args, args.verbose)

    if args.verbose:
        info('Loading phylogenetic tree...\n')

    tree = Phylo.read(args.input, args.tree_format)

    if args.verbose:
        info('Loading mutation rates table...\n')

    mut_rates = pd.read_csv(args.mutation_rates, sep='\t', header=0, index_col=0)
    mut_rates.index = mut_rates.index.map(str)
    mut_rates.columns = mut_rates.columns.map(str)

    mydict = dict([((r, c), float(mut_rates.at[r, c]))
                   for ir, r in enumerate(mut_rates.index)
                   for ic, c in enumerate(mut_rates.columns) if ir < ic])
    mydict.update(dict([((c, r), float(mut_rates.at[r, c]))
                        for ir, r in enumerate(mut_rates.index)
                        for ic, c in enumerate(mut_rates.columns) if (ir < ic) and (c, r) not in mydict]))
    tested_valid = []

    for l in tree.get_terminals():
        if any(l in x.get_terminals() for x in tested_valid):
            continue

        r = (check_thr(get_parent(tree, l), l, tree, mydict, args.phylo_thr, args.mutrate_thr, args.verbose))
        tested_valid.append(r)

    if args.verbose:
        info('Creating output...\n')

    f = args.output

    if isinstance(args.output, str):
        f = open(args.output, 'w')

    print('#phylogenetic_threshold{}{}'.format(args.separator, args.phylo_thr), file=f)
    print('#mutation_rate_threshold{}{}'.format(args.separator, args.mutrate_thr), file=f)
    print('#total_branch_length{}{}'.format(args.separator, tree.total_branch_length()), file=f)
    print(args.separator.join(['#subtree', 'min_dist', 'mean_dist', 'max_dist', 'min_mut',
                               'mean_mut', 'max_mut', 'distances', 'mutation_rates']), file=f)

    for test in tested_valid:
        subtree = Phylo.Newick.Tree(test).format(args.tree_format).strip()
        sons = [s for s in test.get_terminals()]

        if len(sons) == 1:
            continue

        tup = list(itertools.combinations(sons, 2))
        m_rate = []
        distances = []
        m_min = d_min = tree.total_branch_length()
        m_max = m_mean = d_max = d_mean = count = 0

        for t in tup:
            m = mydict[(t[0].name, t[1].name)]
            d = tree.distance(t[0], t[1])
            distances.append(t[0].name + ',' + t[1].name + ':' + str(d))
            m_rate.append(t[0].name + ',' + t[1].name + ':' + str(m))
            m_min = min(m_min, m)
            m_max = max(m_max, m)
            m_mean = m_mean + m
            d_min = min(d_min, d)
            d_max = max(d_max, d)
            d_mean = d_mean + d
            count = count + 1

        d_mean = d_mean / count
        m_mean = m_mean / count

        print(args.separator.join([subtree, str(d_min), str(d_mean), str(d_max), str(m_min), str(m_mean), str(m_max),
                                   '|'.join(distances), '|'.join(m_rate)]), file=f)

    if isinstance(args.output, str):
        if not f.closed:
            f.close()


if __name__ == '__main__':
    t0 = time.time()
    phylophlan_strain_finder()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
