#!/usr/bin/env python

import sys

try:
    import argparse as ap
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

try:
    import pyphlan as ppa
except ImportError:
    sys.stderr.write( "pyphlan.py not found\n" )
    sys.exit(-1)


def read_params( args ):
    p = ap.ArgumentParser(description='Reroot trees')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('outtree', nargs='?', default=None, type=str,
            help=   "the output PhyloXML tree [stdout if not present]")
    p.add_argument('-f', metavar="File containing clade names",
            default=None, type=str )
    p.add_argument( '-n', metavar="N leaves for longest_internal_edge_n", 
                    default=1, type=int)
    p.add_argument( '--name', metavar="Name of the new root (the name must already exists)", 
                    default=None, type=str)

    reroot_st = ['name','lca','ltcs','longest_edge','longest_internal_edge','longest_internal_edge_n']
    #reroot_st = ['lca','ltcs','midpoint','longest_edge','longest_internal_edge','longest_internal_edge_n']
    p.add_argument('-s', choices=reroot_st, default='longest_internal_edge', type=str,
            help=  "select rerooting strategy")

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    tree = ppa.PpaTree( args['intree'] )
    tree.reroot( strategy = args['s'], tf = args['f'], n = args['n']) # , name = args['name'] )
    tree.export( args['outtree'] )

