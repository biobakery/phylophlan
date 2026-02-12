#!/usr/bin/env python


__description__ = "A script to calculate the pairwise patristic (leaf-to-leaft) distance matrix from a phylogenetic" \
                  " tree. It runs efficiently in O^2 time and space and can be parallelized."
__author__ = ", ".join((
    'Michal Puncochar',
))
__version__ = '3.2.0'
__date__ = '8 October 2024'


import argparse as ap
import pathlib
import multiprocessing as mp

import dendropy
from dendropy.datamodel.treemodel import Node
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd

from tqdm.auto import tqdm

from .utils import info, ArgumentType, openr


DTYPE = np.float64

condensed_dist_array: mp.RawArray | None = None
all_leaves: list[Node] = []


def condensed_index(i: int, j: int, n: int):
    if i > j:
        i, j = j, i  # upper triangle only
    c = (n * i) - (i * (i + 1) // 2) - 1 - i + j

    return c


def condensed_index_vect(i, j, n):
    i = np.array(i)
    j = np.array(j)
    m = i > j
    i[m], j[m] = j[m], i[m]
    c = (n * i) - (i * (i + 1) // 2) - 1 - i + j
    return c


class DistMat:
    def __init__(self, cmat, idx, diagonal_value=0):
        self.cmat = cmat
        self.index = pd.Index(idx)
        self.n = len(idx)
        self.idx_to_i = pd.Series(index=idx, data=np.arange(self.n))
        self.diagonal_value = diagonal_value

    @classmethod
    def from_npy_file(cls, f, in_memory=False):
        f = str(f)
        if in_memory:
            mmap_mode = None
        else:
            mmap_mode = 'r'
        cmat = np.load(f, allow_pickle=False, mmap_mode=mmap_mode)
        idx = np.loadtxt(f + '.index', dtype=str)
        return cls(cmat, idx)

    def condensed_index(self, i, j):
        if i >= self.n or j >= self.n or i < 0 or j < 0:
            raise Exception("i and j have to be integers between 0 and n-1")

        if i == j:
            return -1

        return condensed_index(i, j, self.n)

    def condensed_index_vect(self, i, j):
        i = np.array(i)
        j = np.array(j)
        if (i >= self.n).any() or (j >= self.n).any() or (i < 0).any() or (j < 0).any():
            raise Exception("i and j have to be integers between 0 and n-1")

        r = condensed_index_vect(i, j, self.n)
        r[i == j] = -1
        return r

    def get_i(self, i, j):
        ci = self.condensed_index(i, j)
        if ci == -1:
            return self.diagonal_value
        return self.cmat[ci]

    def get(self, idx_i, idx_j):
        i, j = self.idx_to_i[idx_i], self.idx_to_i[idx_j]
        return self.get_i(i, j)

    def get_slice_i(self, ii, jj):
        iv = np.repeat(ii, len(jj))
        jv = np.tile(jj, len(ii))
        ci = self.condensed_index_vect(iv, jv)

        x = self.cmat[ci]
        x[ci == -1] = self.diagonal_value
        return x.reshape((len(ii), len(jj)))

    def get_slice(self, idx_i, idx_j):
        ii = self.idx_to_i.loc[idx_i]
        jj = self.idx_to_i.loc[idx_j]
        a = self.get_slice_i(ii, jj)
        df = pd.DataFrame(a, index=idx_i, columns=idx_j)
        return df

    def get_subset_i(self, ii):
        iv = np.repeat(ii, len(ii))
        jv = np.tile(ii, len(ii))
        m = iv < jv
        iv, jv = iv[m], jv[m]

        ci = self.condensed_index_vect(iv, jv)
        ci.sort()

        return self.cmat[ci]

    def get_pairs_i(self, ii, jj):
        ci = self.condensed_index_vect(ii, jj)
        x = self.cmat[ci]
        x[ci == -1] = self.diagonal_value
        return x

    def get_pairs(self, idx_i, idx_j):
        ii = self.idx_to_i.loc[idx_i]
        jj = self.idx_to_i.loc[idx_j]
        return self.get_pairs_i(ii, jj)

    def get_subset(self, idx):
        idx_to_i_subset = self.idx_to_i.loc[idx]
        cmat_subset = self.get_subset_i(idx_to_i_subset.values)
        index_subset = idx_to_i_subset.sort_values().index

        return DistMat(cmat_subset, index_subset, self.diagonal_value)

    def to_pandas(self):
        a = spd.squareform(self.cmat)
        if self.diagonal_value != 0:
            np.fill_diagonal(a, self.diagonal_value)
        return pd.DataFrame(a, index=self.index, columns=self.index)


class DistCalculator:
    def __init__(self, n_leaves, label_to_id, dtype):
        self.label_to_id: dict[str, int] = label_to_id
        self.n_leaves: int = n_leaves
        self.dtype = dtype

    def dfs(self, node, path_len, orig_leaf_id, condensed_dist_matrix, skip_child=None):
        """

        :param Node node:
        :param float path_len:
        :param int orig_leaf_id:
        :param np.ndarray condensed_dist_matrix:
        :param Node | None skip_child:
        :return:
        """
        if node.is_leaf():
            leaf_id = self.label_to_id[node.taxon.label]
            idx = condensed_index(orig_leaf_id, leaf_id, self.n_leaves)
            condensed_dist_matrix[idx] = path_len
            return 1
        else:
            ret = 0
            for child in node.child_node_iter():
                if skip_child is not None and child == skip_child:
                    # continue  # can be break here?
                    break
                child_length = child.edge_length
                if child_length is None:
                    child_length = 0
                ret += self.dfs(child, path_len + child_length, orig_leaf_id, condensed_dist_matrix)
            return ret

    def calc_for_leaf(self, leaf_i):
        """

        :param int leaf_i:
        :return:
        """

        global condensed_dist_array, all_leaves

        leaf_1 = all_leaves[leaf_i]
        leaf_1_id = self.label_to_id[leaf_1.taxon.label]
        condensed_dist_matrix = np.frombuffer(condensed_dist_array, self.dtype)

        x = leaf_1
        d = 0
        ret = 0
        while x.parent_node is not None:
            if x.edge_length is not None:
                d += x.edge_length

            ret += self.dfs(x.parent_node, d, leaf_1_id, condensed_dist_matrix, x)
            x = x.parent_node
        return ret


class ArgumentTypes:
    input_tree: pathlib.Path
    output_file: pathlib.Path
    output_format: str
    nproc: int


def read_params():
    p = ap.ArgumentParser(description=__description__, formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input_tree', type=ArgumentType.existing_file, required=True,
                   help="Path to the tree file in Newick format")
    p.add_argument('-o', '--output_file', type=ArgumentType.creatable_file, required=True,
                   help="Output file to write the distance matrix to")
    p.add_argument('-f', '--output_format', type=str, choices=['tsv', 'parquet', 'npy', 'tri'], default=None,
                   help="Output format, if not specified will be determined from the output_file extension.\n"
                        "  - tsv: tab-separated values. A full matrix will be generated, gzip or bzip2 compression"
                        " is supported and is inferred from the output_file extension\n"
                        " - parquet: apache parquet file. A full matrix will be generated, export is handled by pandas."
                        "\n"
                        " - npy: The condensed distance matrix will be written as npy file and the index will be "
                        "written to the [output_file].index file\n"
                        " - tri: The condensed matrix will be written as a lower triangular matrix file, gzip or bzip2 "
                        "compression is supported")
    p.add_argument('--nproc', type=ArgumentType.positive_int, default=1,
                   help="Number of processes")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_patristic_distances.py version {} ({})'.format(__version__, __date__),
                   help="Print the version and exit")

    return p


def check_params(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=ArgumentTypes)

    if args.output_format is None:
        output_name = args.output_file.name
        if output_name.endswith('.parquet'):
            args.output_format = 'parquet'
        elif output_name.endswith('.npy'):
            args.output_format = 'npy'
        elif output_name.endswith('.tri') or output_name.endswith('.tri.gz') or output_name.endswith('.tri.bz2'):
            args.output_format = 'tri'
        else:
            args.output_format = 'tsv'
        info(f'Inferred output format as {args.output_format}')

    return args


def update_bar(q, total):
    pbar = tqdm(total=total)

    while True:
        x = q.get()
        if x == 2:
            pbar.close()
            return
        pbar.update(x)


def main():
    argp = read_params()
    args = check_params(argp)

    info('Loading the tree')
    dtree: dendropy.Tree = dendropy.Tree.get(path=args.input_tree, schema='newick', preserve_underscores=True)

    info('Calculating the distance matrix')

    global all_leaves, condensed_dist_array
    all_leaves: list[Node] = dtree.leaf_nodes()


    n_leaves = len(all_leaves)
    mat_size = n_leaves * (n_leaves - 1) // 2

    label_to_id = dict(zip((x.taxon.label for x in all_leaves), range(len(all_leaves))))

    ctype = np.ctypeslib.as_ctypes_type(DTYPE)
    condensed_dist_array = mp.RawArray(ctype, mat_size)

    dc = DistCalculator(n_leaves, label_to_id, DTYPE)

    with mp.Pool(args.nproc) as pool, tqdm(total=mat_size) as pb:
        for r in pool.imap_unordered(dc.calc_for_leaf, range(len(all_leaves))):
            pb.update(r)

    condensed_dist_matrix = np.frombuffer(condensed_dist_array, DTYPE)
    assert len(condensed_dist_matrix) == mat_size

    s_label_to_id = pd.Series(label_to_id)
    index = s_label_to_id.index


    if args.output_format in ['tsv', 'parquet']:
        info('Converting to square matrix')
        dm = DistMat(condensed_dist_matrix, index)
        df = dm.to_pandas()

        info('Saving the distance matrix')
        if args.output_format == 'tsv':
            df.to_csv(args.output_file, sep='\t')
        elif args.output_format == 'parquet':
            df.to_parquet(args.output_file)
    else:
        info('Saving the distance matrix')
        if args.output_format == 'npy':
            np.save(str(args.output_file), condensed_dist_matrix, allow_pickle=False)
            np.savetxt(str(args.output_file) + '.index', index, fmt='%s')
        elif args.output_format == 'tri':
            with openr(args.output_file, 'wt') as f:
                n = len(index)
                f.write(f'{n}\n')
                for i, idx in enumerate(index):
                    f.write(str(idx))
                    for j in range(i):
                        d = condensed_dist_matrix[condensed_index(i, j, n)]
                        f.write(f'\t{d}')
                    f.write('\n')
        else:
            assert False

    info('Done.')


if __name__ == '__main__':
    main()
