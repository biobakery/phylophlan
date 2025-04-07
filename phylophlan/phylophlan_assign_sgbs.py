#!/usr/bin/env python


# TODO: update description once the script behavior is finalized
__description__ = "The phylophlan_assign_sgbs.py script assigns the SGB and taxonomy to a given set of input genomes." \
                  " Outputs can be of three types: (1) for each input genomes  returns the list of the closest " \
                  "-n/--how_many SGBs sorted by average Mash  distance; (2) for each input genomes returns the  " \
                  "closest SGB, GGB, FGB, and  reference genomes; (3) returns a all vs. all matrix with all the " \
                  "pairwise mash distances"
__author__ = ", ".join((
    'Francesco Asnicar (f.asnicar@unitn.it)',
    'Michal Puncochar',
    'Katarina Mladenovic (katarina.mladenovic@unitn.it)',
    'Francesco Beghini (francesco.beghini@unitn.it)',
    'Fabio Cumbo (fabio.cumbo@unitn.it)',
    'Mattia Bolzan (mattia.bolzan@unitn.it)',
    'Paolo Manghi (paolo.manghi@unitn.it)',
    'Nicola Segata (nicola.segata@unitn.it)',
))
__version__ = '3.2.0'
__date__ = '8 October 2024'


import argparse as ap
import bz2
import gzip
import os
import pathlib
import shlex
import shutil
import subprocess as sp
import sys
import tarfile
import urllib.parse
from collections import Counter

import numpy as np
import pandas as pd

from .utils import info, ArgumentType, mash_sketch, skani_sketch, load_sgb_txt, mash_paste, mash_dist_block, \
    run_parallel, load_pwd_pandas, load_pandas_series, skani_paste, skani_dist_block, tqdm, load_skani_as_pwd, \
    skani_triangle_big, load_big_triangle_skani, cluster_skani_pwd, error, fix_mash_id, path_get_lock, \
    run_command_with_lock, mash_sketch_aa, load_mash_dist_block, download, exists_with_lock, \
    get_threads_per_run, warning, decompress_file, fix_skani_id

if sys.version_info[0] < 3:
    raise Exception("PhyloPhlAn {} requires Python 3, your current Python version is {}.{}.{}"
                    .format(__version__, sys.version_info[0], sys.version_info[1], sys.version_info[2]))


DEPENDENCY_CHECK_CMDS = {
    'mash': 'mash',
    'skani': 'skani -V',
    'prodigal': 'prodigal -v',
}

DB_LIST_URL = "http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan_SGB_databases_v32.txt"
DEFAULT_DATABASE_FOLDER = 'phylophlan_databases/'

PREFILTER_DIST_THRESHOLD = 0.2
SGB_ASSIGNMENT_ANI = 95
GGB_ASSIGNMENT_DIST = 0.15
FGB_ASSIGNMENT_DIST = 0.3
SGB_CLOSEST_PROFILABLE_MIN_ANI = 85
WORK_DIR_NAME = 'phylophlan_assign_sgbs_tmp'
NOVEL_SGB_PREFIX = 'SGBX'
MAX_OVERLAP_WARNING = 30

SMALL_SGB_THRESHOLD = 100 * 100
CHUNK_SIZE_IN = 1000
CHUNK_SIZE_PER_SGB_IN = 10000
CHUNK_SIZE_PER_SGB_DB = 10000
CHUNK_SIZE_AA = 10000


class ArgumentTypes:
    input_directory: pathlib.Path
    input_extension: str
    output_directory: pathlib.Path
    database_folder: pathlib.Path
    database: str
    list_databases: bool
    nproc_cpu: int
    nproc_io: int
    clean: bool


def read_params():
    p = ap.ArgumentParser(description=__description__, formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input_directory', type=ArgumentType.existing_dir, default=None,
                   help="Input folder containing genomes to be indexed")
    p.add_argument('-e', '--input_extension', type=str, default=None,
                   help="Specify the extension of the input file(s) specified via -i/--input. If not specified will "
                        "use all files in the directory")
    p.add_argument('-o', '--output_directory', type=ArgumentType.creatable_dir, default=None,
                   help="Output directory")
    p.add_argument('--database_folder', type=ArgumentType.existing_dir, default=DEFAULT_DATABASE_FOLDER,
                   help="Path to the folder that contains the database file")
    p.add_argument('-d', '--database', type=str, default=None,
                   help="Database name, available options can be listed using the --database_list parameter")
    p.add_argument('--list_databases', action='store_true', default=False,
                   help="List of all the available databases that can be specified with the -d/--database option")
    p.add_argument('--nproc_cpu', type=ArgumentType.positive_int, default=1,
                   help="The number of threads to use for CPU intensive jobs")
    p.add_argument('--nproc_io', type=ArgumentType.positive_int, default=1,
                   help="The number of CPUs to use for I/O intensive jobs")
    p.add_argument('--clean', action='store_true', default=False,
                   help="Whether to clean the temporary directory (don't reuse intermediate results)")
    p.add_argument('--citation', action='version',
                   version=('Asnicar, F., Thomas, A.M., Beghini, F. et al. '
                            'Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using'
                            ' PhyloPhlAn 3.0. '
                            'Nat Commun 11, 2500 (2020). '
                            'https://doi.org/10.1038/s41467-020-16366-7'),
                   help="Print citation and exit")
    p.add_argument('-v', '--version', action='version',
                   version='phylophlan_assign_sgbs.py version {} ({})'.format(__version__, __date__),
                   help="Print the version and exit")

    return p


def check_params(argp: ap.ArgumentParser):
    args = argp.parse_args(namespace=ArgumentTypes)

    if not args.list_databases and (args.input_directory is None or args.output_directory is None):
        argp.error('input_directory and output_directory required')

    if args.output_directory is not None:
        args.output_directory.mkdir(exist_ok=True)

    return args


def check_dependencies():
    for program, cmd in DEPENDENCY_CHECK_CMDS.items():
        r = sp.run(shlex.split(cmd), capture_output=True)
        if r.returncode != 0:
            error(f"{program} is not installed or not present in the system path", do_exit=True)


def process_input(input_directory, input_extension, work_dir, nproc_io):
    genome_files = list(input_directory.iterdir())
    if len(genome_files) == 0:
        error('Empty input directory', do_exit=True)

    if input_extension is None:
        assert len(genome_files) > 0
        one_genome_name = genome_files[0].name
        if one_genome_name.endswith('.gz'):
            extension = '.gz'
        elif one_genome_name.endswith('.bz2'):
            extension = '.bz2'
        else:
            extension = ''
        if len(extension) > 0:
            one_genome_name = one_genome_name[:-len(extension)]
        input_extension = '.' + one_genome_name.rsplit('.', maxsplit=1)[-1] + extension
        info(f'Using inferred genome extension: {input_extension}')
    # elif not input_extension.startswith('.'):
    #     input_extension = '.' + input_extension

    genome_names = [x.name[:-len(input_extension)] for x in genome_files if x.name.endswith(input_extension)]
    files_ignored = [x for x in genome_files if x.is_file() and not x.name.endswith(input_extension)]
    info(f'Found {len(genome_names)} genomes in the input directory')
    if len(files_ignored) > 0:
        warning(f'  Ignoring {len(files_ignored)} files in the input directory '
                f'not having {input_extension} extension')

    if input_extension.endswith('.gz') or input_extension.endswith('.bz2'):
        info('Detected compressed files in the input, will decompress temporarily')
        input_decompressed_dir = work_dir / 'input_decompressed'
        input_decompressed_dir.mkdir(exist_ok=True)

        if input_extension.endswith('.gz'):
            open_f = gzip.open
            new_ext = input_extension[:-len('.gz')]
        elif input_extension.endswith('.bz2'):
            open_f = bz2.open
            new_ext = input_extension[:-len('.bz2')]
        else:
            assert False, "extension must be .gz or .bz2"


        f_args = [(input_directory / f'{g}{input_extension}',
                   input_decompressed_dir / f'{g}{new_ext}',
                   open_f) for g in genome_names]


        run_parallel(decompress_file, f_args, nproc_io, star=True, processes=True)

        input_directory = input_decompressed_dir
        input_extension = new_ext
    else:
        input_directory = input_directory
        input_extension = input_extension

    return genome_names, input_directory, input_extension


def prefiltering(work_dir, nproc_io, nproc_cpu, database_path, centroids, centroid_to_sgb, input_mash_sketches,
                 genome_names):
    genome_centroid_pairs_file = work_dir / 'genome_centroid_pairs.tsv'
    dists_dir_centroids = work_dir / 'mash_dists_input_centroids'
    if genome_centroid_pairs_file.exists():
        info('  Loading existing genomes-to-closest-centroids file')
        df_genome_centroid_pairs = pd.read_csv(genome_centroid_pairs_file, sep='\t')
        dist_files_centroids = sorted(dists_dir_centroids.iterdir())
    else:
        centroid_pastes_dir = database_path / 'centroid_mash_pastes'
        input_pastes_dir = work_dir / 'mash_pastes_input'


        mash_pastes_centroids = list(centroid_pastes_dir.iterdir())
        info('  Pasting input MASH sketches')
        input_pastes_dir.mkdir(exist_ok=True)
        mash_pastes_in = mash_paste(input_mash_sketches, input_pastes_dir, nproc_io, chunk_size=CHUNK_SIZE_IN)

        info('  Calculating distances from input genomes to SGB centroids')
        dists_dir_centroids.mkdir(exist_ok=True)
        dist_files_centroids = mash_dist_block(mash_pastes_in, mash_pastes_centroids, dists_dir_centroids, nproc_cpu)

        info('  Processing distance files')
        genome_centroid_pairs = []
        processed_genomes = []
        for df in run_parallel(lambda dist_f: load_pwd_pandas(dist_f), dist_files_centroids, nproc=nproc_io,
                               ordered=False, return_gen=True):

            df.index = df.index.map(fix_mash_id)  # centroid M-ids

            assert set(df.index) == set(centroids)
            processed_genomes.extend(df.columns.to_list())
            a, b = np.where(df.values <= PREFILTER_DIST_THRESHOLD)
            for c, g, d in zip(df.index[a], df.columns[b], df.values[a, b]):
                genome_centroid_pairs.append((g, c, d))

        assert set(processed_genomes) == set(genome_names)

        df_genome_centroid_pairs = pd.DataFrame.from_records(genome_centroid_pairs,
                                                             columns=['genome', 'centroid', 'distance'])

        df_genome_centroid_pairs['sgb_id'] = df_genome_centroid_pairs['centroid'].map(centroid_to_sgb)

        df_genome_centroid_pairs.to_csv(genome_centroid_pairs_file, sep='\t', index=False)


    return df_genome_centroid_pairs, dist_files_centroids


def dist_against_sgbs(nproc_io, nproc_cpu, work_dir, database_path, sgb_to_genomes_in, sgb_to_genomes_db,
                      input_skani_sketches_dir, sgb_to_mp, input_extension):
    input_per_sgb_pastes_dir = work_dir / 'skani_pastes_input_per_sgb'
    db_per_sgb_pastes_dir = work_dir / 'skani_pastes_db_per_sgb'
    dists_dir_db = work_dir / 'skani_dists_input_db_per_sgb'
    long_results_file = work_dir / 'ani_results_long.tsv'

    if exists_with_lock(long_results_file):
        info('Loading DB ANI long results')
        df_long_results = pd.read_csv(long_results_file, sep='\t')
        return df_long_results


    info('Calculating average distances from input genomes to prefiltered SGBs')

    info('  Creating per-SGB skani pastes of input sketches')
    input_per_sgb_pastes_dir.mkdir(exist_ok=True)
    sgb_to_pastes_in = {}
    for sgb_id, sgb_genomes in tqdm(sgb_to_genomes_in.items()):
        sgb_genome_sketches = [input_skani_sketches_dir / f'{g}.sketch' for g in sgb_genomes]
        sgb_to_pastes_in[sgb_id] = skani_paste(sgb_genome_sketches, input_per_sgb_pastes_dir, CHUNK_SIZE_PER_SGB_IN,
                                               file_prefix=f'paste_sgb{sgb_id}')

    info('  Creating per-SGB skani pastes of DB sketches')
    db_per_sgb_pastes_dir.mkdir(exist_ok=True)
    sgb_to_pastes_db = {}
    for sgb_id in tqdm(sgb_to_genomes_in.keys()):
        sgb_db_sketches = sorted((database_path / 'per_sgb_skani_pastes' / f'{sgb_id}').iterdir())
        sgb_to_pastes_db[sgb_id] = skani_paste(sgb_db_sketches, db_per_sgb_pastes_dir, CHUNK_SIZE_PER_SGB_DB,
                                               file_prefix=f'paste_sgb{sgb_id}')

    info('  Running skani distance calculation')
    small_sgbs = []
    big_sgbs = []
    for sgb_id in sgb_to_genomes_in.keys():
        genomes_in = sgb_to_genomes_in[sgb_id]
        genomes_db = sgb_to_genomes_db[sgb_id]
        size = len(genomes_in) * len(genomes_db)
        if size < SMALL_SGB_THRESHOLD:
            small_sgbs.append(sgb_id)
        else:
            big_sgbs.append(sgb_id)

    info(f'  Will run {len(small_sgbs)} small and {len(big_sgbs)} big SGBs')
    dists_dir_db.mkdir(exist_ok=True)

    def f(sgb_id_):
        r = skani_dist_block(sgb_to_pastes_in[sgb_id_], sgb_to_pastes_db[sgb_id_], dists_dir_db,
                             nproc=nproc_per_run, progress_bar=False)
        return sgb_id_, r

    info('  Running small SGBs')
    nproc_per_run = get_threads_per_run(nproc_cpu, nproc_io)
    results_small = run_parallel(f, small_sgbs, nproc=nproc_io, ordered=False)

    if len(big_sgbs) > 0:
        info('  Running big SGBs')
        nproc_per_run = nproc_cpu
        results_big = run_parallel(f, big_sgbs, nproc=1, ordered=False)
    else:
        results_big = []

    sgb_to_dist_files_db = dict(results_small + results_big)

    info('  Processing distance files')
    dfs_partial = []
    for sgb_id, dist_files_matrix in sgb_to_dist_files_db.items():
        for dist_files_row in dist_files_matrix:
            # compose input genome batch (dist files row) from all the DB genome batches
            # but then input genomes are columns and DB genomes are rows
            dfs = [load_skani_as_pwd(dist_file, fix_refs=fix_mash_id, fix_queries=fix_skani_id(input_extension),
                                     refs=sgb_to_genomes_db[sgb_id])
                   for dist_file in dist_files_row]

            dfs_ani, dfs_afr = zip(*dfs)
            df_ani = pd.concat(dfs_ani, axis=0).fillna(0)
            df_afr = pd.concat(dfs_afr, axis=0).fillna(0)

            closest_anis = df_ani.max().rename('closest_genome_ani')
            avg_anis = df_ani.mean().rename('sgb_avg_ani')
            max_afrs = df_afr.max().rename('max_afr')

            df_partial_results = pd.concat([closest_anis, avg_anis, max_afrs], axis=1).rename_axis('genome')\
                .reset_index()
            df_partial_results['sgb_id'] = sgb_id
            df_partial_results = df_partial_results.loc[df_partial_results['closest_genome_ani'] > 0]
            df_partial_results = df_partial_results[['genome', 'sgb_id', 'closest_genome_ani', 'sgb_avg_ani',
                                                     'max_afr']]

            dfs_partial.append(df_partial_results)

    df_long_results = pd.concat(dfs_partial, axis=0)

    if sgb_to_mp is not None:
        df_long_results['sgb_is_profilable'] = df_long_results['sgb_id'].isin(set(sgb_to_mp.keys()))
    df_long_results = df_long_results.sort_values('sgb_avg_ani', ascending=False)

    long_results_file_lock = path_get_lock(long_results_file)
    long_results_file_lock.touch()
    df_long_results.to_csv(long_results_file, sep='\t', index=False)
    long_results_file_lock.unlink()

    return df_long_results


def cluster_unassigned_genomes(nproc_cpu, work_dir, genomes_unassigned, input_skani_sketches_dir, input_extension):
    clustering_tmp_file = work_dir / 'clustering_tmp1.tsv'
    if exists_with_lock(clustering_tmp_file):
        info('Reusing clustering file of unassigned genomes')
        s_sgb_clustering = load_pandas_series(clustering_tmp_file)
    else:
        info('Clustering of unassigned genomes')
        genomes_unassigned = sorted(genomes_unassigned)
        input_unassigned_pastes_dir = work_dir / 'skani_pastes_in_unassigned'
        dists_dir_in = work_dir / 'skani_dists_in_unassigned'

        info('  Skani pasting the unassigned genomes')
        input_unassigned_pastes_dir.mkdir(exist_ok=True)
        input_unassigned_sketches = [input_skani_sketches_dir / f'{g}.sketch' for g in genomes_unassigned]
        input_unassigned_pastes = skani_paste(input_unassigned_sketches, input_unassigned_pastes_dir,
                                              chunk_size=CHUNK_SIZE_IN)

        info('  Calculating the big input-vs-input distance matrix')
        dists_dir_in.mkdir(exist_ok=True)
        skani_triangle_big(input_unassigned_pastes, dists_dir_in, nproc_cpu)

        info('  Loading the big input-vs-input distance matrix')
        pwd = load_big_triangle_skani(input_unassigned_pastes, dists_dir_in, input_extension)
        pwd = pwd.loc[genomes_unassigned, genomes_unassigned]

        info('  Running clustering of unassigned genomes into new SGBs')
        s_clusters = cluster_skani_pwd(pwd, SGB_ASSIGNMENT_ANI)
        s_sgb_clustering = (NOVEL_SGB_PREFIX + s_clusters.map(str)).rename('sgb_id').rename_axis('genome')

        clustering_tmp_file_lock = path_get_lock(clustering_tmp_file)
        clustering_tmp_file_lock.touch()
        s_sgb_clustering.to_csv(clustering_tmp_file, sep='\t')
        clustering_tmp_file_lock.unlink()

    return s_sgb_clustering


def assign_ggb_fgb(nproc_io, s_sgb_clustering, dist_files_centroids, centroid_to_sgb, sgb_to_ggb, ggb_to_fgb):
    dfs = []
    for df in run_parallel(lambda dist_f: load_pwd_pandas(dist_f), dist_files_centroids, nproc=nproc_io, ordered=False,
                           return_gen=True):
        df.index = df.index.map(fix_mash_id)  # centroids M-ids
        dfs.append(df.loc[:, df.columns.intersection(s_sgb_clustering.index)])

    df_genomes_centroid = pd.concat(dfs, axis=1).T

    assert sorted(df_genomes_centroid.columns) == sorted(centroid_to_sgb.index)
    assert sorted(df_genomes_centroid.index) == sorted(s_sgb_clustering.index)

    df_sgb_centroid = df_genomes_centroid.groupby(s_sgb_clustering).mean()

    idx_ggbs = df_sgb_centroid.columns.map(centroid_to_sgb).map(sgb_to_ggb)
    idx_fgbs = idx_ggbs.map(ggb_to_fgb)

    df_sgb_ggb = df_sgb_centroid.T.groupby(idx_ggbs).mean().T
    df_sgb_fgb = df_sgb_centroid.T.groupby(idx_fgbs).mean().T
    assignment_sgb_to_ggb = df_sgb_ggb.idxmin(axis=1)
    assignment_sgb_to_ggb_d = df_sgb_ggb.min(axis=1)


    assignment_sgb_to_ggb = assignment_sgb_to_ggb[assignment_sgb_to_ggb_d <= GGB_ASSIGNMENT_DIST]

    assignment_sgb_to_fgb_1 = assignment_sgb_to_ggb.map(ggb_to_fgb)  # assigned GGB: use existing FGB

    # Unassigned GGB: get the closest FGB and assign if <= 0.3 MASH

    df_sgb_fgb_u = df_sgb_fgb.loc[df_sgb_fgb.index.difference(assignment_sgb_to_ggb.index)]
    if len(df_sgb_fgb_u) > 0:
        assignment_sgb_to_fgb_2 = df_sgb_fgb_u.idxmin()
        assignment_sgb_to_fgb_d = df_sgb_fgb_u.min()
        assignment_sgb_to_fgb_2 = assignment_sgb_to_fgb_2[assignment_sgb_to_fgb_d <= FGB_ASSIGNMENT_DIST]


        assignment_sgb_to_fgb = pd.concat([assignment_sgb_to_fgb_1, assignment_sgb_to_fgb_2])
    else:
        assignment_sgb_to_fgb = assignment_sgb_to_fgb_1

    return assignment_sgb_to_ggb, assignment_sgb_to_fgb


def assign_phyla(nproc_io, nproc_cpu, input_directory, input_extension, work_dir, database_path, genomes_without_family,
                 centroid_to_sgb, df_sgb_sgb):
    centroids_k = df_sgb_sgb.loc[
        (df_sgb_sgb['Unknown'] == "kSGB") &
        (df_sgb_sgb['Assigned taxonomic ID'].map(lambda x: not pd.isna(x) and (x.split('|')[1] != ""))),
        'SGB centroid']

    prodigal_input_dir = work_dir / 'prodigal_input'
    prodigal_output_dir = work_dir / 'prodigal_output'
    prodigal_input_dir.mkdir(exist_ok=True)
    prodigal_output_dir.mkdir(exist_ok=True)

    info(f'  Running prodigal on {len(genomes_without_family)} input genomes')

    def run_prodigal(g):
        prodigal_output_file = prodigal_output_dir / f'{g}.faa'
        fna_file = input_directory / f'{g}{input_extension}'
        c_ = f'prodigal -c -m -p meta -i {fna_file} -a {prodigal_output_file}'

        run_command_with_lock(c_, prodigal_output_file, shell=False)

    run_parallel(run_prodigal, genomes_without_family, nproc_cpu, ordered=False)

    genome_to_faa = {}
    for genome in genomes_without_family:
        faa_path = prodigal_output_dir / f'{genome}.faa'
        if not faa_path.is_file():
            error(f'Prodigal-predicted protein sequence not found for {genome}: {faa_path}')
            error('  Will skip this genome from phylum assignment')
        else:
            genome_to_faa[genome] = faa_path

    info('  MASH sketching the proteomes')
    sketch_dir = work_dir / 'aa_sketches'
    sketch_dir.mkdir(exist_ok=True)

    genome_to_faa = pd.Series(genome_to_faa).map(pathlib.Path)
    assert set(genome_to_faa.index).issubset(set(genomes_without_family))
    mash_sketch_aa(genome_to_faa, sketch_dir, nproc_cpu, nproc_io)

    genomes_u = genome_to_faa.index

    info('  Pasting the proteome sketches')
    # TODO: we don't really need to paste, we can just use the individual sketches
    sketches_u = str(sketch_dir) + '/' + genomes_u + '.msh'

    pastes_dir_u = work_dir / 'aa_pastes'
    dists_dir = work_dir / 'aa_dists'
    pastes_dir_u.mkdir(exist_ok=True)
    dists_dir.mkdir(exist_ok=True)
    pastes_u = mash_paste(sketches_u, pastes_dir_u, nproc_io, chunk_size=CHUNK_SIZE_AA)

    pastes_k = sorted((database_path / 'centroid_mash_aa_pastes').iterdir())

    info('  Building the pwd genomes without family VS. SGBs centroids with phyla')
    dist_files = mash_dist_block(pastes_k, pastes_u, dists_dir, nproc_cpu)

    info('  Loading the big pwd')
    df_big_pwd = load_mash_dist_block(dist_files, n_rows=len(genomes_u), n_columns=len(centroids_k))
    df_big_pwd.columns = df_big_pwd.columns.map(fix_mash_id)  # centroid M-ids
    assert set(df_big_pwd.index) == set(genomes_u)
    assert set(df_big_pwd.columns) == set(centroids_k)

    info('  Getting the taxonomy assignments by the minimum distance')
    min_dists: pd.Series = df_big_pwd.T.min()
    df_big_mask = (df_big_pwd.T == min_dists).T
    u_taxonomies = df_sgb_sgb.loc[df_big_mask.columns.map(centroid_to_sgb), 'Assigned taxonomy']
    all_taxs = df_big_mask.apply(lambda row_: list(u_taxonomies.values[row_.values]), axis=1).rename('all_taxs')

    def majority_voting_at_other(taxs):
        taxs_kp = ['|'.join(y.split('|')[:2]) for y in taxs]
        taxs_kp_most_common = Counter(taxs_kp).most_common()
        return ','.join([x[0] for x in taxs_kp_most_common])


    majority_taxs = all_taxs.map(majority_voting_at_other)

    df = majority_taxs.rename('closest_phylum').to_frame().join(min_dists.rename('closest_phylum_dist'))

    return df


def get_remote_dbs():
    # TODO: figure out where to download the DB list file
    #  if running with --list_databases, work_dir is not required ==> use some local tmp file and then remove it later?
    db_list_path = os.path.basename(urllib.parse.urlparse(DB_LIST_URL).path)
    download(DB_LIST_URL, db_list_path, overwrite=True, quiet=True)
    return pd.read_csv(db_list_path, sep='\t', index_col=0)


def list_databases(args, df_dbs=None):
    if args.database_folder.is_dir():
        local_databases = [x.name for x in args.database_folder.iterdir() if x.is_dir()]
    else:
        local_databases = []
    info(f'There are {len(local_databases)} databases downloaded in {args.database_folder}:')
    info('  ' + ', '.join(local_databases))

    if df_dbs is None:
        info('Getting the list of databases from the remote server')
        df_dbs = get_remote_dbs()
    remote_databases = [x for x in df_dbs.index if 'tutorial' not in x]

    remote_databases_extra = [x for x in remote_databases if x not in local_databases]

    info(f'There are {len(remote_databases_extra)} additional databases available to download:')
    info('  ' + ', '.join(remote_databases_extra))


def download_db(args):
    df_dbs = get_remote_dbs()
    if args.database not in df_dbs.index:
        error(f'The specified database {args.database} not found')
        error('Listing the available databases')
        list_databases(args, df_dbs)
        exit(1)

    error(f'Downloading the database {args.database} from the remote server')
    url_db = df_dbs.loc[args.database, 'url_db']
    url_md5 = df_dbs.loc[args.database, 'url_md5']

    args.database_folder.mkdir(exist_ok=True)

    path_db = args.database_folder / os.path.basename(urllib.parse.urlparse(url_db).path)
    path_md5 = args.database_folder / os.path.basename(urllib.parse.urlparse(url_md5).path)

    if not path_db.name != f'{args.database}.tar':
        error(f'The database file {path_db.name} is not recognized', do_exit=True)
    if path_md5.name != f'{args.database}.md5':
        error(f'The database file {path_md5.name} is not recognized', do_exit=True)

    download(url_db, path_db)
    download(url_md5, path_md5)

    info('Extracting the database tar file')
    database_path = args.database_folder / args.database
    database_path.mkdir(exist_ok=True)
    tarfile_handle = tarfile.open(path_db)
    tarfile_handle.extractall(path=database_path)
    tarfile_handle.close()


def phylophlan_assign_sgbs(args):
    # prevent OpenBLAS inside MASH from creating multiple threads https://pythonspeed.com/articles/concurrency-control/
    # we parallelize it ourselves
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['GOTO_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'

    check_dependencies()

    info('Starting PhyloPhlAn assign SGBs')

    work_dir = args.output_directory / WORK_DIR_NAME

    if work_dir.exists() and args.clean:
        info('Cleaning the temporary directory')
        shutil.rmtree(work_dir)

    work_dir.mkdir(exist_ok=True)

    # >>> Check the input #
    genome_names, input_directory, input_extension = process_input(args.input_directory, args.input_extension,
                                                                   work_dir, args.nproc_io)

    # >>> Check the database <<< #
    database_path = args.database_folder / args.database

    if not database_path.is_dir():
        download_db(args)


    info('Loading database files')
    sgb_table_path = database_path / f'{args.database}.txt.bz2'
    profilable_sgbs_file = database_path / 'profilable_sgbs.tsv'
    sgb_to_genomes_path = database_path / 'sgb_to_genomes.tsv'

    df_sgb = load_sgb_txt(sgb_table_path)

    df_sgb_sgb = df_sgb.loc['SGB']
    df_sgb_ggb = df_sgb.loc['GGB']
    df_sgb_fgb = df_sgb.loc['FGB']

    df_sgb_sgb.index = 'SGB' + df_sgb_sgb.index.map(str)
    df_sgb_ggb.index = 'GGB' + df_sgb_ggb.index.map(str)
    df_sgb_fgb.index = 'FGB' + df_sgb_fgb.index.map(str)

    ggb_to_sgbs = dict(zip(df_sgb_ggb.index, df_sgb_ggb['List of reconstructed genomes'].str.split(',')))
    fgb_to_ggbs = dict(zip(df_sgb_fgb.index, df_sgb_fgb['List of reconstructed genomes'].str.split(',')))
    sgb_to_ggb = pd.Series({sgb: ggb for ggb, sgbs in ggb_to_sgbs.items() for sgb in sgbs})
    ggb_to_fgb = pd.Series({ggb: fgb for fgb, ggbs in fgb_to_ggbs.items() for ggb in ggbs})
    centroid_to_sgb = pd.Series(dict(zip(df_sgb_sgb['SGB centroid'], df_sgb_sgb.index)))
    centroids = centroid_to_sgb.index

    if profilable_sgbs_file.exists():
        sgb_to_mp = load_pandas_series(profilable_sgbs_file).to_dict()
    else:
        sgb_to_mp = None

    sgb_to_genomes_db = pd.read_csv(sgb_to_genomes_path, sep='\t', index_col=0).squeeze('columns').str.split(',')


    # >>> Sketching <<< #
    info('Generating MASH sketches for the input genomes')
    input_mash_sketches_dir = work_dir / 'input_sketches_mash'
    input_mash_sketches_dir.mkdir(exist_ok=True)
    input_mash_sketches = mash_sketch(genome_names, input_extension, input_directory,
                                      input_mash_sketches_dir, args.nproc_cpu, args.nproc_io)

    info('Generating skani sketches for the input genomes')
    input_skani_sketches_dir = work_dir / 'input_sketches_skani'
    input_skani_sketches_dir.mkdir(exist_ok=True)
    skani_sketch(genome_names, input_extension, input_directory, input_skani_sketches_dir, args.nproc_cpu,
                 args.nproc_io)

    # >>> Prefiltering <<< #
    info('Prefiltering SGBs')
    df_genome_centroid_pairs, dist_files_centroids = prefiltering(work_dir, args.nproc_io, args.nproc_cpu,
                                                                  database_path, centroids, centroid_to_sgb,
                                                                  input_mash_sketches, genome_names)
    genomes_near = df_genome_centroid_pairs['genome'].unique()
    genomes_far = set(genome_names).difference(set(genomes_near))

    all_centroid_hits = df_genome_centroid_pairs['centroid'].unique()
    info(f'  There are {len(genomes_far)} genomes far and {len(genomes_near)} genomes '
         f'near {len(all_centroid_hits)} centroids')

    sgb_to_genomes_in = df_genome_centroid_pairs.set_index('genome').groupby('sgb_id').groups


    # >>> Distancing against prefiltered SGBs <<< #
    df_long_results = dist_against_sgbs(args.nproc_io, args.nproc_cpu, work_dir, database_path, sgb_to_genomes_in,
                                        sgb_to_genomes_db, input_skani_sketches_dir, sgb_to_mp, input_extension)

    genome_to_max_overlap = df_long_results.groupby('genome').agg({'max_afr': 'max'})['max_afr']


    df_top_results_single = df_long_results.drop_duplicates(subset=['genome'])  # the df is sorted by decreasing avg ANI
    df_genome_assignment = df_top_results_single[df_top_results_single['sgb_avg_ani'] >= SGB_ASSIGNMENT_ANI]\
        .set_index('genome')

    if sgb_to_mp is not None:
        df_closest_profilable = df_long_results[df_long_results['sgb_is_profilable']]\
            .drop_duplicates(subset=['genome']).set_index('genome')
        df_closest_profilable = df_closest_profilable[df_closest_profilable['sgb_avg_ani']
                                                      >= SGB_CLOSEST_PROFILABLE_MIN_ANI]
    else:
        df_closest_profilable = None


    df_genome_assignment = df_genome_assignment \
        .join(sgb_to_ggb.rename('ggb_id'), on='sgb_id') \
        .join(ggb_to_fgb.rename('fgb_id'), on='ggb_id')

    assigned_genomes_set = set(df_genome_assignment.index)
    genomes_near_unassigned = set(genomes_near).difference(set(df_genome_assignment.index))

    info(f'  There are {len(assigned_genomes_set)} assigned genomes,'
         f' {len(genomes_near_unassigned)} unassigned but near genomes, {len(genomes_far)} far genomes')

    genomes_unassigned = set(genome_names).difference(assigned_genomes_set)

    if len(genomes_unassigned) > 0:
        # >>> Clustering unassigned genomes <<< #
        assignment_genomes_to_new_sgbs = cluster_unassigned_genomes(args.nproc_cpu, work_dir, genomes_unassigned,
                                                                    input_skani_sketches_dir, input_extension)
        novel_sgbs = set(assignment_genomes_to_new_sgbs.unique())
        info(f'  There are {len(novel_sgbs)} SGBs '
             f'created by {len(assignment_genomes_to_new_sgbs)} genomes')



        # >>> Assigning GGBs and SGBs of newly created SGBs <<< #
        info('Assigning GGBs and FGBs to novel SGBs')
        assignment_sgb_to_ggb, assignment_sgb_to_fgb = assign_ggb_fgb(args.nproc_io, assignment_genomes_to_new_sgbs,
                                                                      dist_files_centroids, centroid_to_sgb, sgb_to_ggb,
                                                                      ggb_to_fgb)

        df_genome_assignment_unassigned = assignment_genomes_to_new_sgbs.rename('sgb_id').to_frame() \
            .join(assignment_sgb_to_ggb.rename('ggb_id'), on='sgb_id') \
            .join(assignment_sgb_to_fgb.rename('fgb_id'), on='sgb_id')

        df_genome_assignment = pd.concat([df_genome_assignment, df_genome_assignment_unassigned])

        assert set(df_genome_assignment.index) == set(genome_names)


        # >>> Assign phyla to genomes without family <<< #
        genomes_without_family = df_genome_assignment.index[df_genome_assignment['fgb_id'].isna()]
        if len(genomes_without_family) > 0:
            info(f'Will assign phyla to {len(genomes_without_family)} genomes without family')
            df_phyla = assign_phyla(args.nproc_io, args.nproc_cpu, input_directory, input_extension, work_dir,
                                    database_path, genomes_without_family, centroid_to_sgb, df_sgb_sgb)

            df_genome_assignment = df_genome_assignment.join(df_phyla)
    else:
        novel_sgbs = set()


    # >>> Assign taxonomy <<< #


    def get_tax(row):
        if not pd.isna(row['sgb_id']) and row['sgb_id'] not in novel_sgbs:
            return df_sgb_sgb.loc[row['sgb_id'], 'Assigned taxonomy']
        if not pd.isna(row['ggb_id']):
            return df_sgb_ggb.loc[row['ggb_id'], 'Assigned taxonomy']
        if not pd.isna(row['fgb_id']):
            return df_sgb_fgb.loc[row['fgb_id'], 'Assigned taxonomy']
        return row['closest_phylum']

    df_genome_assignment['taxonomy'] = df_genome_assignment.apply(get_tax, axis='columns')

    if df_closest_profilable is not None:
        df_genome_assignment = df_genome_assignment.join(df_closest_profilable[['sgb_id', 'sgb_avg_ani']],
                                                         rsuffix='_metaphlan')
        df_genome_assignment['sgb_id_metaphlan'] = df_genome_assignment['sgb_id_metaphlan'].map(sgb_to_mp)

    df_genome_assignment = df_genome_assignment.join(genome_to_max_overlap.rename('max_overlap'))
    genomes_with_low_overlap = genome_to_max_overlap.index[genome_to_max_overlap < MAX_OVERLAP_WARNING]
    if len(genomes_with_low_overlap) > 0:
        warning(f'There are {len(genomes_with_low_overlap)} genomes that have low maximum overlap'
                f' (< {MAX_OVERLAP_WARNING}).'
                f' This could indicate low completeness of the genomes and their assignment might be incorrect.'
                f' We recommend running the assignment on genomes at least 50% complete.')

    df_genome_assignment.to_csv(args.output_directory / 'genome_assignment.tsv', sep='\t')
    info(f"Output table written to {args.output_directory / 'genome_assignment.tsv'}")

    info('End.')


def main():
    argp = read_params()
    args = check_params(argp)

    if args.list_databases:
        list_databases(args)
        return

    phylophlan_assign_sgbs(args)


if __name__ == '__main__':
    main()
