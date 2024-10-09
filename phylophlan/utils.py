import gzip
import os
import pathlib
import argparse as ap
import sys
import shlex
import shutil
import subprocess as sp
import bz2
from datetime import datetime
import itertools as it
from multiprocessing.pool import ThreadPool, Pool

import numpy as np
import pandas as pd
import scipy.spatial.distance as spd
import scipy.cluster.hierarchy as spch
from tqdm.auto import tqdm as tqdm_orig


def message(*args, f, new_line=True):
    d = datetime.now().strftime("%d/%m/%y %H:%M:%S")
    m = ' '.join(args)
    f.write("[%s] %s" % (d, m))
    if new_line:
        f.write('\n')
    f.flush()


def info(*args):
    message("Info:", *args, f=sys.stdout)


def warning(*args):
    message("Warning:", *args, f=sys.stdout)


def error(*args, do_exit=False):
    message("Error:", *args, f=sys.stderr)
    if do_exit:
        exit(1)


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / 1048576


# This is to avoid warning in pycharm that tqdm.asyncio is not an iterable
class tqdm(tqdm_orig):
    def __iter__(self):
        return tqdm_orig.__iter__(self)


class ArgumentType:
    @staticmethod
    def existing_file(path):
        path = pathlib.Path(path).resolve()
        if not os.path.isfile(path):
            raise ap.ArgumentTypeError('The file does not exist (%s).' % path)
        return path

    @classmethod
    def list_in_file(cls, path):
        if not path:
            return []

        path = cls.existing_file(path)
        with open(path) as f:
            r = [x.rstrip('\r\n') for x in f]
            r = [x for x in r if x]
            return r

    @staticmethod
    def existing_dir(path):
        path = pathlib.Path(path).resolve()
        if not path.is_dir():
            raise ap.ArgumentTypeError('The directory does not exist (%s).' % path)
        return path

    @staticmethod
    def creatable_dir(path):
        path = pathlib.Path(path).resolve()
        if path.is_dir() or path.parent.resolve().is_dir():
            return path
        raise ap.ArgumentTypeError('Neither the directory nor its parent exist (%s).' % path)

    @staticmethod
    def creatable_file(path):
        path = pathlib.Path(path).resolve()
        if path.is_file() or path.parent.resolve().is_dir():
            return path
        raise ap.ArgumentTypeError('Neither the file nor its parent exist (%s).' % path)

    @staticmethod
    def percentage(x):
        x = float(x)
        if x < 0 or x > 100:
            raise ap.ArgumentTypeError('The percentage must be between 0.0 and 100.0')

    @staticmethod
    def fraction(x):
        x = float(x)
        if x < 0 or x > 1:
            raise ap.ArgumentTypeError('The fraction must be between 0.0 and 1.0')

    @staticmethod
    def positive_int(x):
        x = int(x)
        if x <= 0:
            raise ap.ArgumentTypeError('The number must be greater than 0')
        return x


def get_threads_per_run(nproc_cpu, nproc_io):
    return max(1, nproc_cpu // nproc_io)


def openr(filepath, mode="rt"):
    """
    Wrapper for "open" and "bz2.open" used to open files in read only mode
    """
    filepath = str(filepath)
    if filepath.endswith('.bz2'):
        return bz2.open(filepath, mode)
    elif filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)


def load_pandas_series(file_path, sep='\t', index_col=0, **kwargs):
    return pd.read_csv(file_path, sep=sep, index_col=index_col, **kwargs).squeeze('columns')


def load_pandas_tsv(file_path, header=0, index_col=False, **kwargs):
    """
    Skips the lines starting with "#" except the last one which is used as a header
    :param file_path:
    :param header:
    :param bool|Sequence[int]|int index_col:
    :param kwargs: Additional arguments passed to pd.read_csv
    :return pd.DataFrame:
    """
    with openr(file_path, 'rt') as f:
        nl_positions = [f.tell()]
        line = f.readline().strip()
        while line:
            if not line.startswith('#'):
                f.seek(nl_positions[-2])
                break
            nl_positions.append(f.tell())
            line = f.readline().strip()
        else:
            print('Warning: EOL encountered')
            return None

        return pd.read_csv(f, sep='\t', skiprows=0, header=header, index_col=index_col, **kwargs)


def load_sgb_txt(sgb_file_path):
    return load_pandas_tsv(sgb_file_path, index_col=[0, 1], dtype={'ID': int})


def fix_mash_id(x):
    """
    Extract a metaref id from the full path encoded in .msh files
    :param x: The original id (full path)
    :return: The metaref id
    """
    return x.split('/')[-1].rsplit('.', maxsplit=1)[0]
    # return os.path.splitext(os.path.basename(x))[0]


def load_pwd_pandas(file_path, fix_index=False, **kwargs):
    df = pd.read_csv(file_path, sep='\t', index_col=0, header=0, **kwargs)
    if fix_index:
        df.index = df.index.map(fix_mash_id)
        df.columns = df.columns.map(fix_mash_id)

    return df


def load_skani_as_pwd(path, queries=None, refs=None):
    rows = []
    with openr(path) as f:
        header = next(f)
        header = header.split('\t', maxsplit=3)[:3]
        for line in f:
            r, q, a, _ = line.split('\t', maxsplit=3)
            a = float(a)
            rows.append([r, q, a])
    df = pd.DataFrame(rows, columns=header)
    df = df.pivot(index='Ref_file', columns='Query_file', values='ANI')
    df.index = df.index.map(fix_mash_id)
    df.columns = df.columns.map(fix_mash_id)
    df = df.reindex(index=refs, columns=queries).fillna(0)
    return df


def run_command(cmd, shell=True, **kwargs):
    """
    Runs a command and checks for exit code

    :param cmd:
    :param shell:
    """
    if shell:
        if '|' in cmd:
            cmd_s = ['/bin/bash', '-o', 'pipefail', '-c', cmd]
            shell = False
        else:
            cmd_s = cmd
    else:
        cmd_s = shlex.split(cmd)

    r = sp.run(cmd_s, shell=shell, capture_output=True, **kwargs)

    if r.returncode != 0:
        stdout = r.stdout
        stderr = r.stderr
        if isinstance(stdout, bytes):
            stdout = stdout.decode()
        if isinstance(stderr, bytes):
            stderr = stderr.decode()

        error('Execution failed for command', cmd)
        print('stdout: ')
        print(stdout)
        print('stderr: ')
        print(stderr)
        error('Exiting')
        exit(1)

    return r


def path_get_lock(path):
    path = pathlib.Path(path)
    return path.with_suffix(path.suffix + '.lock')


def exists_with_lock(path):
    path = pathlib.Path(path)
    lock_file = path_get_lock(path)
    if lock_file.exists():  # unfinished dist file => remove and run again
        info(f'Removing unfinished file {path}')
        path.unlink(missing_ok=True)
        lock_file.unlink()
    elif path.exists():
        return True
    return False


def run_command_with_lock(cmd, output_file, *args, **kwargs):
    lock_file = path_get_lock(output_file)
    if exists_with_lock(output_file):
        return
    lock_file.touch()
    r = run_command(cmd, *args, **kwargs)
    lock_file.unlink()
    return r


def run_parallel_gen(f, f_args, nproc, f_const=None, chunksize=1, ordered=True, star=False, processes=False):
    if chunksize == 'auto':
        if len(f_args) / nproc < 10:
            chunksize = 1
        else:
            chunksize = 10

    if star:
        def f_starred(a):
            return f(*a)
    else:
        f_starred = f

    if processes:
        pool_class = Pool
    else:
        pool_class = ThreadPool

    def init(fu, f_const_):
        fu.constants = f_const_

    with pool_class(nproc, initializer=init, initargs=(f, f_const)) as tp:
        if ordered:
            map_f = tp.imap
        else:
            map_f = tp.imap_unordered
        results = map_f(f_starred, f_args, chunksize=chunksize)
        if len(f_args) > nproc:
            results = tqdm(results, total=len(f_args))
        for r in results:
            yield r


def run_parallel(f, f_args, nproc, f_const=None, chunksize=1, ordered=True, star=False, return_iter=False,
                 processes=False):
    """
    :param function f: single-argument function
    :param Sequence f_args: sequence of arguments to pass to f
    :param int nproc:
    :param f_const: An object shared across workers accessible in the function as f.constants
    :param int|str chunksize: positive integer or 'auto'
    :param ordered: Whether to use imap (preserve order) or imap_unordered
    :param star: If set to True, each of f_args is a Sequence of positional arguments to pass to f
    :param return_iter: Whether to return an iterator instead of list
    :param processes: Whether to spawn processes instead of just Threads
    :return list: list of return values from f
    """
    if processes and star:
        raise Exception('The arguments processes and star are incompatible (for now)')

    r = run_parallel_gen(f, f_args, nproc, f_const, chunksize, ordered, star, processes)
    if return_iter:
        return r
    else:
        return list(r)


def mash_sketch(genome_ids, args, genomes_dir, sketch_dir, sketch_size=10000):
    threads_per_run = get_threads_per_run(args.nproc_cpu, args.nproc_io)

    if args.input_extension.endswith('.bz2'):
        cat_f = 'bzcat'
    elif args.input_extension.endswith('.gz'):
        cat_f = 'zcat'
    else:
        cat_f = 'cat'

    commands = []
    reused = []
    sketch_files = []
    for g in genome_ids:
        genome_file = genomes_dir / f'{g}{args.input_extension}'
        sketch_file = sketch_dir / f'{g}.msh'
        sketch_files.append(sketch_file)
        if sketch_file.exists():
            reused.append(sketch_file)
            continue

        assert genome_file.exists()

        commands.append((f"{cat_f} {genome_file} | mash sketch -p {threads_per_run} -k 21 -s {sketch_size} "
                        f"-o {sketch_file} -I {g} -", sketch_file))

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(run_command_with_lock, commands, args.nproc_io, chunksize='auto', ordered=False, star=True)
    return sketch_files



def mash_sketch_aa(genome_to_faa, sketch_dir, args, sketch_size=10000):
    threads_per_run = get_threads_per_run(args.nproc_cpu, args.nproc_io)

    commands = []
    reused = []
    sketch_files = []
    for g, faa_path in genome_to_faa.items():
        if faa_path.name.endswith('.bz2'):
            cat_f = 'bzcat'
        elif faa_path.name.endswith('.gz'):
            cat_f = 'zcat'
        else:
            cat_f = 'cat'

        sketch_file = sketch_dir / f'{g}.msh'
        sketch_files.append(sketch_file)
        if sketch_file.exists():
            reused.append(sketch_file)
            continue

        assert faa_path.exists()

        commands.append((f"{cat_f} {faa_path} | mash sketch -p {threads_per_run} -a -k 9 -s {sketch_size} "
                        f"-o {sketch_file} -I {g} -", sketch_file))

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(run_command_with_lock, commands, args.nproc_io, chunksize='auto', ordered=False, star=True)
    return sketch_files


def skani_sketch(genome_ids, args, genomes_dir, sketch_dir):
    threads_per_run = get_threads_per_run(args.nproc_cpu, args.nproc_io)

    def skani_sketch_one(sketch_dir_, genome_file_, sketch_file_, g_):
        tmp_dir_ = sketch_dir_ / g_
        if tmp_dir_.exists():
            shutil.rmtree(tmp_dir_)
        run_command(f"skani sketch -t {threads_per_run} -o {tmp_dir_} {genome_file_}")
        os.rename(tmp_dir_ / f'{g_}{args.input_extension}.sketch', sketch_file_)
        os.unlink(tmp_dir_ / "markers.bin")
        os.rmdir(tmp_dir_)


    f_args = []
    reused = []
    sketch_files = []
    for g in genome_ids:
        genome_file = genomes_dir / f'{g}{args.input_extension}'
        sketch_file = sketch_dir / f'{g}.sketch'
        sketch_files.append(sketch_file)
        if sketch_file.exists():
            reused.append(sketch_file)
            continue

        assert genome_file.exists()

        if args.input_extension.endswith('.gz') or args.input_extension.endswith('.bz2'):
            raise NotImplemented('We need to implement decompression for skani sketch')


        f_args.append([sketch_dir, genome_file, sketch_file, g])

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(skani_sketch_one, f_args, args.nproc_io, star=True, ordered=False)
    return sketch_files


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def mash_paste_prepare(genome_sketches, pastes_dir, chunk_size, file_prefix):
    genome_sketches = sorted(genome_sketches)  # for reproducibility
    paste_files = []
    commands_with_output = []
    reused = []
    for k, sketches_chunk in enumerate(chunks(genome_sketches, chunk_size)):
        sketches_txt = pastes_dir / f'{file_prefix}{k}.txt'
        paste_prefix = pastes_dir / f'{file_prefix}{k}'
        with open(sketches_txt, 'w') as f:
            f.write('\n'.join(str(s) for s in sketches_chunk) + '\n')
        paste_file = pastes_dir / f'{file_prefix}{k}.msh'
        paste_files.append(paste_file)
        if paste_file.exists():
            reused.append(paste_file)
            continue

        commands_with_output.append((f"mash paste -l {paste_prefix} {sketches_txt}", paste_file))

    return commands_with_output, paste_files, reused



def mash_paste(genome_sketches, pastes_dir, nproc, chunk_size, file_prefix='paste_pt'):
    """
    Paste sketches of input genomes.
    Assumes sketches are in place.
    :param Sequence[pathlib.Path|str] genome_sketches:
    :param pathlib.Path pastes_dir:
    :param int nproc:
    :param int chunk_size:
    :param file_prefix:
    :returns list: List of mash pastes created
    """

    commands, paste_files, reused = mash_paste_prepare(genome_sketches, pastes_dir, chunk_size, file_prefix)

    if reused:
        info('Reused %d existing paste files' % len(reused))

    run_parallel(run_command, commands, nproc, ordered=False)
    return paste_files


def skani_paste(sketches, pastes_dir, chunk_size, file_prefix='paste'):
    paste_files = []
    for i, sketches_chunk in enumerate(chunks(sketches, chunk_size)):
        paste_file = pastes_dir / f'{file_prefix}_pt{i}.txt'
        with open(paste_file, 'w') as f:
            f.write('\n'.join(map(str, sketches_chunk)))

        paste_files.append(paste_file)

    return paste_files


def mash_dist_block(mash_pastes_1, mash_pastes_2, dists_dir, nproc, progress_bar=True):
    """

    :param mash_pastes_1:
    :param mash_pastes_2:
    :param dists_dir:
    :param nproc:
    :param progress_bar:
    """

    if progress_bar:
        mash_pastes_1 = tqdm(mash_pastes_1)

    dist_files = []
    for mp1 in mash_pastes_1:
        dist_filename = mp1.stem + '.tsv.gz'
        dist_file = dists_dir / dist_filename
        cmd = f"mash dist -p {nproc} -t {mp1} {' '.join(map(str, mash_pastes_2))} | gzip > {dist_file}"
        run_command_with_lock(cmd, output_file=dist_file, shell=True)
        dist_files.append(dist_file)

    return dist_files


def skani_dist_block(skani_pastes_1, skani_pastes_2, dists_dir, nproc, progress_bar=True):
    """

    :param skani_pastes_1:
    :param skani_pastes_2:
    :param dists_dir:
    :param nproc:
    :param progress_bar:
    """

    if progress_bar:
        skani_pastes_1 = tqdm(skani_pastes_1)

    dist_files = []
    for p1 in skani_pastes_1:
        for p2 in skani_pastes_2:
            dist_filename = p1.stem + '__' + p2.stem + '.tsv.gz'
            dist_file = dists_dir / dist_filename
            cmd = f"skani dist -t {nproc} --ql {p1} --rl {p2} | gzip > {dist_file}"
            run_command_with_lock(cmd, output_file=dist_file, shell=True)
            dist_files.append(dist_file)

    return dist_files


def skani_triangle_big(skani_pastes, dists_dir, nproc):
    for mp in skani_pastes:  # diagonal triangles
        dist_file = dists_dir / (mp.stem + '.tri.gz')
        cmd = f'skani triangle -t {nproc} -l {mp} | gzip > {dist_file}'
        run_command_with_lock(cmd, dist_file)

    for mp1, mp2 in it.combinations(skani_pastes, 2):  # pairwise blocks
        dist_file = dists_dir / (mp1.stem + '__' + mp2.stem + '.tsv.gz')
        cmd = f'skani dist -t {nproc} --ql {mp1} --rl {mp2} | gzip > {dist_file}'
        run_command_with_lock(cmd, dist_file)


def load_big_triangle_common(index_pastes, parse_values, pastes, dists_dir, paste_to_offset):
    total_n = len(index_pastes)
    if parse_values:
        matrix = np.zeros((total_n, total_n), dtype=float)
    else:
        matrix = np.zeros((total_n, total_n), dtype=object)

    # load the diagonal
    index_diagonal = []
    for p_i, paste in enumerate(pastes):
        dist_file = dists_dir / (paste.stem + '.tri.gz')
        offset = paste_to_offset[paste]
        with gzip.open(dist_file, 'rt') as f:
            first_line = next(f)
            n_lines = int(first_line.strip())

            for i, line in enumerate(f):
                if i == 0:  # first line does not have any numbers
                    # index_pastes.append(line.rstrip('\n'))
                    idx = line.rstrip('\n')
                    index_diagonal.append(idx)
                    # assert idx == index_pastes[offset + i], f"{idx} != {index_pastes[offset + i]}, {mash_pastes}"
                    continue

                ix = line.index('\t')
                idx = line[:ix]
                index_diagonal.append(idx)
                # assert index_pastes[offset + i] == idx, f"{idx} != {index_pastes[offset + i]}"

                if parse_values:
                    values = np.fromstring(line[ix + 1:], sep='\t', dtype=float)
                else:
                    values = line[ix + 1:].rstrip('\n').split('\t')

                matrix[offset + i, offset: offset + i] = values
                matrix[offset: offset + i, offset + i] = values

            offset += n_lines

    index_diagonal = [fix_mash_id(x) for x in index_diagonal]

    assert set(index_pastes) == set(index_diagonal)

    return matrix, index_diagonal


def load_big_triangle_skani(skani_pastes, dists_dir, parse_values=True):
    skani_pastes = list(skani_pastes)  # fix the ordering

    index_pastes = []
    paste_to_offset = {}
    offset = 0
    for mp in skani_pastes:

        assert str(mp).endswith('.txt')
        with open(mp) as f:
            mp_index = f.read().strip().split('\n')

        index_pastes.extend(mp_index)
        paste_to_offset[mp] = offset
        offset += len(mp_index)

    index_pastes = [fix_mash_id(x) for x in index_pastes]

    matrix, _ = load_big_triangle_common(index_pastes, parse_values, skani_pastes, dists_dir, paste_to_offset)

    idx_to_offset = dict(zip(index_pastes, range(len(index_pastes))))

    # load the blocks
    for mp1, mp2 in it.combinations(skani_pastes, 2):  # pairwise blocks
        dist_file = dists_dir / (mp1.stem + '__' + mp2.stem + '.tsv.gz')

        with openr(dist_file) as f:
            next(f)
            for line in f:
                r, q, a, _ = line.split('\t', maxsplit=3)
                r = fix_mash_id(r)
                q = fix_mash_id(q)
                a = float(a)
                i = idx_to_offset[r]
                j = idx_to_offset[q]
                matrix[i, j] = a
                matrix[j, i] = a

    np.fill_diagonal(matrix, 100.0)
    df = pd.DataFrame(matrix, index=index_pastes, columns=index_pastes)
    return df


def load_dist_block(dist_files, n_rows, n_columns, parse_values=True, progress_bar=False):
    if parse_values:
        a = np.full((n_rows, n_columns), fill_value=-1, dtype=float)
    else:
        a = np.full((n_rows, n_columns), fill_value=-1, dtype=object)

    if progress_bar:
        dist_files = tqdm(dist_files)

    index = None
    header = []
    column_offset = 0
    for dist_file in dist_files:
        with openr(dist_file) as f:
            current_header = next(f).rstrip('\n').split('\t')
            current_header = current_header[1:]
            header.extend(current_header)
            n_columns_current = len(current_header)

            current_index = []
            for row_i, line in enumerate(f):
                i = line.index('\t')
                idx = line[:i]
                current_index.append(idx)
                if parse_values:
                    values = np.fromstring(line[i + 1:], sep='\t', dtype=float)
                else:
                    values = line[i + 1:].rstrip('\n').split('\t')
                assert len(values) == n_columns_current
                a[row_i, column_offset:(column_offset+n_columns_current)] = values
            current_index = pd.Index(current_index)
            # current_index = current_index.map(fix_mash_id)
            if index is None:
                index = current_index
            assert (index == current_index).all()
            column_offset += n_columns_current

    header = pd.Index(header)
    header = header.map(fix_mash_id)
    df = pd.DataFrame(a, index=index, columns=header)
    assert not (df < 0).any().any()  # all values were filled
    return df


def cluster_skani_pwd(pwd, t=5):
    genomes_to_be_clustered = pwd.index
    condensed_dist_mat = spd.squareform(100 - pwd.values)

    if len(genomes_to_be_clustered) > 1:
        dendro = spch.linkage(condensed_dist_mat, method='average')
        clusters = spch.fcluster(dendro, t, criterion='distance')
    else:  # either one or zero genomes (the distance matrix is undefined)
        # assign all to a cluster "1"
        clusters = [1] * len(genomes_to_be_clustered)

    s_clusters = pd.Series(clusters, index=genomes_to_be_clustered)
    return s_clusters
