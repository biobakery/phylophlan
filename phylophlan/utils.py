import contextlib
import functools
import hashlib
import os
import sys
import tarfile
from datetime import datetime
import gzip
import bz2
import pathlib
import argparse as ap
import shlex
import shutil
import subprocess as sp
import urllib.request
from multiprocessing.pool import ThreadPool, Pool
from typing import Sequence, Generator, IO
import itertools as it

import numpy as np
import pandas as pd
import scipy.spatial.distance as spd
import scipy.cluster.hierarchy as spch
from tqdm.auto import tqdm as tqdm_orig


MD5_CHUNK_SIZE = 4096


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


def tqdm_bytes(*args, **kwargs):
    return tqdm(*args, unit='B', unit_scale=True, unit_divisor=1024, **kwargs)


class TqdmDownload:
    """
    To use as progress bar when downloading a file
    Modified from https://github.com/tqdm/tqdm#hooks-and-callbacks
    """
    def __init__(self):
        self.pb = tqdm_bytes(miniters=1)


    def __call__(self, n_blocks, block_size, tot_size):
        downloaded = n_blocks * block_size
        if tot_size is not None and tot_size > 0:
            self.pb.total = tot_size
            if downloaded > tot_size:  # counting the last block the downloaded goes beyond the total size, so shrink it
                downloaded = tot_size

        delta = downloaded - self.pb.n
        self.pb.update(delta)  # also sets self.n = downloaded


class ArgumentType:
    @staticmethod
    def existing_file(path):
        path = pathlib.Path(path).resolve()
        if not os.path.isfile(path):
            raise ap.ArgumentTypeError(f'The file does not exist ({path})')
        return path

    @classmethod
    def list_in_file(cls, path):
        if path is None:
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
            raise ap.ArgumentTypeError(f'The directory does not exist ({path})')
        return path

    @staticmethod
    def creatable_dir(path):
        path = pathlib.Path(path).resolve()
        if path.is_dir() or path.parent.is_dir():
            return path
        raise ap.ArgumentTypeError(f'Neither the directory nor its parent exist ({path})')

    @staticmethod
    def creatable_file(path):
        path = pathlib.Path(path).resolve()
        if path.is_file() or path.parent.is_dir():
            return path
        raise ap.ArgumentTypeError(f'Neither the file nor its parent exist ({path})')

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


def delete_file_or_dir(path):
    """

    :param pathlib.Path path:
    :return:
    """
    if path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
    else:
        raise Exception(f"Not a file neither directory {path}")



class FileInProgress(contextlib.AbstractContextManager):
    def __init__(self, path_target):
        """

        :param pathlib.Path|str path_target:
        """
        self.path_target = pathlib.Path(path_target)
        self.path_tmp = self.path_target.with_name('~' + self.path_target.name)

    def __enter__(self):
        self.clean()
        return self.path_tmp

    def __exit__(self, exception_type, exception_value, exception_traceback):
        if exception_type is None:
            self.path_tmp.rename(self.path_target)
        else:
            self.clean()


    def clean(self):
        if self.path_target.exists():
            delete_file_or_dir(self.path_target)

        if self.path_tmp.exists():
            delete_file_or_dir(self.path_tmp)


def get_threads_per_run(nproc_cpu, nproc_io):
    """
    Assuming we run a program nproc_io times in parallel, how many threads we should give to each call not to exceed
    nproc_cpu
    """
    return max(1, nproc_cpu // nproc_io)


def openr(filepath, mode="rt"):
    """
    Wrapper for "open", "bz2.open" and "gzip.open" used to open files in read text mode

    :param pathlib.Path|str filepath:
    :param str mode:
    :return:
    :rtype: IO
    """
    filepath = str(filepath)
    if filepath.endswith('.bz2'):
        return bz2.open(filepath, mode)
    elif filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)


def decompress_file(input_file, output_file, open_f):
    """
    Decompresses the `input_file` to `output_file` using `open_f` which will typically be gzip.open or bz2.open

    :param (pathlib.Path|str) input_file:
    :param (pathlib.Path|str) output_file:
    :param function open_f:
    :return:
    """
    with open_f(input_file, 'rb') as f_in, open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def download(url, target_path, overwrite=False, quiet=False):
    """

    :param str url: URL to download from
    :param pathlib.Path|str target_path: the target file path to download to
    :param bool overwrite: If set to False and the file is already present, it will skip the download
    :param bool quiet: If set to True it will not show progress bar or messages unless there's an error
    :return:
    """


    if os.path.exists(target_path) and not overwrite:
        if not quiet:
            info(f'Skipping already downloaded file {target_path}')
        return

    if not quiet:
        info(f'Downloading {target_path}')
        if os.path.exists(target_path):
            info(f'Will overwrite existing file {target_path}')

    try:
        if not quiet:
            rh = TqdmDownload()
        else:
            rh = None

        with FileInProgress(target_path) as target_path_tmp:
            urllib.request.urlretrieve(url, target_path_tmp, reporthook=rh)

    except EnvironmentError:
        error(f'Error downloading {target_path} from "{url}"', do_exit=True)


def check_md5(target_file, md5_file):
    """

    :param pathlib.Path target_file:
    :param pathlib.Path md5_file:
    :return:
    """
    with open(md5_file) as f:
        md5_md5, fname = f.read().strip().split()

    if target_file.name != fname:
        error(f'The file name {target_file.name} does not correspond to the name in the md5 file {fname}')
        return False

    target_file_size = target_file.stat().st_size
    n_chunks = (target_file_size - 1) // MD5_CHUNK_SIZE + 1

    # compute MD5 of .tar
    hash_md5 = hashlib.md5()
    with open(target_file, "rb") as f:
        for chunk in tqdm_bytes(iter(lambda: f.read(MD5_CHUNK_SIZE), b""), total=n_chunks):
            hash_md5.update(chunk)

    md5_target = hash_md5.hexdigest()
    assert len(md5_target) == 32

    return md5_target == md5_md5


def extract_tar(tar_path, target_path):
    """

    :param pathlib.Path tar_path:
    :param pathlib.Path target_path:
    :return:
    """
    tar_file_size = tar_path.stat().st_size

    with FileInProgress(target_path) as target_path_tmp, \
            tarfile.open(tar_path) as tf, \
            tqdm_bytes(total=tar_file_size) as pb:
        target_path_tmp.mkdir(exist_ok=True)
        for m in tf:
            tf.extract(m, target_path_tmp)
            pb.update(m.size)


def load_pandas_series(file_path, sep='\t', index_col=0, **kwargs) -> pd.Series:
    """
    Load a two-column file as pandas.Series

    :param file_path:
    :param sep:
    :param index_col:
    :param kwargs:
    :return:
    """
    return pd.read_csv(file_path, sep=sep, index_col=index_col, **kwargs).squeeze('columns')


def load_pandas_tsv(file_path, header=0, index_col=False, **kwargs) -> pd.DataFrame:
    """
    Skips the lines starting with "#" except the last one which is used as a header

    :param file_path:
    :param header:
    :param bool|Sequence[int]|int index_col:
    :param kwargs: Additional arguments passed to pd.read_csv
    :return:
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
            return pd.DataFrame()

        return pd.read_csv(f, sep='\t', skiprows=0, header=header, index_col=index_col, **kwargs)


def load_sgb_txt(sgb_file_path):
    return load_pandas_tsv(sgb_file_path, index_col=[0, 1], dtype={'ID': int})


def fix_mash_id(x):
    """
    Extract a metaref id from the full path encoded in .msh files

    :param x: The original id (full path)
    :return: The metaref id
    """
    x = x.split('/')[-1]
    assert x.startswith('M1')
    x = x.rsplit('.', maxsplit=1)[0]
    return x
    # return os.path.splitext(os.path.basename(x))[0]


def fix_skani_id(genome_extension):
    def fix_skani_id_inner(x):
        x = x.split('/')[-1]
        if len(genome_extension) > 0:
            assert x.endswith(genome_extension)
            return x[:-len(genome_extension)]
        else:
            return x

    return fix_skani_id_inner


def run_command(cmd, shell=True, **kwargs):
    """
    Runs a command and checks for exit code and throws an error in case of failure

    :param str cmd: The command to execute
    :param bool shell: Whether to execute using shell: allows for pipes, glob and other advance features
    :param kwargs: Additional arguments passed to subprocess.run function
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


def run_command_with_output_file(cmd, output_file, *args, **kwargs):
    """
    Runs a command, but checks for tries to skip if the output file is already present, checking for lock files.

    :param str cmd:
    :param pathlib.Path|str output_file:
    :param args:
    :param kwargs:
    :return:
    """
    with FileInProgress(output_file) as output_file_tmp:
        cmd_filled = cmd.format(output_file_tmp)
        return run_command(cmd_filled, *args, **kwargs)



def expand_star_args(args, f):
    return f(*args)


def run_parallel_gen(f, f_args, nproc, f_const=None, chunksize=1, ordered=True, star=False, processes=False):
    if chunksize == 'auto':
        if len(f_args) / nproc < 10:
            chunksize = 1
        else:
            chunksize = 10

    if star:
        f_starred = functools.partial(expand_star_args, f=f)
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


def run_parallel(f, f_args, nproc, f_const=None, chunksize=1, ordered=True, star=False, return_gen=False,
                 processes=False):
    """
    :param function f: function taking single argument or multiple arguments if using star=True
    :param Sequence f_args: sequence of arguments to pass to f
    :param int nproc: the number of threads/processes to use
    :param f_const: An object shared across workers accessible in the function as f.constants
    :param int|str chunksize: positive integer or 'auto'
    :param bool ordered: Whether to preserve the order of outputs
    :param bool star: If set to True, each of f_args is a Sequence of positional arguments to pass to f
    :param bool return_gen: Whether to return a generator instead of list
    :param bool processes: Whether to spawn processes instead of threads (when the f is cpu intensive / hasa lot of
     python code)
    :return: list or generator of return values from f
    :rtype: list | Generator
    """
    # if processes and star:
    #     raise Exception('The arguments processes and star are incompatible (for now)')

    r = run_parallel_gen(f, f_args, nproc, f_const, chunksize, ordered, star, processes)
    if return_gen:
        return r
    else:
        return list(r)


def mash_sketch(genome_names, genome_extension, genomes_dir, sketch_dir, nproc_cpu, nproc_io, sketch_size=10000):
    threads_per_run = get_threads_per_run(nproc_cpu, nproc_io)

    if genome_extension.endswith('.bz2'):
        cat_f = 'bzcat'
    elif genome_extension.endswith('.gz'):
        cat_f = 'zcat'
    else:
        cat_f = 'cat'

    commands = []
    reused = []
    sketch_files = []
    for g in genome_names:
        genome_file = genomes_dir / f'{g}{genome_extension}'
        sketch_file = sketch_dir / f'{g}.msh'
        sketch_files.append(sketch_file)
        if sketch_file.exists():
            reused.append(sketch_file)
            continue

        assert genome_file.exists()

        commands.append((f"{cat_f} {genome_file} | mash sketch -p {threads_per_run} -k 21 -s {sketch_size} "
                        f"-o {{}} -I {g} -", sketch_file))

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(run_command_with_output_file, commands, nproc_io, chunksize='auto', ordered=False, star=True)
    return sketch_files


def mash_sketch_aa(genome_to_faa, sketch_dir, nproc_cpu, nproc_io, sketch_size=10000):
    threads_per_run = get_threads_per_run(nproc_cpu, nproc_io)

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
                         f"-o {{}} -I {g} -", sketch_file))

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(run_command_with_output_file, commands, nproc_io, chunksize='auto', ordered=False, star=True)
    return sketch_files


def skani_sketch(genome_ids, genome_extension, genomes_dir, sketch_dir, nproc_cpu, nproc_io):
    threads_per_run = get_threads_per_run(nproc_cpu, nproc_io)

    def skani_sketch_one(sketch_dir_, genome_file_, sketch_file_, g_):
        tmp_dir_ = sketch_dir_ / g_
        if tmp_dir_.exists():
            shutil.rmtree(tmp_dir_)
        run_command(f"skani sketch -t {threads_per_run} -o {tmp_dir_} {genome_file_}")
        os.rename(tmp_dir_ / f'{g_}{genome_extension}.sketch', sketch_file_)
        os.unlink(tmp_dir_ / "markers.bin")
        os.rmdir(tmp_dir_)


    f_args = []
    reused = []
    sketch_files = []
    for g in genome_ids:
        genome_file = genomes_dir / f'{g}{genome_extension}'
        sketch_file = sketch_dir / f'{g}.sketch'
        sketch_files.append(sketch_file)
        if sketch_file.exists():
            reused.append(sketch_file)
            continue

        assert genome_file.exists()

        if genome_extension.endswith('.gz') or genome_extension.endswith('.bz2'):
            raise NotImplemented('We need to implement decompression for skani sketch')


        f_args.append([sketch_dir, genome_file, sketch_file, g])

    if reused:
        info(f'  Reused {len(reused)} existing sketch files')

    run_parallel(skani_sketch_one, f_args, nproc_io, star=True, ordered=False)
    return sketch_files


def chunks(seq, chunk_size):
    for i in range(0, len(seq), chunk_size):
        yield seq[i:i + chunk_size]


def mash_paste_prepare(genome_sketches, pastes_dir, chunk_size, file_prefix):
    genome_sketches = sorted(genome_sketches)  # for reproducibility
    paste_files = []
    commands_with_output = []
    reused = []
    for k, sketches_chunk in enumerate(chunks(genome_sketches, chunk_size)):
        sketches_txt = pastes_dir / f'{file_prefix}{k}.txt'
        paste_prefix = pastes_dir / f'{file_prefix}{k}'
        with open(sketches_txt, 'w') as f:
            f.write('\n'.join(map(str, sketches_chunk)) + '\n')
        paste_file = pastes_dir / f'{file_prefix}{k}.msh'
        paste_files.append(paste_file)
        if paste_file.exists():
            reused.append(paste_file)
            continue

        commands_with_output.append((f"mash paste -l {paste_prefix} {sketches_txt}", paste_file))

    return commands_with_output, paste_files, reused



def mash_paste(genome_sketches, pastes_dir, nproc, chunk_size, file_prefix='paste_pt'):
    """
    Paste MASH sketches of genomes.

    :param Iterable[pathlib.Path|str] genome_sketches:
    :param pathlib.Path pastes_dir:
    :param int nproc:
    :param int chunk_size:
    :param str file_prefix:
    :return: List of mash pastes created
    """

    commands, paste_files, reused = mash_paste_prepare(genome_sketches, pastes_dir, chunk_size, file_prefix)

    if reused:
        info(f'  Reused {len(reused)} existing paste files')

    run_parallel(run_command, commands, nproc, ordered=False)
    return paste_files


def skani_paste(sketches, pastes_dir, chunk_size, file_prefix='paste'):
    """
    Analog to MASH paste, but with skani we just write a txt file as a list of per-genome sketches

    :param Iterable[pathlib.Path|str] sketches:
    :param pathlib.Path pastes_dir:
    :param int chunk_size:
    :param str file_prefix:
    :return:
    """
    paste_files = []
    for i, sketches_chunk in enumerate(chunks(sketches, chunk_size)):
        paste_file = pastes_dir / f'{file_prefix}_pt{i}.txt'
        with open(paste_file, 'w') as f:
            f.write('\n'.join(map(str, sketches_chunk)))

        paste_files.append(paste_file)

    return paste_files


def mash_dist_block(mash_pastes_1, mash_pastes_2, dists_dir, nproc, progress_bar=True):
    """
    Runs a big rectangle of MASH distances

    :param Iterable[pathlib.Path] mash_pastes_1:
    :param Sequence[pathlib.Path] mash_pastes_2:
    :param pathlib.Path dists_dir:
    :param int nproc:
    :param bool progress_bar:
    :return The created distance files
    """

    if progress_bar:
        mash_pastes_1 = tqdm(mash_pastes_1)

    dist_files = []
    for mp1 in mash_pastes_1:
        dist_filename = mp1.stem + '.tsv.gz'
        dist_file = dists_dir / dist_filename
        cmd = f"mash dist -p {nproc} -t {mp1} {' '.join(map(str, mash_pastes_2))} | gzip > {{}}"
        run_command_with_output_file(cmd, dist_file, shell=True)
        dist_files.append(dist_file)

    return dist_files


def skani_dist_block(skani_pastes_1, skani_pastes_2, dists_dir, nproc, progress_bar=True):
    """
    Runs a big rectangle of skani distances

    :param Iterable[pathlib.Path] skani_pastes_1:
    :param Sequence[pathlib.Path] skani_pastes_2:
    :param pathlib.Path dists_dir:
    :param int nproc:
    :param bool progress_bar:
    :return The created distance files
    """

    if progress_bar:
        skani_pastes_1 = tqdm(skani_pastes_1)

    dist_files_matrix = []
    for p1 in skani_pastes_1:
        dist_files_row = []
        for p2 in skani_pastes_2:
            dist_filename = p1.stem + '__' + p2.stem + '.tsv.gz'
            dist_file = dists_dir / dist_filename
            cmd = f"skani dist -t {nproc} --ql {p1} --rl {p2} | gzip > {{}}"
            run_command_with_output_file(cmd, dist_file, shell=True)
            dist_files_row.append(dist_file)
        dist_files_matrix.append(dist_files_row)

    return dist_files_matrix


def skani_triangle_big(skani_pastes, dists_dir, nproc):
    """
    Runs a big triangle of skani distances

    :param Sequence[pathlib.Path] skani_pastes:
    :param pathlib.Path dists_dir:
    :param int nproc:
    :return:
    """
    for mp in skani_pastes:  # diagonal triangles
        dist_file = dists_dir / (mp.stem + '.tri.gz')
        cmd = f'skani triangle -t {nproc} -l {mp} | gzip > {{}}'
        run_command_with_output_file(cmd, dist_file)

    for mp1, mp2 in it.combinations(skani_pastes, 2):  # pairwise blocks
        dist_file = dists_dir / (mp1.stem + '__' + mp2.stem + '.tsv.gz')
        cmd = f'skani dist -t {nproc} --ql {mp1} --rl {mp2} | gzip > {{}}'
        run_command_with_output_file(cmd, dist_file)


def load_pwd_pandas(file_path, **kwargs):
    return pd.read_csv(file_path, sep='\t', index_col=0, header=0, **kwargs)


def load_skani_as_pwd(path, fix_queries=None, fix_refs=None, queries=None, refs=None, dtype=np.float32):
    """
    Loads skani dist output as distance matrix. You should pass the queries and references you used for the distancing
    as the skani output contains only hits above certain allelic fraction threshold, otherwise they might be missing in
    the output matrix.

    :param pathlib.Path|str path:
    :param fix_queries:
    :param fix_refs:
    :param Sequence[str] queries:
    :param Sequence[str] refs:
    :param dtype:
    :return: two DataFrames with references as index and queries as columns, one for ANI one for ref. AF
    """

    column_names = ['ref', 'query', 'ani', 'af_ref']
    rows = []
    with openr(path) as f:
        next(f)  # header row
        for line in f:
            r, q, ani, af_r, _ = line.split('\t', maxsplit=4)
            ani = dtype(ani)
            rows.append([r, q, ani, af_r])
    df = pd.DataFrame(rows, columns=column_names)
    df_ani = df.pivot(index='ref', columns='query', values='ani').astype(dtype)
    df_afr = df.pivot(index='ref', columns='query', values='af_ref').astype(dtype)

    result = []
    for df_a in [df_ani, df_afr]:
        if fix_refs is not None:
            df_a.index = df_a.index.map(fix_refs)
        if fix_queries is not None:
            df_a.columns = df_a.columns.map(fix_queries)
        df_a = df_a.reindex(index=refs, columns=queries).fillna(0).astype(dtype)
        result.append(df_a)

    return result


def load_big_triangle_skani(skani_pastes, dists_dir, genome_extension, dtype=np.float32):
    """
    Loads a big skani triangle into a pandas DataFrame. It takes the paste files and assumes you ran the big triangle
    distancing

    :param Iterable[pathlib.Path] skani_pastes:
    :param pathlib.Path dists_dir:
    :param str genome_extension:
    :param dtype:
    :return:
    """
    skani_pastes = list(skani_pastes)  # fix the ordering

    paste_indexes_orig = []
    for paste in skani_pastes:
        assert str(paste).endswith('.txt')
        with open(paste) as f:
            paste_index = f.read().strip().split('\n')
            paste_index = list(map(fix_skani_id('.sketch'), paste_index))
        paste_indexes_orig.append(paste_index)

    n = sum(map(len, paste_indexes_orig))

    matrix = np.zeros(shape=(n, n), dtype=dtype)

    # Load the diagonal
    index_all = []
    offset = 0
    for paste, paste_index_orig in zip(skani_pastes, paste_indexes_orig):
        dist_file = dists_dir / (paste.stem + '.tri.gz')

        paste_index = []
        with openr(dist_file) as f:
            first_line = next(f)
            n_lines = int(first_line.strip())
            assert n_lines == len(paste_index_orig)

            for i, line in enumerate(f):
                if i == 0:  # first line does not have any numbers in triangle
                    paste_index.append(fix_skani_id(genome_extension)(line.rstrip('\n')))
                    continue

                ix = line.index('\t')
                idx = fix_skani_id(genome_extension)(line[:ix])
                paste_index.append(idx)

                values = np.fromstring(line[ix + 1:], sep='\t', dtype=dtype)

                assert len(values) == i
                matrix[offset + i, offset:offset + i] = values
                matrix[offset:offset + i, offset + i] = values

            assert len(paste_index) == n_lines
            assert sorted(paste_index) == sorted(paste_index_orig)
            offset += n_lines

        index_all.extend(paste_index)

    assert len(index_all) == len(set(index_all))
    assert len(index_all) == n
    idx_to_pos = dict(zip(index_all, range(len(index_all))))

    # load the blocks
    for mp1, mp2 in it.combinations(skani_pastes, 2):  # pairwise blocks
        dist_file = dists_dir / (mp1.stem + '__' + mp2.stem + '.tsv.gz')

        with openr(dist_file) as f:
            next(f)
            for line in f:
                r, q, a, _ = line.split('\t', maxsplit=3)
                r = fix_skani_id(genome_extension)(r)
                q = fix_skani_id(genome_extension)(q)
                a = dtype(a)
                i = idx_to_pos[r]
                j = idx_to_pos[q]
                assert matrix[i, j] == 0
                assert matrix[j, i] == 0
                matrix[i, j] = a
                matrix[j, i] = a

    np.fill_diagonal(matrix, 100.0)
    df = pd.DataFrame(matrix, index=index_all, columns=index_all)

    return df


def load_mash_dist_block(dist_files, n_rows, n_columns, progress_bar=False, dtype=np.float32):
    """
    Loads a big MASH distance rectangle into a pandas DataFrame. You need to remember the size of the block from when
     you were distancing.

    :param Iterable[pathlib.Path|str] dist_files:
    :param int n_rows:
    :param int n_columns:
    :param bool progress_bar:
    :param dtype:
    :return:
    """
    matrix = np.full((n_rows, n_columns), fill_value=-1, dtype=dtype)

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
                values = np.fromstring(line[i + 1:], sep='\t', dtype=dtype)
                assert len(values) == n_columns_current
                matrix[row_i, column_offset:(column_offset+n_columns_current)] = values
            current_index = pd.Index(current_index)
            if index is None:
                index = current_index
            assert (index == current_index).all()
            column_offset += n_columns_current

    df = pd.DataFrame(matrix, index=index, columns=header)
    assert not (df < 0).any().any()  # all values were filled
    return df


def cluster_skani_pwd(pwd, t_ani):
    """

    :param pd.DataFrame pwd:
    :param float t_ani:
    :return: pd.Series mapping genomes (index of pwd) to the cluster numbers
    """
    assert (pwd.index == pwd.columns).all(), "The matrix is not a pair-wise matrix"

    genomes_to_be_clustered = pwd.index
    condensed_dist_mat = spd.squareform(100 - pwd.values)

    if len(genomes_to_be_clustered) > 1:
        dendro = spch.linkage(condensed_dist_mat, method='average')
        clusters = spch.fcluster(dendro, 100 - t_ani, criterion='distance')
    else:  # either one or zero genomes (the distance matrix is undefined)
        # assign all to a cluster "1"
        clusters = [1] * len(genomes_to_be_clustered)

    s_clusters = pd.Series(clusters, index=genomes_to_be_clustered)
    return s_clusters
