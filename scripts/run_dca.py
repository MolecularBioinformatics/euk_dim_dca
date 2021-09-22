#!/usr/bin/env python3
"""
run_dca.py

Runs pydca for a given joint alignment

1. First trim MSA with pydca's trimmer
2. Then choose either plmdca or mfdca
"""

import sys
sys.path.append("/cluster/projects/nn9795k/yin/pydca-master/pydca")

import time
import subprocess
from pathlib import Path

from io_utils import does_target_exist
from meanfield_dca import meanfield_dca
from sequence_backmapper import sequence_backmapper


def run_pydca_mfdca(jointaln_path, redo):
    """
    Spawns subprocess to run pydca mfdca.

    :param jointaln_path: pathlib.PosixPath
        - needs to be just the string due to pydca's code
    :param redo: bool

    :returns mfdca_FN_APC: list
    """

    mfdca_inst = meanfield_dca.MeanFieldDCA(str(jointaln_path),'protein', pseudocount = 0.5, seqid = 0.8)

    start = time.perf_counter()
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    stop = time.perf_counter()

    return mfdca_FN_APC


def writeout_scores(dcalist, jointalnpath, outfilepath, method='mfdca'):
    """
    Writes out dca scores in 3 column file.

    :param dcalist: list of tuples [((i,j),score),...]
    :param outpath: pathlib.PosixPath
    """

    with open(outfilepath, 'w') as outf:
        for scorepair in dcalist:
            outf.write(f'{scorepair[0][0]}\t{scorepair[0][1]}\t{scorepair[1]}\n')


def run_dca(jointaln_path, outpath, redo, method='mfdca'):
    """
    Runs dca method (default mfdca) on a joint alignment.

    Writes out scores to a []_scores.dat file.

    :param jointaln_path: pathlib.PosixPath
    :param outpath: pathlib.PosixPath
    :param redo: bool

    :returns scorefile_path: pathlib.PosixPath
    """

    outfilename = f'{jointaln_path.stem}_{method}_scores.dat'
    outfilepath = outpath / outfilename

    if not does_target_exist(jointaln_path, 'file'):
        raise FileNotFoundError(f'JOINT ALN FILE MISSING: Could not find {jointaln_path}')
    elif does_target_exist(outfilepath, 'file') and redo == False:
        print(f'DCA scores files: ({outfilepath}) already exists in {outfilepath.parent}')
        return outfilepath

    dcascores = run_pydca_mfdca(jointaln_path, redo)
    if not dcascores:
        raise ValueError('DCA run unsuccessful!')
    writeout_scores(dcascores, jointaln_path, outfilepath)
    print(f'DCA scores written into {outfilepath}')
    return outfilepath

def run_mfdca_saga(jointaln_path, outpath, redo):
    """
    Runs dca method on a joint alignment.
    Spawns subprocess to call method from the command line.
    Written to be run on saga (errors with msa_numerics and numba)

    :param jointaln_path: pathlib.PosixPath
    :param outpath: pathlib.PosixPath
    :param redo: bool

    :returns outfilepath: pathlib.PosixPath
    """

    outfilename = f'MFDCA_apc_fn_scores_{jointaln_path.stem}.txt'
    outfilepath = outpath / outfilename
    cmdargs = ['python3', 
               'mfdca_main.py',
               'compute_fn',
               'protein',
               f'{jointaln_path}',
               '--apc',
               '--pseudocount',
               '0.5', 
               '--verbose',
               '--output_dir',
               f'{outpath}']

    if not does_target_exist(jointaln_path, 'file'):
        raise FileNotFoundError(f'JOINT ALN FILE MISSING: Could not find {jointaln_path}')
    elif does_target_exist(outfilepath, 'file') and redo == False:
        print(f'DCA scores files: ({outfilepath}) already exists in {outfilepath.parent}')
        return outfilepath

    proc = subprocess.run(cmdargs)
    if proc.returncode != 0:
        raise ValueError(f'DCA run unsuccessful!')
    print(f'DCA scores written into {outfilepath}')
    return outfilepath
