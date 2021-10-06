#!/usr/bin/env python3
"""
run_dca.py

Runs dca for a given joint alignment

1. First trim MSA with pydca's trimmer
2. Then run plmdca, mfdca or gaussdca on aln"""

import sys
import time
import subprocess
from io_utils import does_target_exist
from meanfield_dca import meanfield_dca

sys.path.append("/cluster/projects/nn9795k/yin/pydca-master/pydca")


def run_pydca_mfdca(jointaln_path):
    """
    Spawns subprocess to run pydca mfdca.

    :param jointaln_path: pathlib.PosixPath
        - needs to be just the string due to pydca's code
    :param redo: bool

    :returns mfdca_FN_APC: list
    """

    mfdca_inst = meanfield_dca.MeanFieldDCA(str(jointaln_path), 'protein', pseudocount=0.5, seqid=0.8)

    start = time.perf_counter()
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    stop = time.perf_counter()
    print(f'mfDCA ran in {stop - start:0.4f} seconds')

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


def run_dca(jointaln_path, outpath, redo, method):
    """
    Runs dca method (default mf) on a joint alignment.

    Writes out scores to a []_scores.dat file.

    :param jointaln_path: pathlib.PosixPath
    :param outpath: pathlib.PosixPath
    :param redo: bool

    :returns scorefile_path: pathlib.PosixPath
    """

    outfilename = f'{jointaln_path.stem}_{method}dca_scores.dat'
    outfilepath = outpath / outfilename

    if not does_target_exist(jointaln_path, 'file'):
        raise FileNotFoundError(f'JOINT ALN FILE MISSING: Could not find {jointaln_path}')
    elif does_target_exist(outfilepath, 'file') and redo is False:
        print(f'DCA scores files: ({outfilepath}) already exists in {outfilepath.parent}')
        return outfilepath

    print(f"DCA approach: {method}")

    if method == "mf":
        # mean-fiel DCA approach by pydca
        dcascores = run_pydca_mfdca(jointaln_path)
        if not dcascores:
            raise ValueError('DCA run unsuccessful!')
        writeout_scores(dcascores, jointaln_path, outfilepath)

    elif method == "plm":
        # pseudo-likelihood maximization DCA approach by CCMPred
        pass  # TODO CCMPred
    else:
        # gaussian DCA approach by gaussDCA
        start = time.perf_counter()
        cmd = ["julia", "-t 8", "run_gaussdca.jl", str(jointaln_path), outfilepath]  # TODO threads?
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            raise ValueError(f'GaussDCA run unsuccessful!')
        stop = time.perf_counter()
        print(f'gaussDCA ran in {stop - start:0.4f} seconds')

        # re edit format of scoring file
        new_file_content = ""
        for line in open(outfilepath, "r"):
            new_file_content += line.replace(" ", "\t")

        # overwrite existing scoring file with tab-separated scores
        new_file = open(outfilepath, "w")
        new_file.write(new_file_content)
        new_file.close()

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
