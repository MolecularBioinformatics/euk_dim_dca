#!/usr/bin/env python3
"""
run_dca.py

Runs pydca for a given joint alignment

1. First trim MSA with pydca's trimmer
2. Then choose either plmdca or mfdca
"""

import time
import subprocess
from pathlib import Path

from io_utils import does_target_exist
from pydca.meanfield_dca import meanfield_dca
from pydca.sequence_backmapper import sequence_backmapper


def run_pydca_mfdca(jointaln_path, phmmerpath, redo):
    """
    Spawns subprocess to run pydca mfdca.

    :param jointaln_path: pathlib.PosixPath
    :param phmmerpath: pathlib.PosixPath
    :param redo: bool

    :returns mfdca_FN_APC: list
    """

    ## jointaln cannot be a path? just a string to the path?
    ## TODO: check this out in the source code
    mfdca_inst = meanfield_dca.MeanFieldDCA(jointaln_path,'protein', pseudocount = 0.5, seqid = 0.8)

    start = time.perf_counter()
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    stop = time.perf_counter()

    return mfdca_FN_APC


def writeout_scores(dcalist, jointalnpath, outpath, method='mfdca'):
    """
    Writes out dca scores in 3 column file.

    :param dcalist: list of tuples [((i,j),score),...]
    :param outpath: pathlib.Posix
    """

    outfilename = f'{jointalnpath.stem}_{method}_scores.dat'
    outfilepath = outpath / outfilename

    with open(outfilepath, 'w') as outf:
        for scorepair in dcalist:
            outf.write(f'{scorepair[0][0]}\t{scorepair[0][1]}\t{scorepair[1]}\n')

