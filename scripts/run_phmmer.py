#!/usr/bin/env python3
"""
run_phmmer.py

Runs phmmer for a given sequence and database

Requirements:
1. Phmmer installed and on system path
2. Input of seq, db, and path to output.
3. redo boolean flag if phmmer needs to be rerun or not
"""

import time
import subprocess
from pathlib import Path

from io_utils import does_target_exist, phmmerlog_formatter


def run_phmmer(databasepath, seqpath, phmmerpath, redo):
    """
    Spawns subprocess to run phmmer.

    :param databasepath: pathlib.PosixPath, database to search
    --> should be given as path to fasta with database
    --> does this also require SSI index to be made in same folder?
    :param seqpath: pathlib.PosixPath, input seqfile
    :param phmmerpath: pathlib.PosixPath
    
    :returns: outpath or None
    """
    filename = phmmerlog_formatter(seqpath)
    outpath = phmmerpath.joinpath(filename)
    cmdargs = ["phmmer",
               "-o",
               f'{outpath}',
               "--noali",
               "--cpu",
               "6",
               f'{seqpath}',
               f'{databasepath}']

    if not does_target_exist(seqpath, 'file'):
        raise FileNotFoundError(f'REFSEQ FILE MISSING: Could not find {seqpath}!')
    elif does_target_exist(outpath, 'file') and redo == False:
        print(f'Phmmer logfile: ({outpath.name}) already exists in {outpath.parent}') 
        return outpath

    start = time.perf_counter()
    proc = subprocess.run(cmdargs)
    stop = time.perf_counter()
    if proc.returncode != 0:
        raise ValueError(f'Phmmer run unsuccessful for {seqpath}')
    print(f'Phmmer ran in {stop-start:0.4f} seconds')
    print(f'Phmmer log stored in {outpath}')
    return outpath
