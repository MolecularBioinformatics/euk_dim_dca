#!/usr/bin/env python3
"""
run_phmmer.py

Runs phmmer for a given sequence and database

Requirements:
1. Phmmer installed and on system path
2. Input of seq, db, and path to output.
"""

import time
import subprocess
from pathlib import Path


def phmmerlog_formatter(seqpath):
    """
    Formats phmmer output files for given seqfile.
    
    :param seqpath: pathlib.PosixPath
    :returns: str, outfile name
    """
    return f'{seqpath.stem}_phmmer.log'


def run_phmmer(databasepath, seqpath, outpath):
    """
    Spawns subprocess to run phmmer.

    :param databasepath: pathlib.PosixPath, database to search
    --> should be given as path to fasta with database
    --> does this also require SSI index to be made in same folder?
    :param seqpath: pathlib.PosixPath, input seqfile
    :param outpath: pathlib.PosixPath, outpath logfile
    :returns: None
    """
    cmdargs = ["phmmer",
               "-o",
               f'{outpath}',
               "--noali",
               "--cpu",
               "4",
               f'{seqpath}',
               f'{databasepath}']

    start = time.perf_counter()
    proc = subprocess.run(cmdargs)
    stop = time.perf_counter()

    print(f'Phmmer ran in {stop-start:0.4f} seconds')

    return proc
