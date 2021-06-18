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

from io_utils import does_target_exist

def phmmerlog_formatter(seqpath):
    """
    Formats phmmer output files for given seqfile.
    
    :param seqpath: pathlib.PosixPath
    :returns: str, outfile name
    """
    return f'{seqpath.stem}_phmmer.log'


def run_phmmer(databasepath, seqpath, phmmerpath):
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
               "4",
               f'{seqpath}',
               f'{databasepath}']

    if not does_target_exist(seqpath, 'file'):
        raise FileNotFoundError(f'Could not find {seqpath}!')
    else:
        if does_target_exist(outpath, 'file'):
            print(f'Phmmer logfile: ({outpath.name}) already exists in {outpath.parent}') 
            return outpath
        else:
            try:
                start = time.perf_counter()
                proc = subprocess.run(cmdargs)
                stop = time.perf_counter()
                if proc.returncode == 0:
                    print(f'Phmmer ran in {stop-start:0.4f} seconds')
                    print(f'Phmmer log stored in {outpath}')
                return outpath
            except:
                raise Exception(f'Phmmer run unsuccessful for {seqpath}')
