#!/usr/bin/env python3
"""
run_easel_getseqs.py

Runs HMMER easel tool to extract sequences from database.

Takes in a matched keyfile.

Returns fasta file with extracted sequences.
"""

import subprocess
from pathlib import Path

from io_utils import does_target_exist, easeled_seq_formatter

def run_easel(easelpath, databasepath, phmmerpath, keyfilepath, redo):
    """
    Spawns subprocess to run esl-sfetch.

    :param easelpath: pathlib.PosixPath
    :param databasepath: pathlib.PosixPath
    :param keyfilepath: pathlib.PosixPath
    :param redo: bool, whether to re-extract

    :returns outpath: pathlib.PosixPath, path to fasta with extracted seqs
    """
    filename = easeled_seq_formatter(keyfilepath)
    outpath = phmmerpath.joinpath(filename)

    if not does_target_exist(keyfilepath, 'file'):
        raise FileNotFoundError(f'KEYFILE MISSING: Could not find {keyfilepath}!')
    elif does_target_exist(outpath, 'file') and redo == False:
        print(f'Easel-fetched Fasta file: ({outpath}) already exists in {outpath.parent}')
        return outpath

    f = open(outpath, 'w+')
    cmdargs = [f"{easelpath}/esl-sfetch",
                "-f",
               f"{databasepath}",
               f"{keyfilepath}"] 
    proc = subprocess.run(cmdargs, stdout=f)
    f.close()

    if proc.returncode != 0:
        raise ValueError('Easel extract unsuccessful!')

    print(f'Retrieved seqs in fasta file: {outpath}')
    return outpath


def run_easel_iterate(easelpath, databasepath, phmmerpath, keyfilepath, redo):
    """Easel stops if it cannot find a sequence, no way
    to get it to continue sequence extraction.
    
    This function iterates over all seq accession ids,
    extracting one sequence at a time

    :param easelpath: pathlib.PosixPath
    :param databasepath: pathlib.PosixPath
    :param keyfilepath: pathlib.PosixPath
    :param redo: bool, whether to re-extract
    """
    filename = easeled_seq_formatter(keyfilepath)
    outpath = phmmerpath.joinpath(filename)

    if not does_target_exist(keyfilepath, 'file'):
        raise FileNotFoundError(f'KEYFILE MISSING: Could not find {keyfilepath}!')
    elif does_target_exist(outpath, 'file') and redo == False:
        print(f'Easel-fetched Fasta file: ({outpath}) already exists in {outpath.parent}')
        return outpath

    # get list of ids

    for item in idlist:
        cmd = [f'{easelpath}/esl-sfetch',
               f'{databasepath}',
               f'{item}']
        proc = subprocess.run(cmd, stdout=subprocess.PIPE)
        if proc.returncode != 0:
            print('bla')  # issue some warning
        seq = proc.stdout.decode("utf-8")
