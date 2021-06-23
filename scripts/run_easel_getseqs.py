#!/usr/bin/env python3
"""
run_easel_getseqs.py

Runs HMMER easel tool to extract sequences from database.

Takes in a matched keyfile.
"""

import subprocess
from pathlib import Path

from io_utils import does_target_exist

def easeled_seq_formatter(keyfilepath):
    """Returns format of FASTA file with seqs
    extracted from easel

    :param keyfilepath: pathlib.PosixPath to keyfile
    :returns: str, outfile name
    """
    return f'{keyfilepath.stem}.fasta'


def run_easel(easelpath, databasepath, phmmerpath, keyfilepath, redo):
    """
    Spawns subprocess to run esl-sfetch.

    :param easelpath:
    :param databasepath:
    :param keyfilepath:
    :param redo:
    """
    filename = easeled_seq_formatter(keyfilepath)
    outpath = phmmerpath.joinpath(filename)

    if not does_target_exist(keyfilepath, 'file'):
        raise FileNotFoundError(f'KEYFILE MISSING: Could not find {keyfilepath}!')
    elif does_target_exist(outpath, 'file') and redo == False:
        print(f'Easel-fetched Fasta file: ({outpath}) already exists in {outpath.parent}')
        return outpath
    else:
        try:
            f = open(outpath, 'w+')
            cmdargs = [f"{easelpath}/esl-sfetch",
                       "-f",
                       f"{databasepath}",
                       f"{keyfilepath}"] 
            proc = subprocess.run(cmdargs, stdout=f)
            f.close()
            if proc.returncode == 0:
                print(f'Retrieved seqs in fasta file: {outpath}')
            return outpath
        except:
            raise Exception('Could not extract seqs with easel')
