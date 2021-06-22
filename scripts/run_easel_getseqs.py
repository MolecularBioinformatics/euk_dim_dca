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

    f = open(outpath, 'w+')
    cmdargs = [f"{easelpath}/esl-sfetch",
               "-f",
               f"{databasepath}",
               f"{keyfilepath}"] 

    proc = subprocess.run(cmdargs, stdout=f)
    f.close()

    return proc

