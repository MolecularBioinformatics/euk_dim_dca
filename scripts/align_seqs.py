#!/usr/bin/env python3
"""
align_seqs.py

Runs muscle alignment for a given set of sequences.

Appends refseq back into seq file before alignment.
"""

import time
import subprocess
from pathlib import Path

from io_utils import does_target_exist, fa_todict, writeout_fasta


def add_refseq(fastafile_path, refseqfile_path):
    """Modifies fasta file with seqs by adding
    reference sequence. 

    :param fastafile_path: pathlib.PosixPath with seqs to align
    :param refseqfile_path: pathlib.PosixPath
    """
    fadict = fa_todict(fastafile_path)
    refseqdict = fa_todict(refseqfile_path)

    writeout_fasta(fastafile_path, fadict, overwrite=True, addict=refseqdict)


def run_muscle(fastafile_path, refseqfile_path, phmmerpath, redo):
    """
    Spawns subprocess to run muscle.

    :param fastafile_path: pathlib.PosixPath
    :param phmmerpath: pathlib.PosixPath

    :returns alnfile_path: pathlib.PosixPath, aligned fasta
    """
    outpath = phmmerpath / Path(f'{fastafile_path.stem}.aln')

    if not does_target_exist(fastafile_path, 'file'):
        raise FileNotFoundError(f'Fasta file with seqs to align not found: {fastafile_path}.')
    elif fastafile_path.stat().st_size == 0:
        raise ValueError(f'EMPTY FILE: {fastafile_path}.')
    elif does_target_exist(outpath, 'file') and redo==False:
        print(f'Alignment file {outpath} already exists. Give --redo True to realign.')

    add_refseq(fastafile_path, refseqfile_path)
    
    cmdargs = ["muscle", 
               "-in",
               f"{fastafile_path}",
               "-out",
               f"{outpath}"]

    start = time.perf_counter()
    proc = subprocess.run(cmdargs)
    stop = time.perf_counter()
    if proc.returncode != 0:
        raise ValueError(f'Muscle could not align {fastafile_path}.')
    print(f'Muscle ran in {stop-start:0.4f} seconds')
    print(f'Alignment stored in {outpath}')
    return outpath
