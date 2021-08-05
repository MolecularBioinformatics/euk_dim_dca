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

def run_easel(easelpath, databasepath, fastapath, keyfilepath, redo):
    """
    Spawns subprocess to run esl-sfetch.

    :param easelpath: pathlib.PosixPath
    :param databasepath: pathlib.PosixPath
    :param keyfilepath: pathlib.PosixPath
    :param redo: bool, whether to re-extract

    :returns outpath: pathlib.PosixPath, path to fasta with extracted seqs
    """
    filename = easeled_seq_formatter(keyfilepath)
    outpath = fastapath.joinpath(filename)

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


def writeout_seqsnotfound(listoferrmessages, keyfilepath):
    """Writes out easel error messages into a .easelerror file.
    This is done per set of keyfiles.

    :param listoferrmessages: list of stderrs, generated in run_easel_iterate 
    :param keyfilepath: pathlib.PosixPath, keyfile input in easel
    """
    # TODO: where should we output the errors
    errorfilepath = Path(f"{keyfilepath.stem.split('_')[0]}.easelerror")

    with open(errorfilepath, 'a') as errf:
        errf.write(f'======= {keyfilepath}\n')
        errf.write(''.join(listoferrmessages))


def run_easel_iterate(easelpath, databasepath, fastapath, keyfilepath, redo):
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
    outpath = fastapath.joinpath(filename)

    if not does_target_exist(keyfilepath, 'file'):
        raise FileNotFoundError(f'KEYFILE MISSING: Could not find {keyfilepath}!')
    elif does_target_exist(outpath, 'file') and redo == False:
        print(f'Easel-fetched Fasta file: ({outpath}) already exists in {outpath.parent}')
        return outpath

    with open(keyfilepath, 'r') as k:
        idlist=k.readlines()
        idlist=[item.strip() for item in idlist]

    seqsnotfound = []
    with open(outpath, 'w+') as f:
        for item in idlist:
            cmd = [f'{easelpath}/esl-sfetch',
                   f'{databasepath}',
                   f'{item}']
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if proc.returncode != 0:
                error = proc.stderr.decode("utf-8")
                print(error) 
                seqsnotfound.append(error)
            seq = proc.stdout.decode("utf-8")
            f.write(seq)

    if seqsnotfound:
        writeout_seqsnotfound(seqsnotfound, keyfilepath)

    print(f'Fasta file with sequences written: {outpath}')

    return outpath
