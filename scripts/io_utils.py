#!/usr/bin/env python3
"""
io_utils.py

Collection of file writing and checking utility functions.

Used with the eukaryotic dimer dca method
"""

import warnings
from pathlib import Path

def refseq_formatter(pdbid):
    """Returns formatted refseq fasta 
    str expression as a regex matching string"""
    return f'{pdbid}*refseq.fasta'


def keyfile_formatter(pathtophmmerlog):
    """Returns formatted keyfile from phmmerlogfile"""
    return f'{pathtophmmerlog.stem}.keyfile'


def get_globbed_list(pathtodir, target):
    """Searches directory for files matching
    a certain target pattern.

    :param pathtodir: pathlib.PosixPath 
    :param target: str with expression to match

    :returns: list of matching paths
    """
    p = pathtodir.expanduser()
    return list(p.glob(target))


def does_target_exist(pathtotarget):
    """Given directory and a target (file or dir),
    returns True if it exists, False if not.

    :param pathtotarget: pathlib.PosixPath

    :returns: bool
    """
    if pathtotarget.is_file() or pathtotarget.is_dir():
        return True
    else:
        return False


def readin_list(filepathoflist):
    """Reads in a file with \n spaced items
    into a list"""
    with open(filepathoflist, 'r') as r:
        lines = r.readlines()
    return [item.strip() for item in lines]


def writeout_list(listofitems, outpath):
    """Writes out contents of a list, newline-separated
    to a file specified by outpath.

    :param listofitems: list to write out
    :param outpath: pathlib.PosixPath
    """
    with open(outpath, 'w+') as f:
        f.write('\n'.join(listofitems))

    print(f'File: {outpath.name} written into dir: {outpath.parent}')
