#!/usr/bin/env python3
"""
io_utils.py

Collection of file writing and checking utility functions.

Used with the eukaryotic dimer dca method
"""

import warnings
from pathlib import Path

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


def writeout_list(listofitems, outpath):
    """Writes out contents of a list, newline-separated
    to a file specified by outpath.

    :param listofitems: list to write out
    :param outpath: pathlib.PosixPath
    """
    with open(outpath, 'w+') as f:
        f.write('\n'.join(listofitems))

    print(f'{outpath} written.')
