#!/usr/bin/env python3
"""
io_utils.py

Collection of file writing and checking utility functions.

Used with the eukaryotic dimer dca method
"""

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
