#!/usr/bin/env python3
"""
io_utils.py

Collection of file writing and checking utility functions.

Used with the eukaryotic dimer dca method
"""

from pathlib import Path

# search for matching files in dir

def get_globbed_list(pathtodir, target):
    """Searches directory for files matching
    a certain target pattern.

    :param pathtodir: pathlib.PosixPath 
    :param target: str with expression to match

    :returns: list of matching paths
    """
    p = pathtodir.expanduser()
    return list(p.glob(target))
