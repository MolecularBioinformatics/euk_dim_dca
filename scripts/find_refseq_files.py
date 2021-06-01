#!/usr/bin/env python3
"""
find_refseq_files.py

Globs in a given directory for refseq files
corresponding to given 4-letter PDB ID
"""

from pathlib import Path

def iscorrect_pdbid(pdbid):
    """Checks if pdbid is a 4-char string
    Returns the lowercase pdbid"""
    if not isinstance(pdbid, str):
        raise TypeError('PDB ID must be str.')
    elif not len(pdbid)==4:
        raise Exception('PDB ID must have 4 characters.')
    else:
        return pdbid.lower()


def refseq_formatter(pdbid):
    """Returns formatted refseq fasta str expression

    :param pdbid: str
    
    :returns: str with matching regex
    """
    return f'{pdbid}*refseq.fasta'


def get_globbed_list(pathtodir, target):
    """Searches directory for matches to target

    :param pathtodir: pathlib.PosixPath 
    :param target: str with expression to match

    :returns: list of matching paths
    """
    p = pathtodir.expanduser()
    return list(p.glob(target))

def find_refseq_files(pdbid, dirpath):
    """Takes a pdbid and directory. 
       Searches directory for matches to the pdbid.
       Written to search for refseq files of form:
       PDBID_[chainID]_refseq.fasta

    :param pdbid: str
    :param pathtodir: pathlib.PosixPath
    :returns: list of pathlib.PosixPaths
    """
    try:
        fourletter=iscorrect_pdbid(pdbid)
    except:
        raise Exception('PDB ID incorrect!')

    if not fourletter:
        raise Exception('REDO')
    else:
        refmt_pdbid=refseq_formatter(fourletter)
        return get_globbed_list(dirpath, refmt_pdbid)
    
