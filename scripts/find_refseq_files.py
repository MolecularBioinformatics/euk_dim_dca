#!/usr/bin/env python3
"""
find_refseq_files.py

Searches in a given directory for refseq files
corresponding to given 4-letter PDB ID
"""

from pathlib import Path

from io_utils import get_globbed_list

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
        raise Exception('Give correct PDBID.')
    else:
        refmt_pdbid = refseq_formatter(fourletter)
        fileslist = get_globbed_list(dirpath, refmt_pdbid)
        if len(fileslist) == 0:
            raise FileNotFoundError('No refseq files found.')
        elif len(fileslist) != 2:
            raise Exception('Number refseq files not equal to 2!')
        else:
            return fileslist

