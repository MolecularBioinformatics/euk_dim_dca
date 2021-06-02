#!/usr/bin/env python3
"""
find_refseq_files.py

Globs in a given directory for refseq files
corresponding to given 4-letter PDB ID
"""

from pathlib import Path
import pytest

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
        raise Exception('Give correct PDBID.')
    else:
        refmt_pdbid = refseq_formatter(fourletter)
        fileslist = get_globbed_list(dirpath, refmt_pdbid)
        if len(fileslist) != 2:
            raise Exception('Number refseq files not equal to 2!')
        else:
            return fileslist


def test_iscorrect_pdbid_notstr():
    with pytest.raises(TypeError):
        iscorrect_pdbid(1111)

def test_iscorrect_pdbid_toolong():
    with pytest.raises(Exception):
        iscorrect_pdbid('aaaaa')

def test_iscorrect_pdbid_lower():
    assert iscorrect_pdbid('1AUD') == '1aud'

def test_refseq_formatter():
    assert refseq_formatter('1aud') == '1aud*refseq.fasta'

def test_get_globbed_list():
    dirpath = Path('../testdata')
    tworefseqs = '1c0f*refseq.fasta'
    onerefseqs = '1d4x*refseq.fasta'
    nonrefseqs = '1euc*refseq.fasta'
    assert len(get_globbed_list(dirpath, tworefseqs)) == 2
    assert len(get_globbed_list(dirpath, onerefseqs)) == 1
    assert len(get_globbed_list(dirpath, nonrefseqs)) == 0

def test_find_refseq_files_twofiles():
    dirpath = Path('../testdata')
    two = '1c0f'
    two_out = [Path('../testdata/1c0f_A_refseq.fasta'),
               Path('../testdata/1c0f_S_refseq.fasta')]
    assert find_refseq_files(two, dirpath) == two_out
     
def test_find_refseq_files_error():
    dirpath = Path('../testdata')
    one = '1d4x'
    non = '1euc'
    with pytest.raises(Exception):
        find_refseq_files(one, dirpath)
        find_refseq_files(non, dirpath)
