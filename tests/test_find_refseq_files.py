#!/usr/bin/env python3
"""
Tests for find_refseq_files.py
"""

import pytest
import sys

sys.path.append("../scripts/")

from find_refseq_files import *
from io_utils import *

def test_iscorrect_pdbid_notstr():
    with pytest.raises(TypeError):
        iscorrect_pdbid(1111)

def test_iscorrect_pdbid_toolong():
    with pytest.raises(ValueError):
        iscorrect_pdbid('aaaaa')

def test_iscorrect_pdbid_lower():
    assert iscorrect_pdbid('1AUD') == '1aud'

def test_refseq_formatter():
    assert refseq_formatter('1aud') == '1aud*refseq.fasta'

def test_find_refseq_files_twofiles():
    dirpath = Path('../testdata')
    two = '1c0f'
    two_out = [Path('../testdata/1c0f_A_refseq.fasta'),
               Path('../testdata/1c0f_S_refseq.fasta')]
    assert find_refseq_files(two, dirpath) == two_out
     
def test_find_refseq_files_valueerror():
    dirpath = Path('../testdata')
    one = '1d4x'
    non = '1euc'
    with pytest.raises(ValueError):
        find_refseq_files(one, dirpath)
        find_refseq_files(non, dirpath)
        
def test_find_refseq_files_nofileserror():
    dirpath = Path('../testdata')
    one = '1111'
    two = '2222'
    with pytest.raises(FileNotFoundError):
        find_refseq_files(one, dirpath)
        find_refseq_files(two, dirpath)
