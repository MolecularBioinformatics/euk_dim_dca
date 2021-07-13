#!/usr/bin/env python3
"""
Tests for process_phmmerhits.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts/")

from process_phmmerhits import get_two_keyfiles, get_hit_list


def test_get_hit_list_res_max():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    res = get_hit_list(testkeyfile, 10000, 'MAXIMUM')
    assert(len(res) == 6462)

def test_get_hit_list_res_min():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    res = get_hit_list(testkeyfile, 1000, 'MINIMUM')
    assert(len(res) == 6462)

def test_get_hit_list_error_max():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 1000, 'MAXIMUM')

def test_get_hit_list_error_min():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 10000, 'MINIMUM')

def test_get_hit_list_error_relation():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 1000, 'DEANRONALD')

def test_get_two_keyfiles_globber():
    testphmmerpath = Path('../testdata')
    testpdbid = '1c0f'
    res = get_two_keyfiles(testphmmerpath, testpdbid)
    assert(res[0].name == '1c0f_A_refseq_phmmer.keyfile') 
    assert(res[1].name == '1c0f_S_refseq_phmmer.keyfile') 

def test_get_two_keyfiles_not_2_error():
    testphmmerpath = Path('../testdata')
    testpdbid = '1d4x'
    with pytest.raises(ValueError):
        get_two_keyfiles(testphmmerpath, testpdbid)
