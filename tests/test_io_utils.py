#!/usr/bin/env python3
"""
Tests for io_utils
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from io_utils import *

def test_refseq_formatter_typeerror():
    with pytest.raises(TypeError):
        refseq_formatter(123)

def test_refseq_formatter():
    assert(refseq_formatter('1c0f')=='1c0f*refseq.fasta')

def test_phmmerlog_formatter_typeerror():
    with pytest.raises(TypeError):
        phmmerlog_formatter('1c0f_A_refseq.fasta')

def test_phmmerlog_formatter():
    filepath=Path('../testdata/1c0f_A_refseq.fasta')
    assert(phmmerlog_formatter(filepath)=='1c0f_A_refseq_phmmer.log')

def test_keyfile_formatter_typeerror():
    with pytest.raises(TypeError):
        keyfile_formatter('1c0f_A_refseq_phmmer.log')

def test_keyfile_formatter():
    filepath=Path('../testdata/1c0f_A_refseq_phmmer.log')
    assert(keyfile_formatter(filepath)=='1c0f_A_refseq_phmmer.keyfile')

def test_matched_keyfile_formatter_typeerror():
    with pytest.raises(TypeError):
        matched_keyfile_formatter('1c0f_A_refseq_phmmer_matched.log')

def test_matched_keyfile_formatter():
    filepath=Path('../testdata/1c0f_A_refseq_phmmer.keyfile')
    assert(matched_keyfile_formatter(filepath)=='1c0f_A_refseq_phmmer_matched.keyfile')

def test_easeled_seq_formatter_typeerror():
    with pytest.raises(TypeError):
        easeled_seq_formatter('1c0f_A_refseq_phmmer_matched.keyfile')

def test_easeled_seq_formatter():
    filepath=Path('../testdata/1c0f_A_refseq_phmmer_matched.keyfile')
    assert(easeled_seq_formatter(filepath)=='1c0f_A_refseq_phmmer_matched.fasta')

def test_get_globbed_list():
    dirpath = Path('../testdata')
    tworefseqs = '1c0f*refseq.fasta'
    onerefseqs = '1d4x*refseq.fasta'
    nonrefseqs = '1euc*refseq.fasta'
    assert len(get_globbed_list(dirpath, tworefseqs)) == 2
    assert len(get_globbed_list(dirpath, onerefseqs)) == 1
    assert len(get_globbed_list(dirpath, nonrefseqs)) == 0

def test_does_target_exist_file():
    filepath1=Path('../testdata/1c0f_A_refseq_phmmer.log')
    filepath2=Path('../testdata/xxxx_A_refseq_phmmer.log')
    assert(does_target_exist(filepath1, 'file')==True)
    assert(does_target_exist(filepath2, 'file')==False)

def test_does_target_exist_dir():
    dirpath1=Path('../testdata')
    dirpath2=Path('../blabla')
    assert(does_target_exist(dirpath1, 'dir')==True)
    assert(does_target_exist(dirpath2, 'dir')==False)

def test_does_target_exist():
    filepath1=Path('../testdata/1c0f_A_refseq_phmmer.log')
    dirpath2=Path('../blabla')
    assert(does_target_exist(filepath1, 'else')==True)
    assert(does_target_exist(dirpath2, 'else')==False)
