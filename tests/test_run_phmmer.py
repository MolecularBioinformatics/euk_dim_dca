#!/usr/bin/env python3
"""
Tests for run_phmmer.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from run_phmmer import *


def test_phmmerlog_formatter():
    seqpath = Path('../testdata/ZZZZ_A_refseq.fasta')
    correctfmt = 'ZZZZ_A_refseq_phmmer.log'
    res = phmmerlog_formatter(seqpath)
    assert(res == correctfmt)

def test_run_phmmer_filenotfound():
    dbpath = Path()
    seqpath = Path('../testdata/ZZZZ_A_refseq.fasta')
    phmmerpath = Path('../testdata')
    redo = False
    with pytest.raises(FileNotFoundError):
        run_phmmer(dbpath, seqpath, phmmerpath, redo)

def test_run_phmmer_run_unsuccessful():
    dbpath = Path()
    seqpath = Path('../testdata/1c0f_A_refseq.fasta')
    phmmerpath = Path('../testdata')
    redo = True
    with pytest.raises(ValueError):
        run_phmmer(dbpath, seqpath, phmmerpath, redo)

def test_run_phmmer_no_redo():
    dbpath = Path()
    seqpath = Path('../testdata/1c0f_A_refseq.fasta')
    phmmerpath = Path('../testdata')
    redo = False
    correctoutpath = Path('../testdata/1c0f_A_refseq_phmmer.log')
    res = run_phmmer(dbpath, seqpath, phmmerpath, redo)
    assert(res == correctoutpath)
