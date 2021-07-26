#!/usr/bin/env python3
"""
Tests for align_seqs.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from align_seqs import *

def test_run_muscle_fnotfound():
    fafilepath = Path('../testdata/2222_A_refseq_phmmer_matched.fasta')
    refseqpath = Path('../testdata/2222_A_refseq.fasta')
    phmmerpath = Path('../testdata/')
    with pytest.raises(FileNotFoundError):
        run_muscle(fafilepath, refseqpath, phmmerpath, False)

def test_run_muscle_emptyfile():
    fafilepath = Path('../testdata/3333_C_refseq_phmmer_matched.fasta')
    refseqpath = Path('../testdata/2222_A_refseq.fasta')
    phmmerpath = Path('../testdata/')
    with pytest.raises(ValueError):
        run_muscle(fafilepath, refseqpath, phmmerpath, False)

def test_run_muscle_aln_exists():
    fafilepath = Path('../testdata/1c0f_A_refseq_phmmer_matched.fasta')
    refseqpath = Path('../testdata/1c0f_A_refseq.fasta')
    phmmerpath = Path('../testdata/')
    res = Path('../testdata/1c0f_A_refseq_phmmer_matched.aln')
    assert(run_muscle(fafilepath, refseqpath, phmmerpath, False)==res)

