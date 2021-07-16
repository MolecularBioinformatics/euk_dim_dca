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
    assert(phmmerlog_formatter(seqpath) == correctfmt)

def test_run_phmmer_filenotfound():
    pass
