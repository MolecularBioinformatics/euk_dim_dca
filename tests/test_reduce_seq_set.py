#!/usr/bin/env python3
"""
Tests for reduce_seq_set.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from reduce_seq_set import *

def test_find_long_seqs():
    testfilepath = Path('../testdata/testremove3.fasta')
    seqsdict, headlist = find_long_seqs(testfilepath, 50)
    assert(headlist[0].split('_')[-1:][0]=='MOUSE')

def test_get_orgset_to_delete():
    pass

def test_reduce_seq_set_orgset_empty():
    pass

def test_reduce_seq_set_orgset_total():
    pass

def test_orgsets_match():
    pass
