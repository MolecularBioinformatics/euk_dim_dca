#!/usr/bin/env python3
"""
Tests for run_easel_getseqs.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from run_easel_getseqs import *

def test_run_easel_no_files():
    easelpath = Path('/home/yin/UiT/Tools/easel/miniapps')
    dbpath = Path()
    phmmerpath = Path('../testdata')
    keyfilepath = Path('../testdata/2222_A_refseq_phmmer.keyfile')
    with pytest.raises(FileNotFoundError):
        run_easel(easelpath, dbpath, phmmerpath, keyfilepath, False)

def test_run_easel_file_already_exists_no_redo():
    easelpath = Path('/home/yin/UiT/Tools/easel/miniapps')
    dbpath = Path()
    phmmerpath = Path('../testdata')
    keyfilepath = Path('../testdata/1c0f_A_refseq_phmmer_matched.keyfile')
    res = run_easel(easelpath, dbpath, phmmerpath, keyfilepath, False)
    assert(res == Path('../testdata/1c0f_A_refseq_phmmer_matched.fasta'))

def test_run_easel_file_easel_unsuccessful():
    easelpath = Path('/home/yin/UiT/Tools/easel/miniapps')
    dbpath = Path()
    phmmerpath = Path('../testdata')
    keyfilepath = Path('../testdata/1c0f_A_refseq_phmmer_matched.keyfile')
    with pytest.raises(ValueError):
        run_easel(easelpath, dbpath, phmmerpath, keyfilepath, True)
