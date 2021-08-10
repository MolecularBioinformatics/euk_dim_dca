#!/usr/bin/env python3
"""
Tests for run_dca.py

Currently run into a segmentation error. Some issue with import
of the pydca module?
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from run_dca import *

def test_run_dca_fnotfound():
    jointalnfilepath = Path('../testdata/Joint_1d4x_A_1d4x_B_aln.fasta')
    with pytest.raises(FileNotFoundError):
        run_dca(jointalnfilepath, Path(), False)

def test_run_dca_scores_file_exists():
    dcascoresfilepath = Path('../testdata/Joint_4ged_B_4ged_A_aln_mfdca_scores.dat')
    res = run_dca(dcascoresfilepath, Path('../testdata'), False)
    assert(res == dcascoresfilepath)
