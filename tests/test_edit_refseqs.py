#!/usr/bin/env python3
"""
Tests for edit_refseqs.py
"""

import pytest
import sys

sys.path.append("../scripts/")

from edit_refseqs import *

def test_remove_nonstandard_aas_filenotfound():
    fapath = Path('../testdata/sample_notexists.fasta')
    with pytest.raises(FileNotFoundError):
        remove_nonstandard_aas(fapath)


def test_remove_nonstandard_aas():
    fapath = Path('../testdata/sample.fasta')
    res = remove_nonstandard_aas(fapath)
    assert res == {'seq1': 'ABDCEGBG'}

