#!/usr/bin/env python3
"""
Tests for parse_accid_phmmerlog.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from parse_accid_phmmerlog import *


def test_get_accidlist():
    minilog = Path('../testdata/1d0d_A_refseq_phmmer.log')
    idlist = ['sp|P17726|TAP_ORNMO','tr|Q8I9U3|Q8I9U3_9ACAR','tr|O44121|O44121_9ACAR','tr|A0A293MND8|A0A293MND8_ORNER']
    phmmerpath = Path('../testdata')
    reslist = get_accidlist(phmmerpath / minilog)
    assert(reslist == idlist)

def test_parse_accid_phmmerlog_filenotfound():
    nonexist_log = Path('../testdata/ZZZZ_A_refseq_phmmer.log')
    outpath = Path('../testdata')
    with pytest.raises(FileNotFoundError):
        parse_accid_phmmerlog(nonexist_log, outpath, False)

def test_parse_accid_phmmerlog_fileempty():
    empty_log = Path('../testdata/2222_A_refseq_phmmer.log')
    outpath = Path('../testdata')
    with pytest.raises(ValueError):
        parse_accid_phmmerlog(empty_log, outpath, False)

def test_parse_accid_phmmerlog_nohits():
    nohits_log = Path('../testdata/2222_B_refseq_phmmer.log')
    outpath = Path('../testdata')
    with pytest.raises(ValueError):
        parse_accid_phmmerlog(nohits_log, outpath, False)

