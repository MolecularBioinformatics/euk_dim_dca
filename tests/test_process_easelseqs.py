#!/usr/bin/env python3
"""
Tests for process_easelseqs.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts/")

from process_easelseqs import *

def test_parse_easelerror_filenotfound():
    eslerror = Path('../testdata/1111.easelerror')
    with pytest.raises(FileNotFoundError):
        parse_easelerror(eslerror)

def test_parse_easelerror_emptyfile():
    eslerror = Path('../testdata/2222.easelerror')
    with pytest.raises(ValueError):
        parse_easelerror(eslerror)

def test_parse_easelerror():
    eslerror = Path('../testdata/1c0f.easelerror')
    res = {'9ACAR', '9AGAR', '9CAEN', 
            '9CETA', '9CNID', '9EUCA',
            '9HYPO', 'BRUMA', 'HYMNN',
            'ORNAN', 'RHIZD', 'SARHA',
            'STRPU', 'TARSY', 'TURTR'}
    assert(parse_easelerror(eslerror) == res)

def test_remove_seqs_of_org_filenotfound():
    fastafile = Path('../testdata/2222_B_refseq_phmmer_matched.fasta')
    dummyset = {'HUMAN', 'MOUSE'}
    with pytest.raises(FileNotFoundError):
        remove_seqs_of_org(fastafile, dummyset)

def test_remove_seqs_of_org():
    fastafile = Path('../testdata/1111_A_refseq_phmmer_matched.fasta')
    dummyset = {'HUMAN', 'MOUSE'}
    resdict = {'sp|_TOXCA BLA BLA':'ddd'}
    assert(remove_seqs_of_org(fastafile, dummyset)==resdict)
