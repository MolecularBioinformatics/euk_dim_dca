#!/usr/bin/env python3
"""
Tests for process_alnseqs.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from process_alnseqs import *


def test_get_orgdict_from_fafile():
    fafilepath = Path('../testdata/1111_A_refseq_phmmer_matched.fasta')
    correctdict = {'HUMAN': ('tr|_HUMAN SOME OTHER USELESS STUFF', 'abc'),
                   'TOXCA': ('sp|_TOXCA BLA BLA', 'def'),
                   'RFSEQ': ('rfseq|1c0f_A_refseq|_RFSEQ XP_636088.1', 'ghi')}
    assert(get_orgdict_from_fafile(fafilepath)==correctdict)


def test_join_two_orgdicts():
    fafilepath1 = Path('../testdata/1111_A_refseq_phmmer_matched.fasta')
    orgdict1 = get_orgdict_from_fafile(fafilepath1)
    fafilepath2 = Path('../testdata/1111_S_refseq_phmmer_matched.fasta')
    orgdict2 = get_orgdict_from_fafile(fafilepath2)
    correctdict = {'tr|_HUMAN SOME OTHER USELESS STUFF||tr|_HUMAN SOME OTHER USELESS STUFF':'abcabc',
                   'sp|_TOXCA BLA BLA||sp|_TOXCA BLA BLA':'defdef',
                   'rfseq|1c0f_A_refseq|_RFSEQ XP_636088.1||rfseq|1c0f_A_refseq|_RFSEQ XP_636088.1':'ghighi'}
    assert(join_two_orgdicts(orgdict1, orgdict2)==correctdict)
    

def test_join_two_orgdicts_valerr():
    fafilepath1 = Path('../testdata/1111_A_refseq_phmmer_matched.fasta')
    orgdict1 = get_orgdict_from_fafile(fafilepath1)
    fafilepath2 = Path('../testdata/1c0f_A_refseq_phmmer_matched.fasta')
    orgdict2 = get_orgdict_from_fafile(fafilepath2)
    with pytest.raises(ValueError):
        join_two_orgdicts(orgdict1, orgdict2)


def test_process_alnseqs_fnotfound():
    alnfilepath1 = Path('../testdata/1c0f_A_refseq_phmmer_matched.aln')
    alnfilepath2 = Path('../testdata/1111_A_refseq_phmmer_matched.aln')
    with pytest.raises(FileNotFoundError):
        process_alnseqs(alnfilepath1, alnfilepath2, Path('../testdata'), False)

def test_process_alnseqs_noredo():
    alnfilepath1 = Path('../testdata/1c0f_A_refseq_phmmer_matched.aln')
    alnfilepath2 = Path('../testdata/1c0f_S_refseq_phmmer_matched.aln')
    res = Path('../testdata/Joint_1c0f_A_1c0f_S_aln.fasta')
    assert(process_alnseqs(alnfilepath1, alnfilepath2, Path('../testdata'), False)==res) 
