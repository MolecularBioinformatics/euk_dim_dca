#!/usr/bin/env python3
"""
Tests for process_phmmerhits.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts/")

from process_phmmerhits import get_two_keyfiles, get_hit_list, get_orgs_from_hitlist, select_seqheader_from_org
from ordered_set import OrderedSet

def test_get_two_keyfiles_globber():
    testphmmerpath = Path('../testdata')
    testpdbid = '1c0f'
    res = get_two_keyfiles(testphmmerpath, testpdbid)
    assert(res[0].name == '1c0f_A_refseq_phmmer.keyfile') 
    assert(res[1].name == '1c0f_S_refseq_phmmer.keyfile') 

def test_get_two_keyfiles_not_2_error():
    testphmmerpath = Path('../testdata')
    testpdbid = '1d4x'
    with pytest.raises(ValueError):
        get_two_keyfiles(testphmmerpath, testpdbid)


def test_get_hit_list_res_max():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    res = get_hit_list(testkeyfile, 10000, 'MAXIMUM')
    assert(len(res) == 6462)

def test_get_hit_list_res_min():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    res = get_hit_list(testkeyfile, 1000, 'MINIMUM')
    assert(len(res) == 6462)

def test_get_hit_list_error_max():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 1000, 'MAXIMUM')

def test_get_hit_list_error_min():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 10000, 'MINIMUM')

def test_get_hit_list_error_relation():
    testkeyfile = Path('../testdata/1c0f_S_refseq_phmmer.keyfile')
    with pytest.raises(ValueError):
        get_hit_list(testkeyfile, 1000, 'DEANRONALD')

def test_get_hit_list(rel='MINIMUM', thresh=0):
    minikeyfile=Path('../testdata/1c0f_A_refseq_phmmer.testkeyfile')
    res_hitlist = get_hit_list(minikeyfile, thresh, rel)
    correct_hitlist = ['tr|_HUMAN', 'sp|_TOXCA', 'tr|_TOXCA', 'sp|_PANTR']
    assert(res_hitlist == correct_hitlist)


def test_get_orgs_from_hitlist():
    hitlist = ['tr|_HUMAN', 'sp|_TOXCA', 'tr|_TOXCA', 'sp|_PANTR']
    res_orgset, res_orgdict = get_orgs_from_hitlist(hitlist)
    correct_orgset = OrderedSet(['HUMAN','TOXCA','PANTR'])
    correct_orgdict = {'HUMAN':['tr|_HUMAN'],
                       'TOXCA':['sp|_TOXCA','tr|_TOXCA'],
                       'PANTR':['sp|_PANTR']}
    assert(res_orgset == correct_orgset)
    assert(res_orgdict == correct_orgdict)


def test_select_seqheader_from_org():
    orgset = OrderedSet(['HUMAN','TOXCA','PANTR'])
    orgdict = {'HUMAN':['tr|_HUMAN'],
               'TOXCA':['sp|_TOXCA','tr|_TOXCA'],
               'PANTR':['sp|_PANTR']}
    res_headerset = select_seqheader_from_org(orgset, orgdict)
    correct_headerset = OrderedSet(['tr|_HUMAN', 'sp|_TOXCA', 'sp|_PANTR'])
    assert(res_headerset == correct_headerset)
