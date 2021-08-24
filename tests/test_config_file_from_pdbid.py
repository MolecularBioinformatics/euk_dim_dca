#!/usr/bin/env python3
"""
Tests for config_file_from_pdbid.py
"""
import sys
from pathlib import Path
import pytest

sys.path.append("../scripts")

from config_file_from_pdbid import *


def test_find_config_file_nofile():
    configpath = Path('../testdata')
    pdbid = '2345'
    correct = Path(f'../testdata/config_{pdbid}.txt')
    assert(find_config_file(configpath, pdbid)==correct)


def test_writeout_config_file_exists():
    configpath = Path('../testdata')
    pdbid = '1euc'
    with pytest.raises(FileExistsError):
        find_config_file(configpath, pdbid)
