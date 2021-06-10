#!/usr/bin/env python3
"""
process_phmmerhits.py

Takes in a phmmer-parsed-out keyfile.

Checks if more than 100 seq hits are found,
puts a cap of total number of sequences at 600,
and matches organisms.
"""

from pathlib import Path

from io_utils import get_globbed_list, does_target_exist
from io_utils import readin_list, writeout_list


def check_two_keyfiles(pathtophmmer, pdbid):
    """Checks dir for exactly 2 phmmer keyfiles for pdbid"""

    keyfilefmt = f'{pdbid}*_refseq_phmmer.keyfile' 
    keyfiles = get_globbed_list(pathtophmmer, keyfilefmt)

    if len(keyfiles) != 2:
        raise Exception(f'Number of keyfiles {len(keyfiles)} not equal to 2:\n{keyfiles}')
    else:
        return keyfiles


def check_hit_nr(pathtokeyfile, minhitnum=100):
    """Checks that num of hits in keyfile is higher
    than minhitnum"""

    hits = readin_list(pathtokeyfile)
    
    if len(hits) <= minhitnum:
        raise Exception(f'Too few hits returned from phmmer: {len(hits)}')
    else:
        return hits 
