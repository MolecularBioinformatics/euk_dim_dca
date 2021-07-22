#!/usr/bin/env python3
"""
process_easelseqs.py

Removes sequences from organisms that were not extractable
from the database during the easel step.

Takes in an easelerror file, two matched fasta files.

Returns two fasta files
"""

import subprocess
from pathlib import Path

from io_utils import does_target_exist, parse_fasta, fa_todict, writeout_fasta 


def parse_easelerror(easelerr_filepath):
    """Parses easel error file.
    Returns a list of orgs that are not extracted
    by easel from a database.
    
    :param easelerr_filepath: pathlib.PosixPath
    
    :returns: set, orgtags to be removed
    """

    if not does_target_exist(easelerr_filepath, 'file'):
        raise FileNotFoundError('Easel error file not found.')
    elif easelerr_filepath.stat().st_size == 0:
        raise ValueError(f'EMPTY FILE: {easelerr_filepath}.')

    orglist = []
    with open(easelerr_filepath, 'r') as e:
        for line in e.readlines():
            if line.startswith('seq'):
                orglist.append(line.split()[1])
    orglist = [item.split("_")[1] for item in orglist]
    return set(orglist)


def remove_seqs_of_org(fasta_filepath, setoforgs):
    """Removes seqs belonging to certain orgs
    from a fasta file. Returns modified fasta file.

    :param fasta_filepath: pathlib.PosixPath
    :param setoforgs: set, organisms for which to remove seqs
    """

    if not does_target_exist(fasta_filepath, 'file'):
        raise FileNotFoundError('Fasta file not found.')

    fadict = fa_todict(fasta_filepath)
    for key in list(fadict.keys()):
        uniprottag = key.strip().split()[0]
        orgtag = uniprottag.split("_")[1]
        if orgtag in setoforgs:
            del fadict[key]
    return fadict
