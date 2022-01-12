#!/usr/bin/env python3
"""
edit_refseqs.py

Processes reference protein sequences for DCA workflow.

1. removes nonstandard amino acids (Xs) from the refseq
"""

from pathlib import Path

from io_utils import parse_fasta, writeout_fasta


def remove_nonstandard_aas(fastafile):
    """Removes Xs from sequence. Done
    on reference protein sequences to avoid
    trimming errors down the line. 

    :param fastafile: pathlib.PosixPath
    :returns fadict: dict of {header:seq}"""

    if not fastafile.is_file():
        raise FileNotFoundError(f'{fastafile} not found!')

    with open(fastafile, 'r') as fafile:
        fadict = parse_fasta(fafile.readlines())
    header = next(iter(fadict))
    seq = fadict[header]

    if 'X' in seq:
        print(f'{fastafile.stem} contains nonstandard residues.')
        fadict[header] = ''.join(seq.split('X'))
    else:
        print(f'No nonstandard residues found in {fastafile.stem}.')
    return fadict


def edit_refseqs(refseqfile1, refseqfile2, redo):
    """Processes refseq protein sequence. Overwrites 
    original fasta file. 

    :param refseqfile: pathlib.PosixPath"""

    edited_refseq1 = remove_nonstandard_aas(refseqfile1)
    edited_refseq2 = remove_nonstandard_aas(refseqfile2)

    writeout_fasta(refseqfile1, edited_refseq1, overwrite=True)
    writeout_fasta(refseqfile2, edited_refseq2, overwrite=True)

    return refseqfile1, refseqfile2
