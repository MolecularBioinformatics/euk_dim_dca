#!/usr/bin/env python3
"""
process_alnseqs.py

Matches aligned sequences by organism,
joins the matched sequences.

Prepares an alignment for DCA.
"""

from pathlib import Path

from io_utils import does_target_exist, parse_fasta, fa_todict, writeout_fasta


def get_orgdict_from_fafile(fafilepath):
    """Takes in fasta file (with uniprot-style headers),
    returns a dict of format {'ORG':('header', 'seq'),...}

    :param fafilepath: pathlib.PosixPath
    :returns fa_orgdict: dict"""

    fadict = fa_todict(fafilepath)
    fa_orgdict = {}
    for header in fadict.keys():
        firstpart = header.strip().split()[0]
        orgtag = firstpart.split("_")[-1:][0]
        fa_orgdict[orgtag] = (header, fadict[header])
    return fa_orgdict


def join_two_orgdicts(orgdict1, orgdict2):
    """Takes in two dicts of format {'ORG':('header', 'seq'),...},
    Joins the headers with '|' in between, and
    joins the sequences. Makes a joint file for DCA.

    :param orgdict: dict
    :param dcadict: dict, joined headers and aligned seqs"""

    dcadict = {}
    
    if len(orgdict1) != len(orgdict2):
        raise ValueError('Fastas you are trying to join have unequal # seqs!')
    for orgtag in orgdict1.keys():
        header1 = orgdict1[orgtag][0]
        header2 = orgdict2[orgtag][0]
        seq1 = orgdict1[orgtag][1]
        seq2 = orgdict2[orgtag][1]
        dcadict[header1+'||'+header2] = seq1+seq2
    return dcadict


def process_alnseqs(alnfile1_path, alnfile2_path, phmmerpath, redo):
    """Prepares aligned sequences from two alignment files for 
    DCA. First matches the sequences based on organism, then 
    zips them up (joins) them into one joint alignment file.

    :param alnfile_path: pathlib.PosixPath, 1 or 2 being either of the fastas
    :param phmmerpath: pathlib.PosixPath
    :param redo: bool

    :returns jointalnfile_path: pathlib.PosixPath"""

    ch1 = alnfile1_path.stem.split('_refseq')[0]
    ch2 = alnfile2_path.stem.split('_refseq')[0]
    outpath = phmmerpath / f"Joint_{ch1}_{ch2}_aln.fasta"

    if not (does_target_exist(alnfile1_path, 'file') and does_target_exist(alnfile2_path, 'file')):
        raise FileNotFoundError('Check that your alignment files exist!')
    if redo == False and does_target_exist(outpath, 'file'):
        print(f'Joint alignment already exists: {outpath}')
        return outpath

    orgdict1 = get_orgdict_from_fafile(alnfile1_path)
    orgdict2 = get_orgdict_from_fafile(alnfile2_path)

    jointdict = join_two_orgdicts(orgdict1, orgdict2)

    writeout_fasta(outpath, jointdict, overwrite=True)
    print(f'Joint alignment written into: {outpath}')

    return outpath
