#!/usr/bin/env python3
"""
prepare_seqs.py

Prepared easel-extracted seqs to be 
aligned. Does this by checking each fasta
file with hits for presence/absence of original
reference sequence, and appending/adjusting accordingly.
"""

from pathlib import Path

from io_utils import parse_fasta, fa_todict, writeout_fasta

def check_if_seq_in_dict(querydict, fastadict):
    """Checks if query seq is in fasta dict,
       returns key if exist, otherwise returns False,
       needs an absolute match.

    :param querydict: dict of refseq
    :param fastadict: dict of all hit seqs"""
    qseq = list(querydict.values())[0]

    for key, seq in fastadict.items():
        if qseq == seq:
            return key
            break
    return False


def remove_orgseq(somedict, orgtag):
    """Removes entry corresponding to certain
       orgtag from somedict.

    :param somedict: dict of seqs
    :param orgtag: 5-letter uniprot organism id"""
    for key in somedict.keys():
        if key.split()[0][-6:] == orgtag:
            print(f'\n     Removing {key}:')
            print(f'{somedict[key]}\n')
            del somedict[key]
            break
        else:
            continue
    return somedict


def rename_seqorg(queryfile, querydict, orgtag=None):
    """Relabels header of seq with RFSEQ tag.

    :param queryfile: path to refseqfile
    :param querydict: dict with hit seqs
    :opt param orgtag: 5-letter uniprot organism id"""
    stem = queryfile.stem
    oldhead = list(querydict.keys())[0]
    if orgtag:
        neworgtag = f'pdb|{stem}|{orgtag} {oldhead}'
    elif not orgtag:
        neworgtag = f'pdb|{stem}|_RFSEQ {oldhead}'
    newdict = {neworgtag:list(querydict.values())[0]}
    return newdict


def prepare_seqs(fastafile1, fastafile2, refseq1, refseq2):
    """Ensures refseqs are in respective fastafiles.
    Relabels refseqs as 'RFSEQ'."""

    fadict1 = fa_todict(fastafile1)
    fadict2 = fa_todict(fastafile2)

    refseqdict1 = fa_todict(refseq1)
    refseqdict2 = fa_todict(refseq2)

    key1 = check_if_seq_in_dict(refseqdict1, fadict1)
    key2 = check_if_seq_in_dict(refseqdict2, fadict2)
    
    if key1 and key2:
        orgtag1 = key1.split()[0][-6:] 
        orgtag2 = key2.split()[0][-6:]
        newrefdict1 = rename_seqorg(refseq1, refseqdict1)
        newrefdict2 = rename_seqorg(refseq2, refseqdict2)
        if not orgtag1 == orgtag2:
            print(f'Different organisms:\n{orgtag1}:{refseq1}\n{orgtag2}:{refseq2}')
            fadict1 = remove_orgseq(fadict1, orgtag1)
            fadict1 = remove_orgseq(fadict1, orgtag2)
            fadict2 = remove_orgseq(fadict2, orgtag1)
            fadict2 = remove_orgseq(fadict2, orgtag2)
        else:
            fadict1 = remove_orgseq(fadict1, orgtag1)
            fadict2 = remove_orgseq(fadict2, orgtag1)
        fadict1[newrefdict1.keys()[0]]=newrefdict1[newrefdict1.keys()[0]]
        fadict2[newrefdict2.keys()[0]]=newrefdict2[newrefdict2.keys()[0]]
#    elif key1 and not key2:
#        continue
#    elif key2 and not key1:
#        continue
#    else: 
#        continue
