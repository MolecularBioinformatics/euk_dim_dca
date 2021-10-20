#!/usr/bin/env python3
"""
reduce_seq_set.py

Reduces sequence set in two fasta files according
to a cutoff length.

Done to make alignments more manageable for Muscle,
which runs into memory issues if sequences are too long.

Current cutoff is set at 1600.
"""

from pathlib import Path
from io_utils import parse_fasta, fa_todict, writeout_fasta

def find_long_seqs(fastafile, cutoffval):
    """Parses a fastafile to return list
    of headers of seqs that are too long,
    and the original seqs dict.

    :param fastafile: pathlib.PosixPath
    :param cutoffval: int, max seq length

    :returns seqsdict: fasta dict {header:seq,}
    :returns headerlis: list of fasta headers
    """

    seqsdict = fa_todict(fastafile)

    headerlist = []
    for header in seqsdict.keys():
        if len(seqsdict[header]) >= cutoffval:
            headerlist.append(header)

    return seqsdict, headerlist


def get_orgset_to_delete(headerlist1, headerlist2):
    """Takes in two lists of fasta headers
    corresponding to seqs that are too long. 
    Returns a set of the organisms contained in 
    the two lists

    :param headerlist: list of fastaheaders
    :returns orgset: set of orgtags
    """
    tags = [header.split()[0].split('_')[-1:] for header in headerlist1+headerlist2]
    return set(tags)


def del_orgseq_from_dict(somedict, orgset):
    """Removes seqs from dict that are from
    organisms found in the orgset"""
    for key in list(somedict.keys()):
        if key.split()[0].split('_')[1] in orgset:
            del somedict[key]
    return somedict

def reduce_seq_set(fastafile1, fastafile2, cutoffval):
    """Takes in two fasta files and a max seqlength cutoff value.
    Removes sequences that are too long, preserves organism
    agreement between the two fasta files.

    :param fastafile: pathlib.PosixPath, path to fasta file
    :param cutoffval: int, maximum sequence length
    """
    seqsdict1, headerlist1 = find_long_seqs(fastafile1, cutoffval)
    seqsdict2, headerlist2 = find_long_seqs(fastafile2, cutoffval)

    orgset_toremove = get_orgset_to_delete(headerlist1, headerlist2)

    if not orgset_toremove:
        raise ValueError(f'No sequences above {cutoffval}. Continuing with original fastas...')
    elif len(orgset_toremove) == len(seqsdict1):
        raise ValueError(f'All sequences are longer than {cutoffval}. Set a more suitable seqlength cutoff. Exiting...') 

    seqsdict_reduced1 = del_orgseq_from_dict(seqsdict1, orgset_toremove)
    seqsdict_reduced2 = del_orgseq_from_dict(seqsdict2, orgset_toremove)

    writeout_fasta(fastafile1, seqsdict_reduced1, overwrite=True)
    print(f'Sequences reduced in {fastafile1}. Original file overwritten.')
    writeout_fasta(fastafile2, seqsdict_reduced2, overwrite=True)
    print(f'Sequences reduced in {fastafile2}. Original file overwritten.')
