#!/usr/bin/env python3
"""
plot_seqlengths.py

Makes basic histogram of sequence lengths in a fasta file.
Useful to determine at which sequence length to exclude sequences.
"""

from pathlib import Path
from io_utils import parse_fasta, fa_todict
import matplotlib.pyplot as plt

def get_lengthsdict(fastafile):
    """Returns a dict of counts of seq lengths

    :param fastafile: pathlib.PosixPath
    """
    seqsdict = fa_todict(fastafile)

    lensdict = {}
    for val in seqsdict.values():
        seqlen = len(val)
        if seqlen not in lensdict.keys():
            lensdict[seqlen] = 1
        else:
            lensdict[seqlen] += 1
    return lensdict

def plot_counts_hist(dictofcounts):
    """Makes histogram of counts of seqlengths"""
    plt.bar(list(dictofcounts.keys()), list(dictofcounts.values()), color='g')
    plt.show()

def plot_seqlengths(fastafile):
    """ """
    lengthsdict = get_lengthsdict(fastafile)
    plot_counts_hist(lengthsdict)
