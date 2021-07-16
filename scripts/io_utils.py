#!/usr/bin/env python3
"""
io_utils.py

Collection of file writing and checking utility functions.

Used with the eukaryotic dimer dca method
"""

import warnings
from pathlib import Path

def refseq_formatter(pdbid):
    """Returns formatted refseq fasta 
    str expression as a regex matching string
    
    :param pdbid: str, 4-letter pdbid
    
    :returns: str
    """
    return f'{pdbid}*refseq.fasta'


def phmmerlog_formatter(seqpath):  
    """Formats phmmer output files for given seqfile.
    
    :param seqpath: pathlib.PosixPath

    :returns: str, outfile name
    """
    return f'{seqpath.stem}_phmmer.log'


def keyfile_formatter(pathtophmmerlog):
    """Returns formatted keyfile from phmmerlogfile
    
    :param pathtophmmerlog: pathlib.PosixPath
    
    :return: str
    """
    return f'{pathtophmmerlog.stem}.keyfile'


def matched_keyfile_formatter(pathtophmmerlog):
    """Returns formatted keyfile from phmmerlogfile
    
    :param pathtophmmerlog: pathlib.PosixPath
    
    :return: str
    """
    return f'{pathtophmmerlog.stem}_matched.keyfile'


def get_globbed_list(pathtodir, target):
    """Searches directory for files matching
    a certain target pattern.

    :param pathtodir: pathlib.PosixPath 
    :param target: str with expression to match

    :returns: list of matching paths
    """
    p = pathtodir.expanduser()
    return list(p.glob(target))


def does_target_exist(pathtotarget, targettype):
    """Given directory and a target (file or dir),
    returns True if it exists, False if not.

    :param pathtotarget: pathlib.PosixPath

    :returns: bool
    """
    if targettype == 'file':
        return pathtotarget.is_file()
    elif targettype == 'dir':
        return pathtotarget.is_dir()
    else:
        return pathtotarget.exists()


def readin_list(filepathoflist):
    """Reads in a file with \n spaced items
    into a list"""
    with open(filepathoflist, 'r') as r:
        lines = r.readlines()
    return [item.strip() for item in lines]


def writeout_list(listofitems, outpath):
    """Writes out contents of a list, newline-separated
    to a file specified by outpath.

    :param listofitems: list to write out
    :param outpath: pathlib.PosixPath
    """
    with open(outpath, 'w+') as f:
        f.write('\n'.join(listofitems))

    print(f'File: {outpath.name} written into dir: {outpath.parent}')


def parse_fasta(lines):
    """Parses a fasta file
    
    :param lines: list of lines
    :returns res: dict of key + seq
    """
    res = {}  
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            label = line[1:]
            res[label] = []
        else:
            res[label].append(line)
    for k, v in res.items():
        res[k] = ''.join(v)
    return res


def fa_todict(fastafile):  
    """Opens fasta file and returns seqdict"""
    with open(fastafile, 'r') as f:
        lines = f.readlines()
        seqsdict = parse_fasta(lines)
    return seqsdict

def writeout_fasta(somepath, somedict, overwrite=False, addict={}):
    """Writes dict of seqs to a fastafile"""
    if not overwrite:
        with open(somepath, 'a+') as f:
            for k, v in somedict.items():
                f.write(''.join(['>', k, '\n']))
                f.write(''.join([v, '\n']))
    elif overwrite and addict:
        with open(somepath, 'w') as f:
            for k, v in somedict.items():
                f.write(''.join(['>', k, '\n']))
                f.write(''.join([v, '\n']))
            for k, v in addict.items():
                f.write(''.join(['>', k, '\n']))
                f.write(''.join([v, '\n']))
    return
