#!/usr/bin/env python3
"""
collect_refseqs.py

Takes in a file with list of pdbids and a path to the refseq files.
Puts all the refseq files for each pdbid in a fasta file.

Intended to be used for clustering
"""

from pathlib import Path

from io_utils import get_globbed_list, parse_fasta, writeout_fasta


def readin_pdbids(pathtofile):
    """Reads in a text file with list of pdbids or files
    beginning with pdbids. Saves into a list and returns list.

    :param pathtofile: pathlib.PosixPath
    :returns pdbidlist: list of 4-letter (lowercase) pdbids
    """
    pdbidlist=[]
    if not pathtofile.is_file():
        raise FileNotFoundError(f'{pathtofile} not a file. Check input')
    elif pathtofile.stat().st_size == 0:
        raise ValueError(f'{pathtofile} is empty')
    with open(pathtofile, 'r') as f:
        for line in f.readlines():
            pdbidlist.append(line.strip()[0:4].lower())
    return pdbidlist


def scan_dir_for_refseq(pathtorefseqs, pdbidlist):
    """Looks through folder for refseq files for each 
    pdbid. Saves the paths to refseq files in a list.

    :param pathtorefseqs: pathlib.PosixPath
    :param pdbidlist: list of 4-letter pdbids

    :returns refseqpathlist: list of pathlib.PosixPath
    """
    refseqpathlist = []
    if not pdbidlist:
        raise ValueError('No pdbids in list. Check input file.')
    for pdbid in pdbidlist:
        targregex = f'{pdbid}*refseq.fasta'
        refseqlist = get_globbed_list(pathtorefseqs, targregex)
        refseqpathlist += refseqlist
    return refseqpathlist

def collect_refseqs(pathto_pdblist, pathto_refseqs, outfilename, outpath, printout=False):
    """Reads in pdbids, searches for matching refseqs and puts
    them into a multisequence fasta file.

    :param pathto_pdblist: pathlib.PosixPath
    :param pathto_refseqs: pathlib.PosixPath
    :param outfilename: str, filename of multiseq fasta file
    :param outpath: pathlib.PosixPath
    """
    pdbs_list = readin_pdbids(pathto_pdblist)
    refseqs_list = scan_dir_for_refseq(pathto_refseqs, pdbs_list)
    collectedseqdict = {}
    if not refseqs_list:
        raise ValueError('No refseq files found for pdbids in {pathto_pdblist}.')
    for refseqfile in refseqs_list:
        with open(refseqfile, 'r') as fa:
            lines = fa.readlines()
        singleseqdict = parse_fasta(lines)
        key = next(iter(singleseqdict))
        collectedseqdict[key] = singleseqdict[key]
    if printout:
        print(f'Creating multiseq fasta from these files: {refseqs_list}')
    writeout_fasta(outpath.joinpath(outfilename), collectedseqdict)
