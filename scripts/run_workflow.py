#!/usr/bin/env python3
"""
run_workflow.py

Runs eukdimerdca workflow.

INPUT:

1. pdbid
2. fastapath
3. databasepath
4. phmmerpath
"""

from pathlib import Path

from find_refseq_files import *
from io_utils import *
from run_phmmer import *
from parse_accid_phmmerlog import *

tasknames = ['findrefseqs', 'phmmer', 'parsephmmer']

tasks = {'findrefseqs': ('find refseq fasta files', find_refseq_files),
         'phmmer': ('run phmmer on refseq', run_phmmer),
         'parsephmmer': ('parse phmmer into keyfile', parse_accid_phmmerlog)} 

def read_in_pathfile(pathfile):
    """Gets fasta, database, and phmmer paths from text file.
    Returns list of pathlib.PosixPaths"""
    with open(pathfile, 'r') as p:
        paths = p.readlines()
        paths = [Path(path.strip()) for path in paths]
    return paths


def run_workflow(pdbid, pathfile, redo=False):
    """Runs eukdimerdca workflow"""

    paths = read_in_pathfile(pathfile)
    fastapath = paths[0]
    databasepath = paths[1]
    phmmerpath = paths[2]

    refseqpaths = find_refseq_files(pdbid, fastapath)
    print(f'Found these refseq files:\n{refseqpaths}\n')

    if refseqpaths:
        for seqpath in refseqpaths:
            print(f'------> Working on {seqpath}')
            outfilepath = run_phmmer(databasepath, seqpath, phmmerpath) 
            if outfilepath.is_file():
                parse_accid_phmmerlog(outfilepath, phmmerpath, overwrite=redo)

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] pdbid pathfile --redo")
    parser.add_argument("pdbid", help="4-letter PDB ID")
    parser.add_argument("pathfile", help="file with fasta, database, phmmer paths")
    parser.add_argument("-r", "--redo", help="overwrite keyfile")
    args = parser.parse_args()

    if not args.redo:
        run_workflow(args.pdbid, args.pathfile)
    else:
        run_workflow(args.pdbid, args.pathfile, args.redo)
