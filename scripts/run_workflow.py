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

tasks = {'findrefseqs': ('find refseq fasta files', findrefseqs),
         'phmmer': ('run phmmer on refseq', run_phmmer),
         'parsephmmer': ('parse phmmer into keyfile', parse_accid_phmmerlog)} 

class InputConfig():
    """Reads input from pathfile and datafile"""

    def __init__(self):
        """Initiates the class"""
        self.fastapath = ''
        self.dbpath = ''
        self.phmmerpath = ''

        self._read_paths()

        self.pdbid = ''
        self.refseq1 = ''
        self.refseq2 = ''
        self.keyfile1 = ''
        self.keyfile2 = ''
        self.matchedkeyfile1 = ''
        self.matchedkeyfile2 = ''
        self.alnfile1 = ''
        self.alnfile2 = ''
        self.jointalnfile = ''
        self.mfdcaoutfile = '' 


    def _read_paths(self):
        """Reads paths from paths.txt"""
        with open('../testdata/paths.txt', 'r') as p:
            paths = p.readlines()
            paths = [Path(path.strip()) for path in paths]
        self.fastapath = paths[0]
        self.dbpath = paths[1]
        self.phmmerpath = paths[2] 

    def _read_inputs(self):
        """Reads file inputs from config.txt"""
        with open('../testdata/config.txt', 'r') as c:
            for line in c.readlines():
                if line.startswith('pdbid'):
                    self.pdbid=line.strip().split('=')[1]

                elif line.startswith('refseq1'):
                    self.refseq1=line.strip().split('=')[1]

                elif line.startswith('refseq2'):
                    self.refseq2=line.strip().split('=')[1]

                elif line.startswith('keyfile1'):
                    self.keyfile1=line.strip().split('=')[1]

                elif line.startswith('keyfile2'):
                    self.keyfile2=line.strip().split('=')[1]

                elif line.startswith('matchedkeyfile1'):
                    self.matchedkeyfile1=line.strip().split('=')[1]

                elif line.startswith('matchedkeyfile2'):
                    self.matchedkeyfile2=line.strip().split('=')[1]

                elif line.startswith('alignment1='):
                    self.alnfile1=line.strip().split('=')[1]

                elif line.startswith('alignment2'):
                    self.alnfile2=line.strip().split('=')[1]

                elif line.startswith('jointalignment'):
                    self.jointalnfile1=line.strip().split('=')[1]

                elif line.startswith('mfdcaoutfile'):
                    self.mfdcaoutfile=line.strip().split('=')[1]


def read_in_datafile(datafile):
    pass


def findrefseqs(pdbid, fastapath):
    """Runs find_refseq_files"""
    refseqpaths = []
    try:
        refseqpaths = find_refseq_files(pdbid, fastapath)
        print(f'Found these refseq files:\n{refseqpaths}\n')
    except Exception as e:
        print(e)
    finally:
        return refseqpaths


def runphmmer():
    pass


def parsephmmer():
    pass


def processphmmer():
    pass


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
    parser.add_argument("-r", "--redo", help="True/False to re-parse out keyfile")
    args = parser.parse_args()

    if not args.redo:
        run_workflow(args.pdbid, args.pathfile)
    else:
        run_workflow(args.pdbid, args.pathfile, args.redo)
