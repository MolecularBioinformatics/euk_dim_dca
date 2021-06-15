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

        self._read_inputs()


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


def findrefseqs(pdbid, fastapath):
    """Runs find_refseq_files"""
    refseqpaths = []
    try:
        refseqpaths = find_refseq_files(pdbid, fastapath)
        print(f'Found these refseq files:\n{refseqpaths}\n')
    except Exception as e:
        print(e)
    finally:
        return refseqpaths[0], refseqpaths[1]


tasknames = ['findrefseqs', 'phmmer', 'parsephmmer']

tasks = {'findrefseqs': ('find refseq fasta files', findrefseqs),
         'phmmer': ('run phmmer on refseq', run_phmmer),
         'parsephmmer': ('parse phmmer into keyfile', parse_accid_phmmerlog)} 

def run_workflow(taskname, redo=False):
    """Runs eukdimerdca workflow"""

    try:
        IC = InputConfig()
    except IOError:
        IC = None


    if taskname not in tasknames:
        raise ValueError(f'{taskname} not a valid task. Try again.')

    elif taskname == 'findrefseqs':
        print(f'{taskname}: {tasks[taskname][0]}')
        IC.refseq1, IC.refseq2 = findrefseqs(IC.pdbid, IC.fastapath)

    elif taskname == 'phmmer': 
        print(f'{taskname}: {tasks[taskname][0]}')
        if not (IC.refseq1 and IC.refseq2):
            IC.refseq1, IC.refseq2 = findrefseqs(IC.pdbid, IC.fastapath)
        else:
            

#    if refseqpaths:
#        for seqpath in refseqpaths:
#            print(f'------> Working on {seqpath}')
#            outfilepath = run_phmmer(databasepath, seqpath, phmmerpath) 
#            if outfilepath.is_file():
#                parse_accid_phmmerlog(outfilepath, phmmerpath, overwrite=redo)

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
