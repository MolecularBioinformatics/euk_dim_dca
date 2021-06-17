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
from process_phmmerhits import *


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
        self.logfile1 = ''
        self.logfile2 = ''
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
                attr = line.strip().split('=')
                if attr[0] in self.__dict__.keys():
                    if attr[0] == 'pdbid':
                        self.__dict__[attr[0]] = attr[1]
                    else:
                        self.__dict__[attr[0]] = Path(attr[1])

    def update_config_var(self):
        """Updates config file with attributes"""
        with open('../testdata/config.txt', 'w+') as c:
            for key, item in self.__dict__.items():
                if not key.endswith('path'):
                    c.write(f'{key}={item}\n') 


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


def runphmmer(databasepath, seqpath, phmmerpath):
    """Runs phmmer on seq"""
    try:
        outfilepath = run_phmmer(databasepath, seqpath, phmmerpath)
        return outfilepath
    except Exception as e:
        print(e)


def parsephmmer(outfilepath, phmmerpath, overwrite):
    """Parses phmmer log to keyfile"""
    try:
        keyfile = parse_accid_phmmerlog(outfilepath, phmmerpath, overwrite)
        return keyfile
    except FileNotFoundError as fnotfound:
        print(fnotfound)


def processphmmer(phmmerpath, pdbid, minhits=100, maxhits=600):
    """Checks total number of hits in keyfile, matches organisms.
       Returns processed keyfile"""
    try:
        file1, file2 = process_phmmerhits(phmmerpath, pdbid, minhits, maxhits)
        return file1, file2
    except:
        raise Exception('Unable to process keyfiles.')

tasknames = ['findrefseqs', 'runphmmer', 'parsephmmer', 'processphmmer']

tasks = {'findrefseqs': ('find refseq fasta files', findrefseqs),
         'runphmmer': ('run phmmer on refseq', runphmmer),
         'parsephmmer': ('parse phmmer into keyfile', parsephmmer),
         'processphmmer': ('process keyfile and match organisms', processphmmer)} 

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

    elif taskname == 'runphmmer': 
        print(f'{taskname}: {tasks[taskname][0]}')
        if not (IC.refseq1 and IC.refseq2):
            IC.refseq1, IC.refseq2 = findrefseqs(IC.pdbid, IC.fastapath)
        IC.logfile1 = runphmmer(IC.dbpath, IC.refseq1, IC.phmmerpath)
        IC.logfile2 = runphmmer(IC.dbpath, IC.refseq2, IC.phmmerpath)

    elif taskname == 'parsephmmer':
        print(f'{taskname}: {tasks[taskname][0]}')
        if not (IC.logfile1.is_file() and IC.logfile2.is_file()):
            raise FileNotFoundError(f'Phmmer logfiles do not exist:\n{IC.logfile1}\n{IC.logfile2}')
        else:
            IC.keyfile1 = parsephmmer(IC.logfile1, IC.phmmerpath, redo)
            IC.keyfile2 = parsephmmer(IC.logfile2, IC.phmmerpath, redo)

    elif taskname == 'processphmmer':
        print(f'{taskname}: {tasks[taskname][0]}')
        if not (IC.keyfile1.is_file() and IC.keyfile2.is_file()):
            raise FileNotFoundError(f'Phmmer unprocessed keyfiles do not exist:\n{IC.keyfile1}\n{IC.keyfile2}')
        else:
            IC.matchedkeyfile1, IC.matchedkeyfile2 = processphmmer(IC.phmmerpath, IC.pdbid)
        
    IC.update_config_var()

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] taskname --redo")
    parser.add_argument("taskname", help="task to run: findrefseqs, runphmmer, parsephmmer")
    parser.add_argument("-r", "--redo", help="True/False to re-parse out keyfile")
    args = parser.parse_args()

    if not args.redo:
        run_workflow(args.taskname)
    else:
        run_workflow(args.taskname, args.redo)
