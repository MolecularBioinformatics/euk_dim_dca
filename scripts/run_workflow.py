#!/usr/bin/env python3
"""
run_workflow.py Runs eukdimerdca workflow.

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


def findrefseqs(ICObject): 
    """Finds two refseqs for pdbid.

    :param ICObject: InputConfig object
    :returns ICObject: InputConfig object with updated attributes.
    """
    try:
        refseqpaths = find_refseq_files(ICObject.pdbid, ICObject.fastapath)
        print(f'Found these refseq files:\n{refseqpaths}\n')
        ICObject.refseq1=refseqpaths[0]
        ICObject.refseq2=refseqpaths[1]
    except Exception as e:
        print(e)
    finally:
        return ICObject


def runphmmer(ICObj):  # test on actually running phmmer
    """Runs phmmer on seq.
    Takes and returns an InputConfigObj."""
    try:
        ICObj.logfile1 = run_phmmer(ICObj.dbpath, ICObj.refseq1, ICObj.phmmerpath)
    except Exception as e:
        print(e)
    try:
        ICObj.logfile2 = run_phmmer(ICObj.dbpath, ICObj.refseq2, ICObj.phmmerpath)
    except Exception as e:
        print(e)
    finally:
        return ICObj


def parsephmmer(ICObj, overwrite=False):
    """Parses phmmer log to keyfile.
    Takes and returns and InputConfigObj.
    Overwrite is a bool."""
    try:
        ICObj.keyfile1 = parse_accid_phmmerlog(ICObj.logfile1, ICObj.phmmerpath, overwrite)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    try:
        ICObj.keyfile2 = parse_accid_phmmerlog(ICObj.logfile2, ICObj.phmmerpath, overwrite)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    finally:
        return ICObj


def processphmmer(ICObj, minhits=100, maxhits=600):
    """Checks total number of hits in keyfile, matches organisms.
       Returns processed keyfile"""
    try:
        ICObj.matchedkeyfile1, ICObj.matchedkeyfile2 = process_phmmerhits(ICObj.phmmerpath, ICObj.pdbid, minhits, maxhits)
    except:
        raise Exception('Unable to process keyfiles.')
    finally:
        return ICObj


tasknames = ['all', 'findrefseqs', 'runphmmer', 'parsephmmer', 'processphmmer']

tasks = {'findrefseqs': ('find refseq fasta files', findrefseqs),
         'runphmmer': ('run phmmer on refseq', runphmmer),
         'parsephmmer': ('parse phmmer into keyfile', parsephmmer),
         'processphmmer': ('process keyfile and match organisms', processphmmer)} 

def run_workflow(tasknamelist, redo=False):
    """Runs eukdimerdca workflow"""

    try:
        IC = InputConfig()
    except IOError:
        IC = None

    for taskname in tasknamelist:  # agreement between taskname as str vs list
        if taskname not in tasknames:
            raise ValueError(f'{taskname} not a valid task. Try again.')
        else:
            torun = tasks[taskname][1] 
            IC = torun(IC)
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
