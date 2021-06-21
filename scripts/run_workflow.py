#!/usr/bin/env python3
"""
run_workflow.py Runs eukdimerdca workflow.

INPUT:

taskname = list of tasknames to run
redo = boolean of whether or not to rerun tasks
"""

from pathlib import Path

from find_refseq_files import *
from io_utils import *
from run_phmmer import *
from parse_accid_phmmerlog import *
from process_phmmerhits import *


class InputConfig():
    """Reads input from pathfile and datafile"""

    def __init__(self, config, paths):
        """Initiates the class"""

        self.fastapath = ''
        self.dbpath = ''
        self.phmmerpath = ''

        self._read_paths(paths)

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

        self._read_inputs(config)


    def _read_paths(self, paths):
        """Reads paths from paths.txt"""
        with open(paths, 'r') as p:
            paths = p.readlines()
            paths = [Path(path.strip()) for path in paths]
        self.fastapath = paths[0]
        self.dbpath = paths[1]
        self.phmmerpath = paths[2] 

    def _read_inputs(self, config):
        """Reads file inputs from config.txt"""
        with open(config, 'r') as c:
            for line in c.readlines():
                attr = line.strip().split('=')
                if attr[0] in self.__dict__.keys():
                    if attr[0] == 'pdbid':
                        self.__dict__[attr[0]] = attr[1]
                    else:
                        self.__dict__[attr[0]] = Path(attr[1])

    def update_config_var(self, config):
        """Updates config file with attributes"""
        with open(config, 'w+') as c:
            for key, item in self.__dict__.items():
                if not key.endswith('path'):
                    c.write(f'{key}={item}\n') 


def findrefseqs(ICObject, redo=False): 
    """Finds two refseqs for pdbid.

    :param ICObject: InputConfig object
    :returns ICObject: InputConfig object with updated attributes.
    """
    try:
        refseqpaths = find_refseq_files(ICObject.pdbid, ICObject.fastapath)
        print(f'Found these refseq files:\n{refseqpaths}')
        ICObject.refseq1=refseqpaths[0]
        ICObject.refseq2=refseqpaths[1]
    except Exception as e:
        print(e)
    finally:
        return ICObject


def runphmmer(ICObj, rerun=False):  # test on actually running phmmer
    """Runs phmmer on seq.
    Takes and returns an InputConfigObj."""
    try:
        ICObj.logfile1 = run_phmmer(ICObj.dbpath, ICObj.refseq1, ICObj.phmmerpath, rerun)
    except Exception as e:
        print(e)
    try:
        ICObj.logfile2 = run_phmmer(ICObj.dbpath, ICObj.refseq2, ICObj.phmmerpath, rerun)
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


def processphmmer(ICObj, overwrite=False):
    """Checks total number of hits in keyfile, matches organisms.
       Returns processed keyfile"""

    minhits = 100
    maxhits = 600

    try:
        ICObj.matchedkeyfile1, ICObj.matchedkeyfile2 = process_phmmerhits(ICObj.phmmerpath, ICObj.pdbid, minhits, maxhits, overwrite)
    except:
        raise Exception('Unable to process keyfiles.')
    finally:
        return ICObj


tasknames = ['all', 'findrefseqs', 'runphmmer', 'parsephmmer', 'processphmmer']

tasks = {'findrefseqs': ('1. find refseq fasta files\n', findrefseqs),
         'runphmmer': ('2. run phmmer on refseq\n', runphmmer),
         'parsephmmer': ('3. parse phmmer into keyfile\n', parsephmmer),
         'processphmmer': ('4. process keyfile and match organisms\n', processphmmer)} 

def run_workflow(configf, pathsf, tasknamelist, redo=False):
    """Runs eukdimerdca workflow"""

    try:
        IC = InputConfig(configf, pathsf)
    except IOError:
        IC = None

    for taskname in tasknamelist:  # agreement between taskname as str vs list
        if taskname not in tasknames:
            raise ValueError(f'{taskname} not a valid task. Try again.')
        else:
            print(f'--- {taskname} --- {tasks[taskname][0]}')
            torun = tasks[taskname][1] 
            IC = torun(IC, redo)
            print('\n')
        IC.update_config_var(configf)

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] taskname [taskname, ...] --redo")
    parser.add_argument("configfile", help="path to config.txt file")
    parser.add_argument("pathfile", help="path to paths.txt file")
    parser.add_argument("taskname", nargs = '+', help="task to run: findrefseqs, runphmmer, parsephmmer")
    parser.add_argument("-r", "--redo", help="True/False to re-parse out keyfile")
    args = parser.parse_args()

    configfile = args.configfile
    pathfile = args.pathfile

    if isinstance(args.taskname, str):
        singletask = [args.taskname]
        if not args.redo:
            run_workflow(configfile, pathfile, singletask)
        else:
            run_workflow(configfile, pathfile, singletask, args.redo)
    elif isinstance(args.taskname, list):
        tasklist = args.taskname
        if not args.redo:
            run_workflow(configfile, pathfile, tasklist)
        else:
            run_workflow(configfile, pathfile, tasklist, args.redo)
    else:
        raise ValueError('Invalid input.')
