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
from run_easel_getseqs import *


class InputConfig():
    """Reads input from pathfile and datafile"""

    def __init__(self, config, paths):
        """Initiates the class"""

        self.fastapath = ''
        self.dbpath = ''
        self.phmmerpath = ''
        self.easelpath = ''

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
        self.eslfastafile1 = ''
        self.eslfastafile2 = ''
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
        self.easelpath = paths[3]

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


def findrefseqs(icObject, redo): 
    """Finds two refseqs for pdbid.

    :param icObject: InputConfig object
    :returns icObject: InputConfig object with updated attributes.
    """
    try:
        refseqpaths = find_refseq_files(icObject.pdbid, icObject.fastapath)
        print(f'Found these refseq files:\n{refseqpaths}')
        icObject.refseq1=refseqpaths[0]
        icObject.refseq2=refseqpaths[1]
    except Exception as e:
        print(e)
    finally:
        return icObject


def runphmmer(icObj, rerun):  # test on actually running phmmer
    """Runs phmmer on seq.
    Takes and returns an InputConfigObj."""
    try:
        icObj.logfile1 = run_phmmer(icObj.dbpath, icObj.refseq1, icObj.phmmerpath, rerun)
    except Exception as e:
        print(e)
    try:
        icObj.logfile2 = run_phmmer(icObj.dbpath, icObj.refseq2, icObj.phmmerpath, rerun)
    except Exception as e:
        print(e)
    finally:
        return icObj


def parsephmmer(icObj, overwrite):
    """Parses phmmer log to keyfile.
    Takes and returns and InputConfigObj.
    Overwrite is a bool."""
    try:
        icObj.keyfile1 = parse_accid_phmmerlog(icObj.logfile1, icObj.phmmerpath, overwrite)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as fileempty:
        print(filempty)
    try:
        icObj.keyfile2 = parse_accid_phmmerlog(icObj.logfile2, icObj.phmmerpath, overwrite)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as fileempty:
        print(filempty)
    finally:
        return icObj


def processphmmer(icObj, overwrite):
    """Checks total number of hits in keyfile, matches organisms.
       Returns processed keyfile"""

    minhits = 100
    maxhits = 600

    try:
        icObj.matchedkeyfile1, icObj.matchedkeyfile2 = process_phmmerhits(icObj.phmmerpath, icObj.pdbid, minhits, maxhits, overwrite)
    except ValueError as valerr: 
        raise ValueError('Unable to process keyfiles.') from valerr
    finally:
        return icObj


def runeasel(icObj, rerun):
    """Runs easel on a keyfile.
    Extracts sequences from a db."""
    try:
        icObj.eslfastafile1 = run_easel(icObj.easelpath, icObj.dbpath, icObj.phmmerpath, icObj.matchedkeyfile1, rerun)
    except Exception as e:
        print(e)
    try:
        icObj.eslfastafile2 = run_easel(icObj.easelpath, icObj.dbpath, icObj.phmmerpath, icObj.matchedkeyfile2, rerun)
    except Exception as e:
        print(e)
    return icObj


TASKNAMES = ['all', 'findrefseqs', 'runphmmer', 'parsephmmer', 'processphmmer', 'runeasel'] 

TASKS = {'findrefseqs': ('1. find refseq fasta files\n', findrefseqs),
         'runphmmer': ('2. run phmmer on refseq\n', runphmmer),
         'parsephmmer': ('3. parse phmmer into keyfile\n', parsephmmer),
         'processphmmer': ('4. process keyfile and match organisms\n', processphmmer),
         'runeasel': ('5. runs easel extract to get seqs from db\n', runeasel)} 

def run_workflow(configf, pathsf, tasknamelist, redo):
    """Runs eukdimerdca workflow"""

    try:
        ic = InputConfig(configf, pathsf) 
    except IOError:
        ic = None

    for taskname in tasknamelist: 
        if taskname not in tasknames:
            raise ValueError(f'{taskname} not a valid task. Try again.')
        else:
            print(f'--- {taskname} --- {tasks[taskname][0]}')
            torun = tasks[taskname][1] 
            ic = torun(ic, redo)
            print('\n')
        ic.update_config_var(configf)

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
    redoflag = args.redo

    if redoflag is None:
        redoflag = False
    elif redoflag == 'False':
        redoflag = False
    elif redoflag == 'True':
        redoflag = True

    if isinstance(args.taskname, str):
        singletask = [args.taskname]
        run_workflow(configfile, pathfile, singletask, redoflag)
    elif isinstance(args.taskname, list):
        tasklist = args.taskname
        run_workflow(configfile, pathfile, tasklist, redoflag)
    else:
        raise ValueError('Invalid input.')
