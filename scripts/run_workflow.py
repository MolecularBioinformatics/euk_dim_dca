#!/usr/bin/env python3
"""
run_workflow.py Runs eukdimerdca workflow.

INPUT:

taskname = list of tasknames to run
redo = boolean of whether or not to rerun tasks
"""

from pathlib import Path
import sys

from find_refseq_files import *
from io_utils import *
from run_phmmer import *
from parse_accid_phmmerlog import *
from process_phmmerhits import *
from run_easel_getseqs import *
from process_easelseqs import *
from align_seqs import *
from process_alnseqs import *
from reduce_seq_set import *
from run_dca import *

class InputConfig():
    """Object to store names of intermediate files.
    Reads input from pathfile and datafile"""

    def __init__(self, config, paths):
        """Initiates the class"""

        # input paths
        self.dbpath = ''
        self.easelpath = ''

        # output paths
        self.fastapath = ''
        self.phmmerpath = ''
        self.keyfilepath = ''
        self.alnpath = ''
        self.dcapath = ''

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
            for line in p.readlines():
                path = line.strip().split('=')
                if path[0] in self.__dict__.keys():
                    self.__dict__[path[0]] = Path(path[1])

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
    :param redo: T/F - not incorporated, function searches every time
    :returns icObject: InputConfig object with updated attributes.
    """
    try:
        refseqpaths = find_refseq_files(icObject.pdbid, icObject.fastapath)
        print(f'Found these refseq files:\n{refseqpaths}')
        icObject.refseq1=refseqpaths[0]
        icObject.refseq2=refseqpaths[1]
    except TypeError as pdberr:
        print(pdberr)
    except FileNotFoundError as nofileerr:
        print(nofileerr)
    except ValueError as valerr:
        print(valerr)
    finally:
        return icObject


def runphmmer(icObj, rerun): 
    """Runs phmmer on seq.
    Takes and returns an InputConfigObj."""
    try:
        icObj.logfile1 = run_phmmer(icObj.dbpath, icObj.refseq1, icObj.phmmerpath, rerun)
    except FileNotFoundError as e:
        print(e)
    except ValueError as valerr:
        print(valerr)
    try:
        icObj.logfile2 = run_phmmer(icObj.dbpath, icObj.refseq2, icObj.phmmerpath, rerun)
    except FileNotFoundError as e:
        print(e)
    except ValueError as valerr:
        print(valerr)
    finally:
        return icObj


def parsephmmer(icObj, overwrite):
    """Parses phmmer log to keyfile.
    Takes and returns and InputConfigObj.
    Overwrite is a bool."""
    try:
        icObj.keyfile1 = parse_accid_phmmerlog(icObj.logfile1, icObj.keyfilepath, overwrite)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as fileempty:
        print(filempty)
    try:
        icObj.keyfile2 = parse_accid_phmmerlog(icObj.logfile2, icObj.keyfilepath, overwrite)
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
        icObj.matchedkeyfile1, icObj.matchedkeyfile2 = process_phmmerhits(icObj.keyfilepath, icObj.keyfile1, icObj.keyfile2, minhits, maxhits, overwrite)
    except ValueError as valerr: 
        raise ValueError('Unable to process keyfiles.') from valerr
    finally:
        return icObj


def runeasel(icObj, rerun):
    """Runs easel on a keyfile.
    Extracts sequences from a db."""
    try:
        icObj.eslfastafile1 = run_easel_iterate(icObj.easelpath, icObj.dbpath, icObj.fastapath, icObj.matchedkeyfile1, rerun)
    except FileNotFoundError as e:
        print(e)
    try:
        icObj.eslfastafile2 = run_easel_iterate(icObj.easelpath, icObj.dbpath, icObj.fastapath, icObj.matchedkeyfile2, rerun)
    except FileNotFoundError as e:
        print(e)
    return icObj


def processeasel(icObj, redo):
    """Processes easel extracted fasta.
    Keeps only seqs with common orgs."""
    easelerrfile = Path(f"{icObj.fastapath}/{icObj.pdbid}.easelerror")
    if not easelerrfile.is_file():
        print('No easel error file found. Continuting with original matched fastas')
    else:
        try:
            orgs = process_easelseqs(easelerrfile, icObj.eslfastafile1, redo)
        except FileNotFoundError as fnotfound:
            print(fnotfound)
        except ValueError as valerr:
            print(valerr)
        try:
            orgs = process_easelseqs(easelerrfile, icObj.eslfastafile2, redo)
            print(f'Orgs removed: {orgs}')
        except FileNotFoundError as fnotfound:
            print(fnotfound)
        except ValueError as valerr:
            print(valerr)
    return icObj


def reduceseqset(icObj, redo):
    """Reduces a set of sequences to remove
    seqs that are too long, hopefully makes
    alignment step more manageable."""

    # TODO: doesn't use redo yet
    # TODO: reduction scheme is brute force right now, only throws error if 0 seqs are left
    
    maxlength = 1600 #kind of arbitrary!

    try:
        icObj.eslfastafile1, icObj.eslfastafile2 = reduce_seq_set(icObj.eslfastafile1, icObj.eslfastafile2, maxlength)
    except ValueError as valerr:
        print(valerr)
    except RuntimeError as rerr:
        print(rerr)
        sys.exit()

    return icObj

def alignseqs(icObj, realign):
    """Runs muscle to align sequences.
    Takes and returns an InputConfigObj."""
    try:
        icObj.alnfile1 = run_muscle(icObj.eslfastafile1, icObj.refseq1, icObj.alnpath, realign)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as valerr:
        print(valerr)
    try:
        icObj.alnfile2 = run_muscle(icObj.eslfastafile2, icObj.refseq2, icObj.alnpath, realign)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as valerr:
        print(valerr)
    return icObj


def processalignment(icObj, redo):
    """Matches aligned sequences by organism.
    Joins (horizontally concatenates) matched sequences
    to make a joint alignment file for DCA"""
    try:
        icObj.jointalnfile = process_alnseqs(icObj.alnfile1, icObj.alnfile2, icObj.refseq1, icObj.refseq2, icObj.alnpath, redo)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as valerr:
        print(valerr)

    return icObj


def rundca(icObj, redo):
    """Runs dca on a joint alignment.
    Deposits scores into a scores.dat file."""
    try:
        icObj.mfdcaoutfile = run_dca(icObj.jointalnfile, icObj.dcapath, redo)
    except FileNotFoundError as fnotfound:
        print(fnotfound)
    except ValueError as valerr:
        print(valerr)
    return icObj


TASKNAMES = ['all', 'findrefseqs', 'runphmmer', 'parsephmmer', 'processphmmer', 'runeasel', 'processeasel', 'reduceseqset', 'alignseqs', 'processalignment', 'rundca'] 

TASKS = {'findrefseqs': ('1. find refseq fasta files\n', findrefseqs),
         'runphmmer': ('2. run phmmer on refseq\n', runphmmer),
         'parsephmmer': ('3. parse phmmer into keyfile\n', parsephmmer),
         'processphmmer': ('4. process keyfile and match organisms\n', processphmmer),
         'runeasel': ('5. runs easel extract to get seqs from db\n', runeasel),
         'processeasel': ('6. process easel extracted seqs based on organisms\n', processeasel),
         'reduceseqset': ('7. removes sequences that are too long from seq set\n', reduceseqset), 
         'alignseqs': ('8. aligns sequences with muscle\n', alignseqs),
         'processalignment': ('9. matches and joins aligned sequences\n', processalignment),
         'rundca': ('10. runs DCA on joint aligned sequences\n', rundca)} 

def run_workflow(configf, pathsf, tasknamelist, redo):
    """Runs eukdimerdca workflow"""
    try:
        ic = InputConfig(configf, pathsf) 
    except IOError:
        ic = None

    if 'all' in tasknamelist:
        tasknamelist = TASKNAMES[1:]

    for taskname in tasknamelist: 
        if taskname not in TASKNAMES:
            raise ValueError(f'{taskname} not a valid task. Try again.')
        print(f'--- {taskname} --- {TASKS[taskname][0]}')
        torun = TASKS[taskname][1] 
        ic = torun(ic, redo)
        print('\n')
        ic.update_config_var(configf)

if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] configfile pathfile taskname [taskname, ...] --redo")
    parser.add_argument("configfile", help="path to config.txt file")
    parser.add_argument("pathfile", help="path to paths.txt file")
    parser.add_argument("taskname", nargs = '+', help="task to run: findrefseqs, runphmmer, parsephmmer, processphmmer, runeasel, processeasel, reduceseqset, alignseqs, processalignment, rundca")
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
        tasklist = [args.taskname]
    else:
        tasklist = args.taskname
    run_workflow(configfile, pathfile, tasklist, redoflag)
