#!/usr/env/bin python3
"""
config_file_from_pdbid.py

Prepares an empty configuration file formatted
for the DCA workflow. Takes in a PDBID and path.
"""

from pathlib import Path

from io_utils import does_target_exist

def find_config_file(somepath, pdbid):
    """Searches for config file.
    If file exists raises error.

    :param somepath: pathlib.PosixPath
    :param pdbid: str
    """
    configf = f'config_{pdbid}.txt'
    searchpath = somepath / configf

    if does_target_exist(searchpath, 'file'):
        raise FileExistsError(f'{searchpath} already exists.\nDo you want to overwrite?')
    return searchpath
    

def writeout_config_file(configpath, pdbid):
    """Writes out configuration file 
    with the correct fields into config_pdbid.txt

    :param configpath: pathlib.PosixPath
    :param pdbid: str
    """
    try:
        filepath = find_config_file(configpath, pdbid)
    except FileExistsError as fee:
        print(fee)

    nameslist = ['refseq1',
                 'refseq2',
                 'logfile1',
                 'logfile2',
                 'keyfile1',
                 'keyfile2',
                 'matchedkeyfile1',
                 'matchedkeyfile2',
                 'eslfastafile1',
                 'eslfastafile2',
                 'alnfile1',
                 'alnfile2',
                 'jointalnfile',
                 'mfdcaoutfile',]

    with open(filepath, 'w+') as f:
        f.write(f'pdbid={pdbid}\n')
        for name in nameslist:
            f.write(f'{name}=\n')
