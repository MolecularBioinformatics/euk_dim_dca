#!/usr/env/bin python3
"""
parse_accid_phmmerlog.py

Parses a phmmerlog file, extracts all accids
of sequences up until the inclusion threshold.
"""

from pathlib import Path

from io_utils import writeout_list


def keyfile_formatter(pathtophmmer):
    """Returns formatted keyfile from phmmerlogfile"""
    return f'{pathtophmmer.stem}.keyfile'


def get_accidlist(pathtophmmer):
    """Returns list of accession IDs from
    phmmer logfile

    :param pathtophmmer: pathlib.PosixPath

    :returns accidlist: list or None
    """

    accidlist=[]

    with open(pathtophmmer, 'r') as ph:
        while True:
            line = ph.readline()
            text = line.strip()
            if "inclusion threshold" in line:
                break 
            if text:
                firstletter = text[0]
                if firstletter.isnumeric():
                    accidlist.append(text.split()[8])
        return accidlist


def parse_accid_phmmerlog(pathtophmmer, filename, outpath, overwrite=False):
    """Parses out accids from phmmerlog
    into a keyfile.

    :param pathtophmmer: pathlib.PosixPath
    :param outpath: pathlib.PosixPath, file and path to output
    """
    keyfilepath = outpath.joinpath(filename)

    if not pathtophmmer.is_file():
        raise FileNotFoundError(f'File {pathtophmmer} not found!')
    else:
        acclist = get_accidlist(pathtophmmer)
        if not keyfilepath.is_file():
            writeout_list(acclist, keyfilepath) 
        elif overwrite:
            writeout_list(acclist, keyfilepath) 
        else:
            print(f'{keyfilepath} exists already. Overwrite with overwrite=True.')
