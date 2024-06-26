#!/usr/env/bin python3
"""
parse_accid_phmmerlog.py

Parses a phmmerlog file, extracts all accids
of sequences up until the inclusion threshold.

Writes these accids into a keyfile.

Returns the path to the keyfile.
"""

from pathlib import Path

from io_utils import writeout_list, keyfile_formatter


def get_accidlist(pathtophmmerlog):
    """Returns list of accession IDs from
    phmmer logfile

    :param pathtophmmerlog: pathlib.PosixPath

    :returns accidlist: list or None
    """

    accidlist=[]

    with open(pathtophmmerlog, 'r') as ph:
        for line in ph:
            text = line.strip()
            if "inclusion threshold" in text:
                break
            elif text and text[0].isnumeric():
                accidlist.append(text.split()[8])
        return accidlist


def parse_accid_phmmerlog(pathtophmmerlog, outpath, overwrite):
    """Parses out accids from phmmerlog
    into a keyfile.

    :param pathtophmmerlog: pathlib.PosixPath
    :param outpath: pathlib.PosixPath, file and path to output

    :returns keyfilepath: pathlib.PosixPath, path to keyfile
    """
    filename = keyfile_formatter(pathtophmmerlog)
    keyfilepath = outpath.joinpath(filename)

    if not pathtophmmerlog.is_file():
        raise FileNotFoundError(f'PHMMERLOG MISSING: File {pathtophmmerlog} not found!')
    elif pathtophmmerlog.stat().st_size == 0:
        raise ValueError(f'PHMMERLOG EMPTY: File {pathtophmmerlog} contains nothing!')
    else:
        acclist = get_accidlist(pathtophmmerlog)
        if not acclist:
            raise ValueError(f'NO HITS: No hits above incl. thresh found in {pathtophmmerlog}')
        if not keyfilepath.is_file():
            writeout_list(acclist, keyfilepath) 
        elif overwrite:
            writeout_list(acclist, keyfilepath) 
        else:
            print(f'Keyfile: {keyfilepath} exists already.\nOverwrite by passing --redo True.')
        return keyfilepath
