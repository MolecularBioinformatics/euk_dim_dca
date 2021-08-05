#!/usr/bin/env python3
"""
process_phmmerhits.py

Takes in a phmmer-parsed-out keyfile.

Checks that a minimum number (100) of seq hits is found,
caps the total number of seq hits (600).

Matches organisms and returns 

Checks if more than 100 seq hits are found,
puts a cap of total number of sequences at 600,
and matches organisms.

Expects a uniprot-like FASTA accession header.

' >db|UniqueIdentifier|EntryName ProteinName ... '

Where the organism is the 5-letter tag at the end of EntryName.
"""

from pathlib import Path

import io_utils as io
from ordered_set import OrderedSet 


def get_two_keyfiles_old(pathtokeyfiles, pdbid):
    """Checks dir for exactly 2 phmmer keyfiles for pdbid.
    Returns list of these two keyfiles if found.

    :param pathtokeyfiles: pathlib.PosixPath
    :param pdbid: str, 4-letter PDB id

    :returns keyfiles: list of pathlib.PosixPaths
    """

    keyfilefmt = f'{pdbid}*_refseq_phmmer.keyfile' 
    keyfiles = io.get_globbed_list(pathtokeyfiles, keyfilefmt)

    if len(keyfiles) != 2:
        raise ValueError(f'Number of keyfiles {len(keyfiles)} not equal to 2:\n{keyfiles}') 
    return keyfiles


def get_two_keyfiles(pathtokeyfiles, keyfile1_path, keyfile2_path):
    """Returns list of two keyfiles if both exist
    in dir indicated by pathtokeyfiles

    :param pathtokeyfiles: pathlib.PosixPath
    :param keyfile1: pathlib.PosixPath
    :param keyfile2: pathlib.PosixPath

    :returns keyfiles: list of pathlib.PosixPaths
    """
    firstfile = io.does_target_exist(keyfile1_path, 'file')
    secondfile = io.does_target_exist(keyfile2_path, 'file')

    if not firstfile:
        raise FileNotFoundError(f'KEYFILE MISSING: {keyfile1_path}')
    if not secondfile:
        raise FileNotFoundError(f'KEYFILE MISSING: {keyfile2_path}')
    return [keyfile1_path, keyfile2_path]


def get_hit_list(pathtokeyfile, hitthresh, relation):
    """Checks that num of hits in keyfile is
    between minhitnum and maxhitnum, if so,
    returns list of hits.

    :param pathtokeyfile: pathlib.PosixPath
    :param hitthresh: int, max num of hits allowed
    :param relation: str, 'MINIMUM' or 'MAXIMUM'

    :returns hits: list of fasta seq ids
    """

    hits = io.readin_list(pathtokeyfile)
    hitnum = len(hits)
    
    if relation == 'MINIMUM' and hitnum <= hitthresh:
        raise ValueError(f'WARNING: {pathtokeyfile.name}: Too FEW hits returned from phmmer: {len(hits)}')
    elif relation == 'MAXIMUM' and hitnum >= hitthresh:
        raise ValueError(f'Too MANY hits returned from phmmer: {len(hits)}')
    elif relation not in ('MINIMUM','MAXIMUM'):
        raise ValueError("Relation has to be MINIMUM or MAXIMUM")
    return hits


def get_orgs_from_hitlist(hitlist):
    """Takes list hits, returns set of all orgs
    and dict with headers from each organism.
    
    :param hitlist: list of fasta seq ids
    
    :returns orgset: OrderedSet of organisms
    :returns orgdict: dict with keys as orgtags, 
        values as a list of fasta seq ids 
    """

    orgdict = {}

    for hit in hitlist:
        org = hit.strip().split('_')[-1]
        if not org in orgdict:
            orgdict[org] = [hit]
        else:
            orgdict[org].append(hit) 
    orgset = OrderedSet(orgdict.keys()) 
    return orgset, orgdict


def select_seqheader_from_org(orgset, orgheaderdict):
    """Selects for each org a fasta header for
    seq from that org. Prioritizes swissprot proteins.

    :param orgset: OrderedSet of organism tags
    :param orgheaderdict: dict of {org: header, ...}

    :returns: OrderedSet of headers
    """

    headerlist = []

    for org in orgset:
        spseqid = None
        headers = orgheaderdict[org]
        for idx, header in enumerate(headers):
            if header.startswith('sp'):
                spseqid = idx
                break
        if not spseqid:
            headerlist.append(headers[0])
        else:
            headerlist.append(headers[spseqid])
    return OrderedSet(headerlist) 


def match_orgtags(orgset, *orgsets):
    """Takes in minimum 1 set of orgtags.
    Returns a set with intersection of orgsets
    
    :param orgset: OrderedSet of organism tags
    
    :returns: set of common organisms"""
    if not orgsets:
        return orgset
    else:
        return orgset.intersection(*orgsets) 


def matched_keyfiles_exist(keyfilepaths, phmmerpath):
    """Takes in two keyfile paths. Formats filenames 
    into matched keyfiles and checks if these exist.
    Returns True if they exist, else returns False.
    """
    matchedkeyfilespaths = [phmmerpath / io.matched_keyfile_formatter(keyfile) for keyfile in keyfilepaths]
    print(f'Matched keyfiles: {matchedkeyfilespaths}')
    if matchedkeyfilespaths[0].exists() and matchedkeyfilespaths[1].exists():
        return True
    else:
        return False


def process_phmmerhits(pathtokeyfiles, keyfile1path, keyfile2path, minhits, maxhits, redo=False):  # TODO: refactor and break up into more functions
    """Performs post-processing of phmmer hits.
    Checks for suitable number of hits returned.
    Matches organisms for the hits and returns
    a matched keyfile.

    :param pathtokeyfiles: pathlib.PosixPath
    :param pdbid: str, 4-letter PDB id 
    :param minhits: int, minimum number of hits
    :param maxhits: int, maximum number of hits

    :returns keylist: list of fasta seq ids
    """
    print(f'MINHITS: {minhits}')
    print(f'MAXHITS: {maxhits}')
    print(f'OVERWRITE: {redo}\n')  # TODO: redo not yet incorporated, it does it automatically every time

    keyfilepaths = get_two_keyfiles(pathtokeyfiles, keyfile1path, keyfile2path)

    if redo==False and matched_keyfiles_exist(keyfilepaths, pathtokeyfiles)==True:
        print(f'Matched keyfiles already exist!')
        return pathtokeyfiles/io.matched_keyfile_formatter(keyfilepaths[0]), pathtokeyfiles/io.matched_keyfile_formatter(keyfilepaths[1])

    hits = {}
    for keyfile in keyfilepaths:
        hits[keyfile] = get_hit_list(keyfile, minhits, 'MINIMUM') 
    if len(hits) != 2:
        raise ValueError("Incorrect number of lists of hits.") 
    for keyfile, hitlist in hits.items():
        orgset, orgdict = get_orgs_from_hitlist(hitlist) 
        hits[keyfile] = (orgset, orgdict)

    orgsets = [entry[0] for entry in hits.values()]
    masterorgset = match_orgtags(*orgsets)

    if len(masterorgset) >= maxhits: 
        print(f'ATTENTION: Number of seqs {len(masterorgset)} exceeded limit of {maxhits}!')
        print(f'           Matched keyfiles are reduced to {maxhits} sequences.')
        masterlist = list(masterorgset)
        masterorgset = OrderedSet(masterlist[0:maxhits])

    for keyfile, entry in hits.items(): 
        hits[keyfile] = select_seqheader_from_org(masterorgset, entry[1])
    
    keylist = []
    for keyfile in hits.keys():
        outfile = io.matched_keyfile_formatter(keyfile)
        outpath = pathtokeyfiles / outfile
        keylist.append(outpath)
        io.writeout_list(list(hits[keyfile]), outpath) 

    return keylist[0], keylist[1]
