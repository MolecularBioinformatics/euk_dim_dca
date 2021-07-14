#!/usr/bin/env python3
"""
process_phmmerhits.py

Takes in a phmmer-parsed-out keyfile.

Checks if more than 100 seq hits are found,
puts a cap of total number of sequences at 600,
and matches organisms.

Expects a uniprot-like FASTA accession header.

' >db|UniqueIdentifier|EntryName ProteinName ... '

Where the organism is the 5-letter tag at the end of EntryName.
"""

from pathlib import Path

import io_utils as io


def get_two_keyfiles(pathtophmmer, pdbid):
    """Checks dir for exactly 2 phmmer keyfiles for pdbid.
    Returns list of these two keyfiles if found.

    :param pathtophmmer: pathlib.PosixPath
    :param pdbid: str, 4-letter PDB id

    :returns keyfiles: list of pathlib.PosixPaths
    """

    keyfilefmt = f'{pdbid}*_refseq_phmmer.keyfile' 
    keyfiles = io.get_globbed_list(pathtophmmer, keyfilefmt)

    if len(keyfiles) != 2:
        raise ValueError(f'Number of keyfiles {len(keyfiles)} not equal to 2:\n{keyfiles}') 
    return keyfiles


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
    
    :returns orgset: set of organisms
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
    orgset = set(orgdict.keys()) 
    return orgset, orgdict


def select_seqheader_from_org(orgset, orgheaderdict):
    """Selects for each org a fasta header for
    seq from that org. Prioritizes swissprot proteins.

    :param orgset: set of organism tags
    :param orgheaderdict: dict of {org: header, ...}

    :returns: set of headers
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
    return set(headerlist) 


def match_orgtags(orgset, *orgsets):
    """Takes in minimum 1 set of orgtags.
    Returns a set with intersection of orgsets
    
    :param orgset: set of organism tags
    
    :returns: set of common organisms"""
    if not orgsets:
        return orgset
    else:
        return orgset.intersection(*orgsets)


def process_phmmerhits(pathtophmmer, pdbid, minhits, maxhits, redo=False):
    """Performs post-processing of phmmer hits.
    Checks for suitable number of hits returned.
    Matches organisms for the hits and returns
    a matched keyfile.

    :param pathtophmmer: pathlib.PosixPath
    :param pdbid: str, 4-letter PDB id 
    :param minhits: int, minimum number of hits
    :param maxhits: int, maximum number of hits

    :returns keylist: list of fasta seq ids
    """
    print(f'MINHITS: {minhits}')
    print(f'MAXHITS: {maxhits}')
    print(f'OVERWRITE: {redo}\n')

    keyfilepaths = get_two_keyfiles(pathtophmmer, pdbid) 

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

    if len(masterorgset) >= maxhits: # TODO list of setsare still unordered
        print(f'ATTENTION: Number of seqs for {pdbid} {len(masterorgset)} exceeded limit of {maxhits}!')
        print(f'           Matched keyfiles are reduced to {maxhits} sequences.')
        masterlist = list(masterorgset)
        masterorgset = set(masterlist[0:maxhits])

    for keyfile, entry in hits.items():
        hits[keyfile] = select_seqheader_from_org(masterorgset, entry[1])
    

    keylist = []
    for keyfile, entry in hits.items(): # TODO do you need two for loops
        outfile = io.matched_keyfile_formatter(keyfile)
        outpath = pathtophmmer / outfile
        keylist.append(outpath)
        io.writeout_list(list(entry), outpath) # TODO check if list is needed

    return keylist[0], keylist[1]
