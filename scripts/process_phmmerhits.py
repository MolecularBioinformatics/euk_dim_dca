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

from io_utils import get_globbed_list
from io_utils import readin_list, writeout_list, matched_keyfile_formatter


def check_two_keyfiles(pathtophmmer, pdbid):
    """Checks dir for exactly 2 phmmer keyfiles for pdbid

    :param pathtophmmer: pathlib.PosixPath
    :param pdbid: str, 4-letter PDB id

    :returns keyfiles: list of pathlib.PosixPaths
    """

    keyfilefmt = f'{pdbid}*_refseq_phmmer.keyfile' 
    keyfiles = get_globbed_list(pathtophmmer, keyfilefmt)

    if len(keyfiles) != 2:
        raise Exception(f'Number of keyfiles {len(keyfiles)} not equal to 2:\n{keyfiles}')
    else:
        return keyfiles


def check_hit_nr(pathtokeyfile, hitthresh, relation):
    """Checks that num of hits in keyfile is higher
    than minhitnum and lower than maxhitnum

    :returns hits: list
    """

    hits = readin_list(pathtokeyfile)
    hitnum = len(hits)
    
    if relation == 'MINIMUM':
        if hitnum <= hitthresh:
            raise Exception(f'WARNING: {pathtokeyfile.name}: Too FEW hits returned from phmmer: {len(hits)}')
        else:
            return hits

    elif relation == 'MAXIMUM':
        if hitnum >= hitthresh:
            raise Exception(f'Too MANY hits returned from phmmer: {len(hits)}')
        else:
            return hits


def get_orgs_from_hitlist(hitlist):
    """Takes list hits, returns set of all orgs
    and dict with headers from each organism"""

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
    Returns a set with intersection of orgsets"""
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
    """
    print(f'MINHITS: {minhits}')
    print(f'MAXHITS: {maxhits}')
    print(f'OVERWRITE: {redo}\n')

    try:
        keyfilepaths = check_two_keyfiles(pathtophmmer, pdbid)
    except Exception as e:
        print(e)

    if keyfilepaths:
        hits = {}
        for keyfile in keyfilepaths:
            try:
                hitlist = check_hit_nr(keyfile, minhits, 'MINIMUM') 
                hits[keyfile] = hitlist
            except Exception as e:
                print(e)
                break
        if hits:
            keylist = []
            for keyfile, hitlist in hits.items():
                orgset, orgdict = get_orgs_from_hitlist(hitlist) 
                hits[keyfile] = (orgset, orgdict)

            orgsets = [entry[0] for entry in hits.values()]
            masterorgset = match_orgtags(*orgsets)

            if len(masterorgset) >= maxhits:
                print(f'ATTENTION: Number of seqs for {pdbid} {len(masterorgset)} exceeded limit of {maxhits}!')
                print(f'           Matched keyfiles are reduced to {maxhits} sequences.')
                masterlist = list(masterorgset)
                masterorgset = set(masterlist[0:maxhits])

            for keyfile, entry in hits.items():
                headers = select_seqheader_from_org(masterorgset, entry[1])
                hits[keyfile] = headers

            for keyfile, entry in hits.items():
                outfile = matched_keyfile_formatter(keyfile)
                outpath = pathtophmmer.joinpath(outfile)
                keylist.append(outpath)
                writeout_list(list(entry), outpath)

            return keylist[0], keylist[1]
