#!/usr/bin/env python3
"""
process_alnseqs.py

Matches aligned sequences by organism,
joins the matched sequences.

Prepares an alignment for DCA.
"""

from pathlib import Path

from io_utils import does_target_exist, parse_fasta, fa_todict, writeout_fasta
from pydca.msa_trimmer import msa_trimmer


def move_refseq_up(aln_path):
    """Moves reference sequence up to first
    sequence in an alignment file. This is done
    to force pydca's trim-by-refseq to choose
    the correct seq as the reference (in case multiple
    copies of the refseq exist in the alignment file)

    Searches for RFSEQ tag. Prints out another fasta
    file with RFSEQ first.

    :param aln_path: pathlib.PosixPath"""
    faorgdict = get_orgdict_from_fafile(aln_path)

    refseq_entry = faorgdict['RFSEQ'] 
    fadict = {}
    fadict[refseq_entry[0]] = refseq_entry[1]

    for org in faorgdict.keys():
        if org != 'RFSEQ':
            entry = faorgdict[org]
            fadict[entry[0]]=entry[1]

    writeout_fasta(aln_path, fadict, overwrite=True)


def trim_msa_by_refseq(aln_path, refseqpath):
    """Trims an alignment based on the reference
    sequence. Writes out trimmed alignment to file.

    :param aln_path: pathlib.PosixPath
    :param refseqpath: pathlib.PosixPath"""
    
    trim_outpath = Path(f"{aln_path.stem}_trimmed.fasta")

    trimmer = msa_trimmer.MSATrimmer(aln_path, 
                                     biomolecule='protein',
                                     refseq_file=refseqpath)

    trim_data = trimmer.get_msa_trimmed_by_refseq(remove_all_gaps=True)

    return trim_data


def trim_list_to_fadict(trimlist):
    """Converts a list of tuples [('header','seq'), ...] to
    a dictionary of {'header':'seq', ... }.

    :param trimlist: list, from pydca's trimmer"""
    fadict = {}
    for pair in trimlist:
        fadict[pair[0]] = pair[1]
    return fadict


def get_orgdict_from_fadict(fadict):
    """Takes in a fasta dict of form 
    {'header':'seq', ... }.
    Returns a dict of format {'ORG':('header','seq'),...}

    :param fastadict: dict

    :returns: dict"""

    fa_orgdict = {}
    for header in fadict.keys():
        firstpart = header.strip().split()[0]
        orgtag = firstpart.split("_")[-1:][0]
        fa_orgdict[orgtag] = (header, fadict[header])
    return fa_orgdict


def get_orgdict_from_fafile(fafilepath):
    """
    Takes in fasta file (with uniprot-style headers),
    returns a dict of format {'ORG':('header', 'seq'),...}

    :param fafilepath: pathlib.PosixPath
    :returns fa_orgdict: dict"""

    fadict = fa_todict(fafilepath)
    fa_orgdict = {}
    for header in fadict.keys():
        firstpart = header.strip().split()[0]
        orgtag = firstpart.split("_")[-1:][0]
        fa_orgdict[orgtag] = (header, fadict[header])
    return fa_orgdict


def join_two_orgdicts(orgdict1, orgdict2):
    """Takes in two dicts of format {'ORG':('header', 'seq'),...},
    Joins the headers with '|' in between, and
    joins the sequences. Makes a joint file for DCA.

    :param orgdict: dict
    :param dcadict: dict, joined headers and aligned seqs"""

    dcadict = {}
    
    if len(orgdict1) != len(orgdict2):
        raise ValueError('Fastas you are trying to join have unequal # seqs!')
    for orgtag in orgdict1.keys():
        header1 = orgdict1[orgtag][0]
        header2 = orgdict2[orgtag][0]
        seq1 = orgdict1[orgtag][1]
        seq2 = orgdict2[orgtag][1]
        dcadict[header1+'||'+header2] = seq1+seq2
    return dcadict


def process_alnseqs(alnfile1_path, alnfile2_path, refseq1_path, refseq2_path, phmmerpath, redo):
    """Prepares aligned sequences from two alignment files for 
    DCA. First matches the sequences based on organism, then 
    zips them up (joins) them into one joint alignment file.

    :param alnfile_path: pathlib.PosixPath, 1 or 2 being either of the fastas
    :param phmmerpath: pathlib.PosixPath
    :param redo: bool

    :returns jointalnfile_path: pathlib.PosixPath"""

    ch1 = alnfile1_path.stem.split('_refseq')[0]
    ch2 = alnfile2_path.stem.split('_refseq')[0]
    outpath = phmmerpath / f"Joint_{ch1}_{ch2}_aln.fasta"

    if not (does_target_exist(alnfile1_path, 'file') and does_target_exist(alnfile2_path, 'file')):
        raise FileNotFoundError('Check that your alignment files exist!')
    if redo == False and does_target_exist(outpath, 'file'):
        print(f'Joint alignment already exists: {outpath}')
        return outpath

    move_refseq_up(alnfile1_path)
    move_refseq_up(alnfile2_path)

    print(f'trimming msa 1 ...')
    trim1list = trim_msa_by_refseq(alnfile1_path, refseq1_path)
    print(f'trimming msa 2 ...')
    trim2list = trim_msa_by_refseq(alnfile2_path, refseq2_path)

    fadict1 = trim_list_to_fadict(trim1list)
    fadict2 = trim_list_to_fadict(trim2list)

    orgdict1 = get_orgdict_from_fadict(fadict1)
    orgdict2 = get_orgdict_from_fadict(fadict2)

    print(f'joining trimmed msas...')
    jointdict = join_two_orgdicts(orgdict1, orgdict2)

    writeout_fasta(outpath, jointdict, overwrite=True)
    print(f'Joint alignment written into: {outpath}')

    return outpath
