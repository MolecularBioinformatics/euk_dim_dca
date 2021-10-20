#!/usr/bin/env python3
"""
run_dca.py

Runs dca for a given joint alignment

1. First trim MSA with pydca's trimmer
2. Then run plmdca, mfdca or gaussdca on aln"""

import sys
import time
import subprocess
from io_utils import does_target_exist
from meanfield_dca import meanfield_dca
import convert_alignment
import top_couplings

sys.path.append("/cluster/projects/nn9795k/yin/pydca-master/pydca")


def run_pydca_mfdca(jointaln_path):
    """
    Spawns subprocess to run pydca mfdca.

    :param jointaln_path: pathlib.PosixPath
        - needs to be just the string due to pydca's code
    :param redo: bool

    :returns mfdca_FN_APC: list
    """

    mfdca_inst = meanfield_dca.MeanFieldDCA(str(jointaln_path), 'protein', pseudocount=0.5, seqid=0.8)

    start = time.perf_counter()
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    stop = time.perf_counter()
    print(f'mfDCA ran in {stop - start:0.4f} seconds')

    return mfdca_FN_APC


def writeout_scores(dcalist, outfilepath):
    """
    Writes out dca scores in 3 column file.

    :param dcalist: list of tuples [((i,j),score),...]
    :param outfilepath: pathlib.PosixPath
    """

    with open(outfilepath, 'w') as outf:
        for scorepair in dcalist: #TODO check residue pair seperation here
            outf.write(f'{scorepair[0][0]}\t{scorepair[0][1]}\t{scorepair[1]}\n')


def run_gaussdca(aln_path, outfile_path):
    """
    Runs GaussDCA over the Julia script run_gaussdca.jl on input alignment and stores the scores in the outfilepath
    :param aln_path: path and filename of joint MSA (pathlib.PosixPath)
    :param outfile_path: path and filename of output DCA score file (pathlib.PosixPath)
    """
    start = time.perf_counter()
    cmd = ["julia", "-t",  "6", "run_gaussdca.jl", str(aln_path), outfile_path]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.stdout:
        print(proc.stdout)
    if proc.stderr:
        print(proc.stderr)
    stop = time.perf_counter()
    if proc.returncode != 0:
        raise ValueError(f'GaussDCA run unsuccessful!')
    print(f'gaussDCA ran in {stop - start:0.4f} seconds')


def reformat_scoring_file(outfilepath):
    """
    Replaces whitespaces with tabs in file stored at outfilepath
    :param outfilepath: path and filename of output DCA score file (pathlib.PosixPath)
    """
    # replace all spaces by tabs and save tab-separated file content in new_file_content
    new_file_content = ""
    for line in open(outfilepath, "r"):
        new_file_content += line.replace(" ", "\t")

    # overwrite existing scoring file with tab-separated scores
    new_file = open(outfilepath, "w")
    new_file.write(new_file_content)
    new_file.close()


def run_plmdca(ccmpredpath, aln, outmtx):
    """
    Runs CCMpred via bash on an input joint alignment which was reformated with convert_alignment.py
    Output scores are stored as matrix
    :param ccmpredpath: path to CCMpred installation (pathlib.PosixPath)
    :param aln: joint alignment (pathlib.PosixPath)
    :param outmtx: output scoring matrix (pathlib.PosixPath)
    """
    run_ccmpred_cmd = ccmpredpath / "bin/ccmpred"

    start = time.perf_counter()
    cmd = [run_ccmpred_cmd, "-t", "6", aln, outmtx]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.stdout:
        print(proc.stdout)
    if proc.stderr:
        print(proc.stderr)
    if proc.returncode != 0:
        raise ValueError(f'CCMpred run unsuccessful!')
    stop = time.perf_counter()
    print(f'CCMpred ran in {stop - start:0.4f} seconds')


def run_dca(jointaln_path, outpath, redo, method, **kwargs):
    """
    Runs dca method (default mf) on a joint alignment.

    Writes out scores to a []_scores.dat file.

    :param jointaln_path: pathlib.PosixPath
    :param outpath: pathlib.PosixPath
    :param redo: bool
    :param method: mf/plm/gauss which DCA method should be used (str)
    :param ccmpredpath: path to CCMpred installation (pathlib.PosixPath)

    :returns scorefile_path: pathlib.PosixPath
    """

    outfilename = f'{jointaln_path.stem}_{method}dca_scores.dat'
    outfilepath = outpath / outfilename

    jointaln_conv = None  # two files that are only used for the plmDCA approach
    outmtx_file = None

    if not does_target_exist(jointaln_path, 'file'):
        raise FileNotFoundError(f'JOINT ALN FILE MISSING: Could not find {jointaln_path}')
    elif does_target_exist(outfilepath, 'file') and redo is False:
        print(f'DCA scores files: ({outfilepath}) already exists in {outfilepath.parent}')
        return outfilepath, jointaln_conv, outmtx_file

    print(f"DCA approach: {method}")

    if method == "mf":
        # mean-field DCA approach by pydca
        dcascores = run_pydca_mfdca(jointaln_path)
        if not dcascores:
            raise ValueError('DCA run unsuccessful!')
        writeout_scores(dcascores, outfilepath)

    elif method == "plm":
        # pseudo-likelihood maximization DCA approach by CCMpred

        ccmpredpath = kwargs.get("ccmpredpath", None)
        alnpath = kwargs.get("alnpath", None)

        jointaln_conv = alnpath / f"{jointaln_path.stem}_conv.fasta"
        outmtx_file = outpath / outfilename.replace(".dat", ".mat")

        # convert alignment that it fits as input for CCMpred
        if redo or not does_target_exist(jointaln_conv, 'file'):
            convert_alignment.main([jointaln_path, "fasta", jointaln_conv])
            print(f"Converted alignment written to {jointaln_conv}")
        else:
            print(f"Converted alignment already exists in {jointaln_conv}\nContinue with it.")

        # run CCMpred on converted alignment. Output is scoring matrix
        if redo or not does_target_exist(outmtx_file, 'file'):
            run_plmdca(ccmpredpath, jointaln_conv, outmtx_file)
            print(f"DCA score matrix written to {outmtx_file}")
        else:
            print(f"DCA score matrix already exists in {outmtx_file}\nContinue with it.")

        # convert matrix to residue pair list with a minimum sequence separation of 5 for all residue pairs
        sep = 5
        top_couplings.main([outmtx_file, outfilepath, sep])

    else:
        # gaussian DCA approach by gaussDCA
        run_gaussdca(jointaln_path, outfilepath)

        # edit format of scoring file
        reformat_scoring_file(outfilepath)

    print(f'DCA scores written into {outfilepath}')
    return outfilepath, jointaln_conv, outmtx_file


def run_mfdca_saga(jointaln_path, outpath, redo):
    """
    Runs dca method on a joint alignment.
    Spawns subprocess to call method from the command line.
    Written to be run on saga (errors with msa_numerics and numba)

    :param jointaln_path: pathlib.PosixPath
    :param outpath: pathlib.PosixPath
    :param redo: bool

    :returns outfilepath: pathlib.PosixPath
    """

    outfilename = f'MFDCA_apc_fn_scores_{jointaln_path.stem}.txt'
    outfilepath = outpath / outfilename
    cmdargs = ['python3', 
               'mfdca_main.py',
               'compute_fn',
               'protein',
               f'{jointaln_path}',
               '--apc',
               '--pseudocount',
               '0.5', 
               '--verbose',
               '--output_dir',
               f'{outpath}']

    if not does_target_exist(jointaln_path, 'file'):
        raise FileNotFoundError(f'JOINT ALN FILE MISSING: Could not find {jointaln_path}')
    elif does_target_exist(outfilepath, 'file') and redo == False:
        print(f'DCA scores files: ({outfilepath}) already exists in {outfilepath.parent}')
        return outfilepath

    proc = subprocess.run(cmdargs)
    if proc.returncode != 0:
        raise ValueError(f'DCA run unsuccessful!')
    print(f'DCA scores written into {outfilepath}')
    return outfilepath
