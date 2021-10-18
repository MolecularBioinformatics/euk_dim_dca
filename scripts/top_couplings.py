#!/usr/bin/env python
"""
Usage:  {prog} matrix_infile score_outfile (num-contacts) (min-separation)

matrix_infile: DCA scores as matrix (.mat), output from CCMpred
score_outfile: name of output file, which will be list of top scoring residue pairs
num-contacts: Set the number of pairs to output [default: 30]
min-separation: Set the minimum sequence separation of pairs to be outputted [default: 7]

copied (1.10.2021) and modified (15.10.2021) from https://github.com/soedinglab/CCMpred
released under liscense: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
"""

import numpy as np


def main(args):

    # check number of input arguments. Has to be between 2 and 4
    if len(args) < 2:
        raise ValueError("Too few arguments. Need at least matrix input file and output filename")
    elif len(args) > 4:
        raise ValueError("Too many arguments. Maximum: 4 arguments")

    # set default values
    num_contacts = 30
    min_separation = 7

    # parse args
    matrix_infile, outfile = args[:2]
    if len(args) == 3:
        num_contacts = int(args[2])
    elif len(args) == 4:
        num_contacts = int(args[2])
        min_separation = int(args[3])

    # load coupling matrix
    mat = np.loadtxt(matrix_infile)

    # find top-scoring pairs with sufficient separation
    top = get_top_pairs(mat, num_contacts, min_separation)

    file = open(outfile, "w")
    for i, j, coupling in zip(top[0], top[1], mat[top]):
        file.write(f"{i}\t{j}\t{coupling}\n")
    file.close()


def get_top_pairs(mat, num_contacts, min_separation):
    """Get the top-scoring contacts
    :param mat: numpy matrix of residue pair scores
    :param num_contacts: number of top residue pairs that should be extracted
    :param min_separation: minimum sequence separation of residue pairs
    """

    idx_delta = np.arange(mat.shape[1])[np.newaxis, :] - np.arange(mat.shape[0])[:, np.newaxis]
    mask = idx_delta < min_separation

    mat_masked = np.copy(mat)
    mat_masked[mask] = float("-inf")

    top = mat_masked.argsort(axis=None)[::-1][:(num_contacts)]
    top = (top % mat.shape[0]).astype(np.uint16), np.floor(top / mat.shape[0]).astype(np.uint16)
    return top


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
