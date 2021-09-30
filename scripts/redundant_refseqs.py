"""Python script to select unique reference sequences.
Usage: python redundant_refseqs.py refseqfiles.txt ../refseqfiles"""

import argparse


def argument_parsing():
    """Parses input arguments and gives help message"""
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] list_refseqfiles path_refseqfiles",
                                     description="Analyses list of reference sequences and marks uniques as included."
                                                 "Redundant sequences are excluded")
    parser.add_argument("list_refseqfiles", help="A list of reference sequence files in every new line", type=str)
    parser.add_argument("path_refseqfiles", help="path where reference sequences are stored", type=str)
    return parser.parse_args()


def first_line_of_file(path, filename):
    """
    return the first line of a file, remove ending line break.
    :param path: path where file is stored (str)
    :param filename: filename (str)
    """
    filepath = "/".join([path, filename])
    with open(filepath) as f:
        header = f.readline().replace("\n", "")
    f.close()
    return header


def write_pdbid2file(pdbid, filename):
    """Adds a four-digit PDB ID in a new line at the end of a file
    :param pdbid: PDB ID (str)
    :param filename: name of file where PDB ID will be appended
    """
    f = open(filename, "a")
    f.write(pdbid)
    f.write("\n")
    f.close()


def main():
    args = argument_parsing()
    list_refseqfiles = args.list_refseqfiles
    path_refseqfiles = args.path_refseqfiles

    incfile = "included_pdbs_by_refseq.txt"
    excfile = "excluded_pdbs_by_refseq.txt"

    count = 0
    header_set = set()  # set that stores both fasta headers of protein dimers as tuple: (header1, header2)

    # iterate over reference sequence files
    for line in open(list_refseqfiles, "r"):

        file = line.replace("\n", "")  # remove linebreak at end of line
        pdbid = file.split("_")[0]
        count += 1

        if count % 2 == 1:  # count is odd: chain 1
            # open reference sequence file and read header (first line)
            header1 = first_line_of_file(path_refseqfiles, file)

        else:  # count is even: chain 2
            header2 = first_line_of_file(path_refseqfiles, file)

            header_tuple = (header1, header2)

            if header_tuple in header_set:
                # Header combination already appeared. Reference sequence is redundant. Exclude PDB entry
                write_pdbid2file(pdbid, excfile)
            else:
                # this header combination appears for the first time. Include PDB ID
                write_pdbid2file(pdbid, incfile)
                # add header tuple to set, cause it's not represented there yet
                header_set.add(header_tuple)


if __name__ == "__main__":
    main()
