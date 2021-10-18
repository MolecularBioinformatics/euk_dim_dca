#!/usr/bin/env python
"""
Usage:  {prog} infile format outfile
        {prog} infile format
        {prog} format

copied (1.10.2021) and modified (7.10.2021) from https://github.com/soedinglab/CCMpred
released under liscense: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
"""

import Bio.SeqIO


def convert(f_in, format, f_out):
    for record in Bio.SeqIO.parse(f_in, format):
        try:
            f_out.write(record.seq.tostring())
        except AttributeError:
            f_out.write(str(record.seq))
        f_out.write("\n")


def main(args):
    import sys

    fn_in = "-"
    fn_out = "-"

    if len(args) == 3:
        fn_in, format, fn_out = args
    elif len(args) == 2:
        fn_in, format = args
    elif len(args) == 1:
        format, = args
    else:
        sys.stderr.write("Need 1-3 arguments!\n")
        usage()
        sys.exit(1)

    if fn_in == "-":
        f_in = sys.stdin
    else:
        f_in = open(fn_in, "r")

    if fn_out == "-":
        f_out = sys.stdout
    else:
        f_out = open(fn_out, "w")

    convert(f_in, format, f_out)

    if f_in != sys.stdin:
        f_in.close()

    if f_out != sys.stdout:
        f_out.close()


def usage():
    import sys
    print(__doc__.format(prog=sys.argv[0]))


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
