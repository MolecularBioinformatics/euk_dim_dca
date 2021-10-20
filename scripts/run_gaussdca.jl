# Julia script to run gaussDCA
# input arg: MSA file, output file, minimum residue pair separation

using GaussDCA
aln = ARGS[1]
outfile = ARGS[2]
sep = ARGS[3]
sep = parse(Int64, sep)  # string to int

# run Gauss DCA on input file
FNR = gDCA(aln, min_separation=sep);

# save output file
printrank(outfile, FNR)