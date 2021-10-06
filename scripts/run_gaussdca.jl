# Julia script to run gaussDCA
# input arg: MSA file

using GaussDCA
aln = ARGS[1]
outfile = ARGS[2]

# run Gauss DCA on input file
FNR = gDCA(aln);

# save output file
printrank(outfile, FNR)
