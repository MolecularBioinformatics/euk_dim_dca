# Julia script to run gaussDCA
# input arg: MSA file, output file, minimum residue pair separation

using GaussDCA
aln = ARGS[1]
outfile = ARGS[2]
sep = ARGS[3]

# run Gauss DCA on input file
FNR = gDCA(aln, min_separation=sep);

# save output file
printrank(outfile, FNR)
