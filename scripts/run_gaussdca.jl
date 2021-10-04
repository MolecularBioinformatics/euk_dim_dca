using GaussDCA
aln = ARGS[1]
# run Gauss DCA on input file
FNR = gDCA(aln);
outfile = aln[7:19]*"_gaussdca_scores.txt"
# save output file
printrank(outfile, FNR)
