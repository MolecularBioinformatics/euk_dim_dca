#!bin/bash
# input1: list of filenames that contain single fasta sequences.
# inout2: path where fasta files are stored
# input3: filename of output file
# output: multi-fasta file that contains joint fasta sequences.
# Sequences are joined pairwise (1st with 2nd, 3rd with 4th ...)

list=$(cat $1)
i=1
path=$2
outfile=$3

# loop over filenames
for line in $list
do
	# iteration count $i is even or odd -> current filename ($line) belongs to chain 1 or 2
	if [ $((i%2)) -eq 0 ];
	then  # chain 2
		# read header and sequence of fasta file
		header2=$(grep "^>" $path$line)
                seq2=$(grep "^>" -A1 -v $path$line)
		
		# write output file
		echo "${header1}||${header2}" >> $outfile
		echo "$seq1$seq2" >> $outfile

		
	else  # chain 1
		# read header and sequence
		header1=$(grep "^>" $path$line)
		seq1=$(grep "^>" -A1 -v $path$line)
	fi

	i=$((i+1))
done
