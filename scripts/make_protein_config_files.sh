#!/bin/bash

# script to create empty protein config files (config_XXXX.txt) to run the pydca_workflow.py on
# takes two input arguments
# 1) default config file config_template.txt 
# 2) list of PDB IDs in every new line

# empty protein config file, usually: config_template.txt
config_tmp=$1

# list of PDB IDs
list=$(cat $2)

for line in $list
do
	pdb_id="$(echo $line | tr '[A-Z]' '[a-z]')"  # PDB ID to lowercase

	cat $config_tmp | sed s/pdbid=/pdbid=$pdb_id/ > config_"$pdb_id".txt
done
