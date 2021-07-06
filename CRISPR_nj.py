#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-7-6
# DESCRIPTION :  Perform neighbour joining on a given group of CRISPR arrays

import sys
import argparse



parser = argparse.ArgumentParser(
	description="Perform neighbour joining on a given group of CRISPR arrays"
	)
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"arrays_to_join", nargs="+",  
	help="Specify the IDs of the arrays you want to join. **Must come at the end of your command after all other arguments.**"
	)



args = parser.parse_args(sys.argv[1:])


array_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

