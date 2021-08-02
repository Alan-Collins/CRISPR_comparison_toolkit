#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-8-2
# DESCRIPTION :  Generate random bifircating tree topologies and calculate parsimony of each then return the best.

import sys
import argparse
import CRISPR_mp


parser = argparse.ArgumentParser(
	description="Generate random bifircating tree topologies and calculate parsimony of each then return the best.")
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"-r",  dest="replicates", type=int, nargs="?", default = 1,
		help="Specify number of replicates of tree building to perform. The more replicates, the greater the chance that a better tree will be found. Default: 1"
	)
parser.add_argument(
	"arrays_to_join", nargs="*",  
	help="Specify the IDs of the arrays you want to join. If none provided, joins all arrays in the provided array representatives file. **If given, must come at the end of your command after all other arguments.**"
	)

args = parser.parse_args(sys.argv[1:])
