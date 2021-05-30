#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-5-30
# DESCRIPTION :  Identify CRISPR arrays that may have undergone recombination with another array in the same genome.
# Approach: for a pair of arrays (A and B) to be identified as possible recombinants they should:
#	Be in the same genome
#	Share no spacers
#	Be Share spacers with 2 arrays (C and D)
#	C and D should be in a different genome
#	C and D should share no spacers with one another
# This is looking for the following pattern (where letters are array names and numbers are spacer IDs):
#
#	A	1 2 3 4
#	B	5 6 7 8
#
#	C	1 2 7 8
#	D	5 6 3 4

import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
	description="Identify CRISPR arrays that may have undergone recombination with another array in the same genome.")
parser.add_argument(
	"-n", dest="network", required = True,
	help="network file produced by CRISPRspacers2network.py"
	)
parser.add_argument(
	"-a", dest="array_reps", required = True,
	help="Array representatives file output by CRISPRspacers2network.py that lists the genomes each array is found in."
	)
parser.add_argument(
	"-i", dest="Array_IDs", required = True,
	help="Array IDs file with 1 array per line, array ID in the first column and spacer IDs or sequences in subsequent columns."
	)


args = parser.parse_args(sys.argv[1:])


network_edges = defaultdict(list)


with open(args.network, 'r') as fin:
	for line in fin.readlines()[1:]:
		array1 = line.split()[0]
		array2 = line.split()[2]
		network_edges[array1].append(array2)
		network_edges[array2].append(array1)

# print(list(network_edges.items())[:1])
