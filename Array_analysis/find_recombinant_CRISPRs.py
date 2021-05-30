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
from itertools import permutations
import subprocess

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
parser.add_argument(
	"-o", dest="outdir", required = True,
	help="Path to location you want output files written."
	)


args = parser.parse_args(sys.argv[1:])

outdir = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir

network_edges = defaultdict(list)


with open(args.network, 'r') as fin:
	for line in fin.readlines()[1:]:
		array1 = line.split()[0]
		array2 = line.split()[2]
		network_edges[array1].append(array2)
		network_edges[array2].append(array1)

# print(list(network_edges.items())[:1])

array_rep_genomes = {}

with open(args.array_reps, 'r') as fin:
	for line in fin.readlines()[1:]:
		array_rep_genomes[line.split()[0]] = line.split()[1:]

# print(list(array_rep_genomes.items())[:5])

possible_hits = []

for array1, array2 in permutations(list(array_rep_genomes.keys()), 2):
	genome_overlap = set(array_rep_genomes[array1]).intersection(set(array_rep_genomes[array2]))
	if len(genome_overlap) > 1:
		if array2 not in network_edges[array1]:
			edge_overlap = set(network_edges[array1]).intersection(set(network_edges[array2]))
			if len(edge_overlap) > 1:
				for hit1, hit2 in permutations(edge_overlap, 2):
					hit_genome_overlap = set(array_rep_genomes[hit1]).intersection(set(array_rep_genomes[hit2]))
					if len(hit_genome_overlap) > 1:
						if hit2 not in network_edges[hit1]:
							with open("{}{}.txt".format(outdir, "_".join([array1, array2,hit1, hit2])), 'w') as fout:
								fout.write("Arrays {} and {} were found in these genomes: {}\n\
									Connected arrays {} and {} were found in genomes {}\n".format(
										array1, array2, " ".join(genome_overlap),
										hit1, hit2, " ".join(hit_genome_overlap)))

							command = "CRISPR_alignment_plot.py -a {} -o {}{}.png {}".format(args.Array_IDs, outdir, "_".join([array1, array2,hit1, hit2]), " ".join([array1, array2,hit1, hit2]))
							run = subprocess.Popen(command, shell=True)
							run.wait()
			


