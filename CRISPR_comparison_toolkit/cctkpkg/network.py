#!/usr/bin/env python3

import sys
import argparse

from . import sequence_operations, file_handling

description = """
usage: cctk network [-h] -i <path> [-t <path>] [-o <path>] [--min-shared <int>]

required argument:
  -i, --input 	Array_IDs.txt or Array_seqs.txt

optional arguments:
  -h, --help  	show this help message and exit
  -t, --types 	array CRISPR subtypes file. 2 columns: Array Type
  -o, --outdir	output directory path. Default ./
  --min-shared  minimum spacers shared to draw an edge in network
"""

def build_parser(parser):
	parser.add_argument(
		"-i", "--input",
		metavar=" ",
		required = True,
		help="Array_IDs.txt or Array_seqs.txt"
		)
	parser.add_argument(
		"-t", "--types",
		metavar=" ",
		help="what CRISPR subtype is each array. 2 columns: Array\tType"
		)
	parser.add_argument(
		"-o", "--outdir",
		metavar=" ",
		default="./",
		help="output directory path"
		)
	parser.add_argument(
		"--min-shared",
		metavar=" ",
		default=1,
		type=int,
		help="Minimum number of spacers shared to draw an edge"
		)
	

	return parser


def main(args):

	if args.outdir[-1] != "/":
		args.outdir+= "/"

	array_dict = file_handling.read_array_file(args.input)

	if args.types:
		array_types = file_handling.read_array_types_file(args.types)

	array_list = []
	# convert to FoundArray class for downstream use
	for array_id, spacers in array_dict.items():
		array = file_handling.FoundArray()
		array.id = array_id
		array.spacers = spacers
		if args.types:
			array.repeat_id = array_types[array.id]
		array_list.append(array)
			

	# Build network of array spacer sharing
	network = sequence_operations.build_network(array_list)

	# Write network file
	file_handling.write_network_file(
		network,
		args.outdir+"Array_network.txt",
		args.min_shared
		)

	# Identify network clusters for cluster file
	network_filt = [
		(edge.a, edge.b) for edge in network if edge.nshared >= args.min_shared
		]

	clusters = sequence_operations.identify_network_clusters(network_filt)

	# need to unpack clustes from dict for file
	cluster_list = [[k for k in clus.keys()] for clus in clusters]

	file_handling.write_array_cluster_file(
		cluster_list,
		args.outdir+"Array_clusters.txt"
		)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Create a network representation of the relationships \
		between arrays."
		)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)	
