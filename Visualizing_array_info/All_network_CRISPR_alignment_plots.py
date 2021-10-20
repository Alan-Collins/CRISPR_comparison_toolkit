#!/usr/bin/env python3

import sys
import argparse
from itertools import combinations
import pickle
import subprocess


parser = argparse.ArgumentParser(
	description="Given an array file, runs CRISPR_alignment_plot.py for all the clusters of arrays sharing at least N spacers. N can be set using the -n option.")
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
# parser.add_argument(
# 	"-i", dest="network_file", required = False,
# 	help="Specify array network file if you already have a network representation of your arrays.\nFormat: 3 columns. Node1\tNode2\t#_shared_spacers"
# 	)
parser.add_argument(
	"-o", dest="out_dir", required = True, 
	help="Specify output directory for plots."
	)
parser.add_argument(
	"-on", dest="network_out_file", required = False,
	help="Specify array network file name to save a network representation of your arrays."
	)
parser.add_argument(
	"-p", dest="plot_py_path", required = False, 
	help="Specify the path to CRISPR_alignment_plot.py if it isn't in your system PATH. Default behaviour is to run CRISPR_alignment_plot.py as if it is executable and in your PATH."
	)
parser.add_argument(
	"-n", dest="n", nargs="?", default = 2, type=int,
	help="(Default = 2) Set the minimum number of spacers that need to be shared between two arrays in order for there to be an edge between them in the network."
	)
parser.add_argument(
	"-m", dest="min_nodes", nargs="?", default = 2, type=int,
	help="(Default = 2) Set the minimum number of arrays that need to be in a group in order to bother plotting their alignment."
	)
parser.add_argument(
	"-c", dest="colour_file", required = False, 
	help="Specify file with custom colour list (Optional). Colours must be hex codes. One colour per line with no header line in file. e.g. #fd5925 or rbg(10,100,5)."
	)
parser.add_argument(
	"-iter", dest="iterations", nargs="?", default = 10, type=int,
	help="(Default = n^2, where n is the number of arrays to be aligned) If you are aligning fewer than 9 arrays, a good order will be found using a repeated shuffling search method. Set the number of replicates of this search you want performed. Higher numbers more likely to find the best possible ordering, but take longer."
	)
parser.add_argument(
	"-l", dest="legend", action='store_true',  
	help="Include a legend in the output plot (Highly recommended to use spacer IDs rather than sequences with this setting). N.B. Spacer order in the legend is the same as the order of first instance of spacers working from bottom to top, right to left along your plotted arrays."
	)


def build_network(dct, min_edge, outfile):
	parent_dict = {} # dict where the key is an array and the value is its parent array. Whenever an array gets absorbed into an existing network, the array it shared spacers with becomes the parent array. Following the chain of parents should terminate at the array that is the namesake of the cluster in network_dict.
	network_dict = {} # dict where the key is a namesake array for the cluster and value is a list of arrays in that cluster.
	print("processing {} array comparisons".format(len(list(combinations(dct.keys(), 2)))))
	network_list = []
	for combo in combinations(dct.keys(), 2):
		array1 = combo[0]
		array2 = combo[1]
		spacers1 = dct[array1]
		spacers2 = dct[array2]
		nshared = len(set(spacers1).intersection(spacers2))
		if nshared >= min_edge:
			network_list.append('\t'.join([array1, array2, str(nshared), str(nshared/len(set(spacers1 + spacers2)))]))
			if array1 in parent_dict.keys(): # If array is already in a cluster, follow the chain of parents to the namesake at the top.
				parent1 = parent_dict[array1]
				while parent1 in parent_dict.keys(): # Follow the parents up the chain to find the top of the pile.
					parent1 = parent_dict[parent1]
			else:
				parent1 = array1
			if array2 in parent_dict.keys():
				parent2 = parent_dict[array2]
				while parent2 in parent_dict.keys():
					parent2 = parent_dict[parent2]
			else:
				parent2 = array2
			if parent1 != parent2:
				if parent1 in network_dict.keys() and parent2 in network_dict.keys():
					#print("Combining cluster {} and cluster {} which have members:\n{}\n{}".format(parent1, parent2,network_dict[parent1], network_dict[parent2]))
					network_dict[parent1] = list(set(network_dict[parent1] + network_dict[parent2])) # Add the list of cluster 2 members to cluster 1
					del network_dict[parent2] # Remove the now absorbed cluster 2 from the dict
					parent_dict[array2] = parent1
				elif parent1 in network_dict.keys():
					network_dict[parent1].append(array2) # If parent2 not in network_dict.keys() then parent2 must be the same as array2 so no need to worry about adding a whole chain to the existing array
					parent_dict[array2] = parent1
				elif parent2 in network_dict.keys():
					network_dict[parent2].append(array1)
					parent_dict[array1] = parent2
				else:
					network_dict[parent1] = [parent1, parent2] # If neither are in the network_dict, then create a new cluster with just the 2 arrays in it
					parent_dict[array2] = parent1
	for k,v in network_dict.items():
		network_dict[k] = list(set(v))

	return network_dict, network_list

		


args = parser.parse_args(sys.argv[1:])

infile = args.array_file
outdir = args.out_dir + '/' if args.out_dir[-1] != '/' else args.out_dir



array_dict = {}
with open(infile, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

array_network_dict, array_network_list = build_network(array_dict, args.n, args.network_out_file)


if args.network_out_file:
	with open(args.network_out_file + '_dict.pkl', 'wb') as pklout:
		pickle.dump(array_network_dict, pklout)

	with open(args.network_out_file + "_clusters.txt", 'w+') as fout:
		fout.write('\n'.join([k + '\t' + ' '.join(v) for k,v in array_network_dict.items()]))

	with open(args.network_out_file + "_network.txt", 'w+') as fout:
		fout.write('\n'.join(array_network_list))

if args.plot_py_path:
	plot_path = "python3 " + args.plot_py_path
else:
	plot_path = "CRISPR_alignment_plot.py"

for k,v in list(array_network_dict.items()):
	if len(v) > args.min_nodes:
		command_components_list = [plot_path, "-a", args.array_file, "-o", outdir + k + "_cluster.png"]
		if args.iterations:
			command_components_list += ["-iter", str(args.iterations)]
		else:
			command_components_list += ["-iter", str(len(v)**2)]
		if args.legend:
			command_components_list += ["-l"]
		if args.colour_file:
			command_components_list += ["-c", args.colour_file]
		command_components_list += v
		command = " ".join(command_components_list)

		print("\nAligning the {} arrays in cluster {}\nCommand: {}\n".format(len(v), k, command))

		run = subprocess.Popen(command, shell=True)
		run.wait()
