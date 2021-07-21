#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021/7/21
# DESCRIPTION :  Run CRISPR_mp.py on all clusters in network of CRISPR arrays.

import sys
import argparse
import subprocess
from itertools import combinations


parser = argparse.ArgumentParser(
	description="")
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"-o", dest="out_dir", required = True, 
	help="Specify output directory for plots."
	)
parser.add_argument(
	"-r",  dest="replicates", type=int, nargs="?", default = 1,
	help="Specify number of replicates of tree building to perform. The more replicates, the greater the chance that a better tree will be found. Default: 1"
	)
parser.add_argument(
	"-t",  dest="num_threads", type=int, nargs="?", default = 1,
		help="Specify number of threads to use for building trees. Using multiple threads will speed up the search for trees when performing many replicates with the -r option. Default: 1"
	)
parser.add_argument(
	"-m",  dest="min_shared", type=int, nargs="?", default = 1,
		help="Specify Minimum number of spacers that must be shared between two arrays for them to be joined into a cluster to be plotted. Default: 1"
	)
parser.add_argument(
	"-p", dest="CRISPR_mp_py_path", required = False, default = "CRISPR_mp.py",
	help="Specify the path to CRISPR_mp.py if it isn't in your system PATH. Default behaviour is to run CRISPR_mp.py as if it were executable and in your PATH."
	)
parser.add_argument(
	"-x", dest="extra_args", required = False, default = "",
	help="If you want to provide any other arguments that are accepted by CRISPR_mp.py then provide them as a quoted string here."
	)

args = parser.parse_args(sys.argv[1:])

def build_network(dct, min_edge):
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


infile = args.array_file
outdir = args.out_dir + '/' if args.out_dir[-1] != '/' else args.out_dir



array_dict = {}
with open(infile, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

array_network_dict, array_network_list = build_network(array_dict, args.min_shared)

for k,v in list(array_network_dict.items()):
	if len(v) > 2:
		out_image = "{}{}_cluster_tree.png".format(outdir,k)
		out_file = "{}{}_cluster.log".format(outdir,k)

		command = "python3 {} -a {} -r {} -o {} -t {} {} {} > {}".format(args.CRISPR_mp_py_path, infile, args.replicates, out_image, args.num_threads, args.extra_args, " ".join(v), out_file)

		print("\nAligning the {} arrays in cluster {}\nCommand: {}\n".format(len(v), k, command))

		run = subprocess.Popen(command, shell=True)
		run.wait()
