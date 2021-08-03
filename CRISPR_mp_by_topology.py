#!/usr/bin/env python3

# AUTHOR	  :  ALAN COLLINS
# VERSION	 :  v1
# DATE		:  2021-8-2
# DESCRIPTION :  Generate random bifircating tree topologies and calculate parsimony of each then return the best.

import sys
import argparse
import dendropy
import CRISPR_mp
import random
from math import ceil, factorial
from itertools import permutations
from string import ascii_lowercase



def make_topologies(leaves, root=False, randomize=False):
	"""
	Given a number of leaves for a tree, generate all possible configurations of the rooted bifurcating tree.
	Args:
		num_leaves (list): Leaves to be placed in the tree.
		root (int): The 0 base integer location at which to place the root and generate all tree topologies with that root.
		randomize (bool): Would you like a single random tree to be returned?
	Yields:
		 Newick string of tree or subtree.
	"""
	if len(leaves) == 1: # When all nodes have been joined, yield the product.
		yield leaves[0]
	else:
		if randomize:
			root = random.randint(1,ceil((len(leaves))/2)) # Root the tree at a random point
			for left in make_topologies(leaves[:root], randomize=True): # Root the subtree to the left at a random point
				for right in make_topologies(leaves[root:], randomize=True): # Root the subtree to the right at a random point
					yield "({},{})".format(left, right)
		elif root:
			for left in make_topologies(leaves[:root]): # Make subtrees to the left of (sub)tree root
				for right in make_topologies(leaves[root:]): # Make subtrees to the right of (sub)tree root
					yield "({},{})".format(left, right)
		else:
			for i in range(1, ceil((len(leaves)+1)/2)):
				for left in make_topologies(leaves[:i]): # Make subtrees to the left of (sub)tree root
					for right in make_topologies(leaves[i:]): # Make subtrees to the right of (sub)tree root
						yield "({},{})".format(left, right)


def main():

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

	array_spacers_dict = {}
	with open(args.array_file, 'r') as fin:
		for line in fin.readlines():
			bits = line.split()
			array_spacers_dict[bits[0]] = bits[2:]

	node_ids = ["Int " + i for i in ascii_lowercase]
	if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
		node_ids += ["Int " + "".join(i) for i in product(ascii_lowercase, repeat=(len(args.arrays_to_join)//26)+1)]

	taxon_namespace = dendropy.TaxonNamespace(args.arrays_to_join + node_ids)
	
	if len(args.arrays_to_join) < 7:
		array_orders = [i for i in permutations(args.arrays_to_join)]
		trees = []
		for order in array_orders:
			trees += [t + ";" for t in make_topologies(order)]
	else:
		array_orders = [random.sample(args.arrays_to_join, len(args.arrays_to_join)) for i in range(args.replicates)]
		trees = []
		for i in range(args.replicates):
			trees.append(next(make_topologies(array_orders[i], randomize=True))+";")

	best_score = 9999999999

	for t in trees[:1]:
		# Initialize tree as dendropy Tree instance
		tree = dendropy.Tree.get(data=t, schema="newick")
		# Add internal nodes and infer their states.
		array_dict = {}
		for node in tree.seed_node.postorder_iter():
			if node.taxon:
				print(node.taxon.label)



if __name__ == '__main__':
	main()
