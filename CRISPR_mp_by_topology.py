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
	if len(leaves) == 1: # When all nodes have been joined, return the list in a tidied format.
		yield leaves[0]
	else:
		if randomize:
			root = random.randint(1,len(leaves)-1) # Root the tree at a random point
			for left in make_topologies(leaves[:root], random.randint(1, root), True): # Root the subtree to the left at a random point
				for right in make_topologies(leaves[root:], random.randint(root, len(leaves)-1), True): # Root the subtree to the right at a random point
					yield "({},{})".format(left, right)
		elif root:
			for left in make_topologies(leaves[:root]): # Make subtrees to the left of (sub)tree root
				for right in make_topologies(leaves[root:]): # Make subtrees to the right of (sub)tree root
					yield "({},{})".format(left, right)
		else:
			for i in range(1,len(leaves)):
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

	array_choices = [random.sample(args.arrays_to_join, len(args.arrays_to_join)) for i in range(args.replicates)]
	
	starting_trees = []
	for order in array_choices:
		t = next(make_topologies(order, randomize=True)) + ";"
		starting_trees.append(dendropy.Tree.get(data=t, schema="newick"))

	for tree in starting_trees:
		print(tree.as_ascii_plot())



if __name__ == '__main__':
	main()
