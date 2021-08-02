#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-8-2
# DESCRIPTION :  Generate random bifircating tree topologies and calculate parsimony of each then return the best.

import sys
import argparse
import dendropy
import CRISPR_mp
from random import randint


def make_topologies(leaves, root=False):
	"""
	Given a number of leaves for a tree, generate all possible configurations of the rooted bifurcating tree.
	Args:
		num_leaves (list): Leaves to be placed in the tree.
		root (int): The 0 base integer location at which to place the root and generate all tree topologies with that root.
	Yields:
		 Newick string of tree or subtree.
	"""
	if len(leaves) == 1: # When all nodes have been joined, return the list in a tidied format.
		yield leaves[0]
	else:
		if root:
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


	print([t for t in make_topologies(args.arrays_to_join, 2)])



if __name__ == '__main__':
	main()
