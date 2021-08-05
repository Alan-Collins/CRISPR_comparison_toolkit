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
from itertools import permutations, product
from string import ascii_lowercase
from copy import deepcopy
from collections import Counter



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


def swap_leaves(tree, num_swaps, seed, array_dict, all_arrays, event_costs):
	""" Swap 2 random leaves on the tree and recalculate ancestral states and total branch lengths. If improved, keep the tree, if not try another 2 random leaves. Try until consecutive failures to improve tree reaches num_swaps.
	Args:
		tree (dendropy Tree class): The tree to be modified.
		num_swaps (int): Number of leaf swaps that need to not yield improvement before the search is aborted.
		seed (int): seed value for random choices of leaves to swap.
		array_dict (dict): The Array class instances of the arrays in the tree.
		all_arrays (list): List of lists of the spacers in each array. [[spacer, in, array1], [spacers, in, array2]]
		event_costs (dict): Parsimony costs of different events.
	
	Returns:
		(dendropy Tree class) Best tree found during leaf-swapping search.
	"""

	swap_count = 0
	best_score = tree.length()
	best_tree = dendropy.Tree(tree)
	best_arrays = deepcopy(array_dict)

	while swap_count < num_swaps:
		random.seed(seed)
		leaf_a = random.choice(tree.leaf_nodes())
		sib_a = leaf_a.sibling_nodes()[0]
		seed += 1
		random.seed(seed)
		leaf_b = random.choice(tree.leaf_nodes())
		seed += 1
		while leaf_a == leaf_b or sib_a == leaf_b:
			random.seed(seed)
			leaf_b = random.choice(tree.leaf_nodes())
			seed += 1
		sib_b = leaf_b.sibling_nodes()[0]
		parent_a = leaf_a.parent_node
		parent_b = leaf_b.parent_node
		parent_a.set_child_nodes([leaf_b, sib_a])
		parent_b.set_child_nodes([leaf_a, sib_b])

		Incomplete_tree = False # Catch instances of sibling nodes sharing 0 spacers and abort the round of leaf swapping.
		for node in [leaf_a, leaf_b]: # For both leaves work back to root and recalculate all ancestors and branch lengths.
			while node.level() != 0:
				node_id = node.taxon.label
				sib = node.sibling_nodes()[0]
				sib_id = sib.taxon.label
				parent = node.parent_node
				parent_id = parent.taxon.label
				results = CRISPR_mp.resolve_pairwise_parsimony(array_dict[node_id], array_dict[sib_id], all_arrays, array_dict, [0], 0, tree, event_costs)
				if results == "No_ID":
					Incomplete_tree = True
					break
				
				array_dict[node_id].reset()
				array_dict[sib_id].reset()
				array_dict[node_id], array_dict[sib_id], ancestor = results
				array_dict[parent_id] = ancestor
				array_dict[parent_id].id = parent_id # Want to change everything except the ID

				# calculate branches including context of arrays elsewhere in the tree.
				array_dict[node_id] = CRISPR_mp.count_parsimony_events(array_dict[node_id], array_dict[parent_id], array_dict, tree, True)
				array_dict[sib_id] = CRISPR_mp.count_parsimony_events(array_dict[sib_id], array_dict[parent_id], array_dict, tree, True)

				# Refresh branch lengths of the children of the newly inferred ancestor
				for n, n_id in [(node, node_id), (sib, sib_id)]:
					for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
						array_dict[n_id].distance += array_dict[n_id].events[k] * v
					n.edge_length = array_dict[n_id].distance

				# Set node to its parent to work up towards the root node.
				node = node.parent_node
			if Incomplete_tree:
				break
		
		if tree.length() < best_score:
			swap_count = 0
			best_score = tree.length()
			best_tree = dendropy.Tree(tree)
			best_arrays = deepcopy(array_dict)
		else:
			swap_count += 1

	return best_tree, best_arrays


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
		"-w",  dest="leaf_swaps", type=int, nargs="?", default = 0,
			help="Specify number of times leaves should be swapped when optimising the leaf placement on a tree. Default: 0"
		)
	parser.add_argument(
		"-q",  dest="acquisition", type=int, nargs="?", default = 1,
			help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
		)
	parser.add_argument(
		"-i",  dest="indel", type=int, nargs="?", default = 10,
			help="Specify the parsimony cost of an indel event involving one or more spacers. Default: 1"
		)
	parser.add_argument(
		"-z",  dest="rep_indel", type=int, nargs="?", default = 50,
			help="Specify the parsimony cost of an indel event involving one or more spacers that is independently acquired in multiple arrays. Default: 50"
		)
	parser.add_argument(
		"-d",  dest="duplication", type=int, nargs="?", default = 1,
			help="Specify the parsimony cost of a duplication event involving one or more spacers. Default: 1"
		)
	parser.add_argument(
		"-l",  dest="trailer_loss", type=int, nargs="?", default = 1,
			help="Specify the parsimony cost of the loss of a spacer from the trailer end of the array. Default: 1"
		)
	parser.add_argument(
		"-o", dest="output_tree", required = False,
		help="Specify filename for the graphical representation of your tree with hypothetical intermediate arrays as a png."
		)
	parser.add_argument(
		"-s",  dest="seed", type=int, nargs="?", default = 0,
			help="Specify the seed to control the generation of the initial tree and leaf order. Use of this option results in the production of a single tree that will always be yielded by this seed value. This overrules the replicates option."
		)
	parser.add_argument(
		"arrays_to_join", nargs="*",  
		help="Specify the IDs of the arrays you want to join. If none provided, joins all arrays in the provided array representatives file. **If given, must come at the end of your command after all other arguments.**"
		)

	args = parser.parse_args(sys.argv[1:])

	event_costs = { 
					"acquisition" : args.acquisition,
					"indel" : args.indel,
					"repeated_indel" : args.rep_indel,
					"duplication": args.duplication,
					"trailer_loss": args.trailer_loss
					}

	array_spacers_dict = {}
	with open(args.array_file, 'r') as fin:
		for line in fin.readlines():
			bits = line.split()
			array_spacers_dict[bits[0]] = bits[2:]

	node_ids = ["Int " + i for i in ascii_lowercase]
	if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
		node_ids += ["Int " + "".join(i) for i in product(ascii_lowercase, repeat=(len(args.arrays_to_join)//26)+1)]

	# hex values from this website http://phrogz.net/css/distinct-colors.html

	Cols_hex_27 = ['#fd5925', '#dbc58e', '#008d40', '#304865', '#934270', '#f7b8a2', '#907500', '#45deb2', '#1f4195', '#d67381', '#8e7166', '#afb200', '#005746', '#a598ff', '#8f0f1b', '#b96000', '#667f42', '#00c7ce', '#9650f0', '#614017', '#59c300', '#1a8298', '#b5a6bd', '#ea9b00', '#bbcbb3', '#00b0ff', '#cd6ec6']

	#hex values from https://mokole.com/palette.html

	Cols_hex_40 = ["#696969","#556b2f","#a0522d","#800000","#006400","#808000","#483d8b","#3cb371","#008080","#bdb76b","#4682b4","#000080","#9acd32","#32cd32","#daa520","#7f007f","#ff4500","#00ced1","#ff8c00","#c71585","#0000cd","#00ff00","#9400d3","#dc143c","#00bfff","#f4a460","#adff2f","#da70d6","#ff00ff","#1e90ff","#db7093","#fa8072","#ffff54","#dda0dd","#7b68ee","#afeeee","#98fb98","#7fffd4","#ffe4c4","#ffc0cb"]

	Cols_tol = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"]

	Cols_hex_12 = ["#07001c", "#ff6f8d", "#4c62ff", "#92ffa9", "#810087", "#bcffe6", "#490046", "#00c8ee", "#b53900", "#ff8cf7", "#5b5800", "#14d625"]

	all_spacers = []
	for array in args.arrays_to_join:
		all_spacers += array_spacers_dict[array]
	non_singleton_spacers = [spacer for spacer, count in Counter(all_spacers).items() if count >1]
	if len(non_singleton_spacers) > 8:
		if len(non_singleton_spacers) > 12: 
			if len(non_singleton_spacers) > 27:
				if len(non_singleton_spacers) > 40:
					print("{} spacers found in multiple arrays. Using fill and outline colour combinations to distinguish spacers.".format(len(non_singleton_spacers)))
					if len(non_singleton_spacers) < 65:
						col_scheme = Cols_tol
					elif len(non_singleton_spacers) < 145:
						col_scheme = Cols_hex_12
					else:
						col_scheme = Cols_hex_27
					colours = []
					for i in range((len(non_singleton_spacers)+len(col_scheme)-1)//len(col_scheme)): # Repeat the same colour scheme.
						for j in col_scheme:
							colours += [(j, col_scheme[i])]

				else:
					colours = [(i, "#000000") for i in Cols_hex_40]
			else:
				colours = [(i, "#000000") for i in Cols_hex_27]
		else:
			colours = [(i, "#000000") for i in Cols_hex_12]
	else:
		colours = [(i, "#000000") for i in Cols_tol]
	# build a dictionary with colours assigned to each spacer.
	spacer_cols_dict  = {}

	for i, spacer in enumerate(sorted(non_singleton_spacers)):
		spacer_cols_dict[spacer] = colours[i]

	taxon_namespace = dendropy.TaxonNamespace(args.arrays_to_join + node_ids)
	
	if len(args.arrays_to_join) < 7 and not args.seed:
		array_orders = [i for i in permutations(args.arrays_to_join)]
		trees = []
		for order in array_orders:
			trees += [t + ";" for t in make_topologies(order)]
			seeds = [0 for _ in trees]
	else:
		if not args.seed:
			seeds = [random.randint(0,9999999999) for i in range(args.replicates)]
			trees = []
			for i in range(args.replicates):
				random.seed(seeds[i])
		
				trees.append(next(make_topologies(random.sample(args.arrays_to_join, len(args.arrays_to_join)), randomize=True))+";")
		else:
			random.seed(args.seed)
			seeds = [args.seed]
			trees = [next(make_topologies(random.sample(args.arrays_to_join, len(args.arrays_to_join)), randomize=True))+";"]



	best_score = 9999999999
	best_seed = []
	best_tree = []
	best_arrays = []
	for seed, t in zip(seeds, trees):
		# Control settings
		Incomplete_tree = False
		node_count = 0

		# Initialize tree as dendropy Tree instance
		tree = dendropy.Tree.get(data=t, schema="newick", taxon_namespace=taxon_namespace)

		# Initialize array dict
		array_dict = {}
		for array in args.arrays_to_join:
			array_dict[array] = CRISPR_mp.Array(array, array_spacers_dict[array], extant=True)
		all_arrays = [array.spacers for array in array_dict.values()]

		# Add internal nodes and infer their states.
		for node in tree.seed_node.postorder_iter():
			if node.taxon:
				a = node.taxon.label 
				if len(node.sibling_nodes()) != 0:
					sister = node.sibling_nodes()[0]
					parent = node.parent_node
					if sister.taxon:
						b = sister.taxon.label
						if not parent.taxon: # Only need to add the ancestor once.
							results = CRISPR_mp.resolve_pairwise_parsimony(array_dict[a], array_dict[b], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
							if results == "No_ID":
								Incomplete_tree = True
								break
							else:
								array_dict[a], array_dict[b], ancestor = results
								node_count+=1

							for n, n_id in [(node, a), (sister, b)]:
								for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
									array_dict[n_id].distance += array_dict[n_id].events[k] * v
								n.edge_length = array_dict[n_id].distance

							array_dict[ancestor.id] = ancestor
							node.parent_node.taxon = taxon_namespace.get_taxon(ancestor.id)
							node.parent_node.edge_length = 0 # Start the ancestor with 0 branch length
		
		if Incomplete_tree: 
			# If there was no identity between neighbouring arrays in the tree then skip this replicate and try again.
			continue

		# Repeat iteration now that tree is built to add repeat indels.
		for node in tree:
			if node.level() != 0:
				parent = node.parent_node				
				array_dict[node.taxon.label].reset()
				array_dict[parent.taxon.label].reset()
				array_dict[node.taxon.label] = CRISPR_mp.count_parsimony_events(array_dict[node.taxon.label], array_dict[parent.taxon.label], array_dict, tree, True)
				for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
					array_dict[node.taxon.label].distance += array_dict[node.taxon.label].events[k] * v
				node.edge_length = array_dict[node.taxon.label].distance
		
		# If user requested leaf swapping to optimize leaf positions, do that.

		if args.leaf_swaps:
			tree, array_dict = swap_leaves(tree, args.leaf_swaps, seed, array_dict, all_arrays, event_costs)
		
		
		tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works					
		score = tree.length()
		if score < best_score:
			best_score = score
			if seed:
				best_seed = [str(seed)]
			best_tree = [dendropy.Tree(tree)]
			best_arrays = [deepcopy(array_dict)]
		elif score == best_score:
			if not any([
						dendropy.calculate.treecompare.weighted_robinson_foulds_distance(good_tree, tree) == 0. for good_tree in best_tree
						]):
				if seed:
					best_seed.append(str(seed))
				best_tree.append(dendropy.Tree(tree))
				best_arrays.append(deepcopy(array_dict))

	if best_tree == []:
		print("Unable to generate trees due to lack of identity between arrays. Try increasing the number of replicates used. Otherwise you may want to consider splitting these arrays into smaller groups that share more spacers.")
		sys.exit()


	print("\n\nThe best score for tree(s) was: {}".format(best_score))
	if best_seed != []:
		print("\n\nThe seed values to recreate the best tree(s) are: {}\n\n".format(", ".join(best_seed)))
	for n, tree in enumerate(best_tree):	
		print(tree.as_ascii_plot(show_internal_node_labels=True))
		print(tree.as_string("newick"))
		if n > 0:
			filename = "{}_{}.png".format(args.output_tree[:-4], n+1)
		else:
			filename = args.output_tree
		CRISPR_mp.plot_tree(tree, best_arrays[n], filename, spacer_cols_dict, branch_lengths=True, emphasize_diffs=True)


if __name__ == '__main__':
	main()
