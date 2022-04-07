#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# DESCRIPTION :  Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship.

import sys
import argparse
import numpy as np
from itertools import product
from string import ascii_uppercase
import random
from collections import Counter, defaultdict
import dendropy
import copy
from itertools import permutations
import matplotlib.pyplot as plt
from matplotlib import rcParams
import multiprocessing
from math import ceil, log
import time
from datetime import timedelta
import json

from . import (
	file_handling,
	colour_schemes,
	array_parsimony,
	sequence_operations,
	tree_operations,
	plotting)

description = """
usage: cctk CRISPRtree [-h] -a [-o] [--output-arrays] [--print-tree] [-x] \
[-r] [--acquisition] [--deletion] [--insertion] [--rep-indel] \
[--duplication] [--trailer-loss] [--no-ident] [-t] [--seed] \
[--colour-file] [--colour-scheme-outfile] [--colour-scheme-infile] [-e] \
[-b] [--brlen-scale] [--no-align-cartoons] [--no-align-labels] [--dpi] \
[--no-fade-anc] [--plot-width] [--plot-height] [--font-override-labels] \
[--font-override-annotations] [arrays_to_join]

optional arguments:
  -h, --help        show this help message and exit

positional arguments:
  arrays_to_join    IDs of the arrays you want to analyse. Default: all

required arguments:
  -a                Array_IDs.txt or Array_seqs.txt

output control:
  set which of the optional outputs you want

  -o, --out-file    output plot file name
  --output-arrays   file to store analyzed arrays and hypothetical ancestors
  --print-tree      print an ascii symbol representation of the tree

running parameters:
  control run behaviour

  -x, --fix-order   only build one tree using the provided order of arrays
  -r, --replicates  number of replicates of tree building. Default: 100
  --acquisition     parsimony cost of a spacer acquisition event. Default: 1
  --deletion        parsimony cost of a deletion event. Default: 10
  --insertion       parsimony cost of an insertion event. Default: 30
  --rep-indel       parsimony cost independently acquiring spacers. Default: 50
  --duplication     parsimony cost of a duplication event. Default: 1
  --trailer-loss    parsimony cost of trailer spacer loss. Default: 1
  --no-ident        parsimony cost of an array having no identity with its \
ancestor. Default: 100
  
  -t, --num-threads 
                    number of threads to use. Default: 1
  --seed            set seed for random processes

colour scheme files:
  set inputs and outputs for optional colour scheme files

  --colour-file     file with custom colour list
  --colour-scheme-outfile
                    output file to store json format colour schemes
  --colour-scheme-infile
                    input file json format colour scheme

plotting parameters:
  control elements of the produced plot

  -b                include branch lengths in tree plot
  --brlen-scale     factor to scale branch length
  --branch-support  Show support at nodes
  --no-emphasize-diffs
                    don't emphasize events in each array since its ancestor
  --no-align-cartoons
                    draw array cartoons next to leaf node
  --no-align-labels
                    draw leaf labels next to leaf node
  --dpi             resolution of the output image. Default = 600
  --no-fade-anc     do not apply transparency to ancestral array depiction
  --plot-width      width of plot in inches. Default = 3
  --plot-height     height of plot in inches. Default = 3
  --font-override-labels
                    set label font size in pts
  --font-override-annotations
                    set annotation font size in pts
"""

class NodeBS(dendropy.Node):
	def __init__(self, **kwargs):
		dendropy.Node.__init__(self, **kwargs)
		self.node_support=0

	def write_newick_bs(self, out, **kwargs):
		# adapted from 
		# https://github.com/jeetsukumaran/DendroPy/blob/29fd294bf05d890ebf6a8d576c501e471db27ca1/src/dendropy/datamodel/treemodel.py#L2439
		edge_lengths = not kwargs.get('suppress_edge_lengths', False)
		edge_lengths = kwargs.get('edge_lengths', edge_lengths)
		child_nodes = self.child_nodes()
		if child_nodes:
			out.write('(')
			f_child = child_nodes[0]
			for child in child_nodes:
				if child is not f_child:
					out.write(',')
				child.write_newick_bs(out, **kwargs)
			out.write(')')

		out.write(self._get_node_token(**kwargs))
		if edge_lengths:
			e = self.edge
			if e:
				sel = e.length
				if sel is not None:
					fmt = kwargs.get('edge_length_formatter', None)
					if fmt:
						out.write(":%s" % fmt(sel))
					else:
						s = ""
						try:
							s = float(sel)
							s = str(s)
						except ValueError:
							s = str(sel)
						if s:
							out.write(":%s" % s)
		if child_nodes and self.parent_node:
			out.write('[%i]' % self.node_support)
		return out.getvalue()


def find_closest_array(array, array_dict, tree, event_costs):
	"""
	Args:
		array (Array class instance): The array you want to find the closest match for.
		array_dict (dict): The dictionary with values of Array class instances of arrays already in your tree
		tree (Deondropy Tree class instance): The tree in which the arrays are located.
		event_costs (dict): Dict to look up event types and their parsimony costs.

	Returns:
		(Array class instance) The array already in your tree that is the most parsimonious match to the query array.
	"""

	best_score = 9999999999
	best_match = False
	for array_id in array_dict.keys():
		comparator_array = copy.deepcopy(array_dict[array_id])
		comparator_array.reset()
		array.reset()
		comparator_array = array_parsimony.count_parsimony_events(comparator_array, array, array_dict, tree, False)
		# calc comarator_array_dist now. It's reset during next comparison
		comparator_array_dist = 0
		for k,v in event_costs.items():
			comparator_array_dist += comparator_array.events[k] * v
		array = array_parsimony.count_parsimony_events(array, comparator_array, array_dict, tree, True)
		for k,v in event_costs.items():
			array.distance += array.events[k] * v
		if any([i.type == "shared" for i in comparator_array.modules]):
			if array.distance+comparator_array_dist < best_score:
				best_score = array.distance+comparator_array_dist
				best_match = comparator_array
			elif array.distance+comparator_array_dist == best_score:
				if best_match.extant and not comparator_array.extant:
					# Prefer to join arrays to existing ancestors rather than to leaves if they are the same distance.
					best_match = comparator_array
				if not best_match.extant:
					# Check if this is an internal polytomy.
					# If it is then set best match to ultimate parent
					node = tree.find_node_with_taxon_label(best_match.id)
					if node.edge_length == 0:
						while node.edge_length == 0 and node.parent_node is not tree.seed_node:
							node = node.parent_node
						best_match = copy.deepcopy(array_dict[node.taxon.label])
	if not best_match:
		return "No_ID"


	if not any([i.type == "shared" for i in best_match.modules]):
		return "No_ID"


	else:
		return best_match


def replace_existing_array(existing_array, new_array, current_parent, tree,
	all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed,
	event_costs, tree_namespace):
	"""
	Given a array already in the tree and an array to be added, finds an ancestral state for the two arrays and replaces the existing array with the following tree structure:
			   /- existing array
	- ancestor|
	           \- new array
	Args:
		existing_array (Array class instance): The array that has already been added to the tree.
		new_array (Array class instance): The array you want to add to the tree.
		current_parent (Array class instance): ID of the parent node of the existing array in the tree
		tree (dendropy Tree instance): The tree object to be modified
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): A list of names to be used to name internal nodes.
		node_count (int): The index of the name to use in the node_ids list to name this internal node.
		array_dict (dict): Dict of arrays in the tree so far.
		tree_child_dict (dict): Dict of dendopy tree node class instances in the tree.
		seed (bool): boolean indicating if the parent of the existing array the seed node in the tree.
		event_costs (dict): Dict to look up event types and their parsimony costs.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
	Returns:
		(tuple) of the following:
			(dendropy Tree instance) Tree with the existing array replaced with the newly formed node.
			(dict) Dict of arrays in the tree so far.
			(dict) Dict of dendopy tree node class instances in the tree.
	"""

	# Make a hypothetical ancestor for the pair of arrays and calculate the distances
	results = array_parsimony.resolve_pairwise_parsimony(existing_array, new_array, all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	
	if results == "No_ID":
		return results

	else:
		existing_array, new_array, ancestor = results
		tree_child_dict[existing_array.id].edge_length = existing_array.distance
		# If the current parent is not seed then do this
		if not seed:
			# Calculate distance from this hypothetical ancestor to its new parent in the tree
			ancestor.reset()
			ancestor = array_parsimony.count_parsimony_events(ancestor, current_parent, array_dict, tree, True)
			for k,v in event_costs.items():
				ancestor.distance += ancestor.events[k] * v

		# modify the edge length for the existing array and update its distance in the array_dict
		tree_child_dict[existing_array.id].edge_length = existing_array.distance
		array_dict[existing_array.id] = existing_array
		for a in [new_array, ancestor]:
			# Create tree nodes
			tree_child_dict[a.id] = NodeBS(edge_length=a.distance)
			tree_child_dict[a.id].taxon = tree_namespace.get_taxon(a.id)
			#Store arrays for further comparisons
			array_dict[a.id] = a
		# Build new tree node
		tree_child_dict[ancestor.id].set_child_nodes([tree_child_dict[existing_array.id], tree_child_dict[new_array.id]])

		if not seed:
			# figure out which were the existing children of the current_parent
			for child in tree_child_dict[current_parent.id].child_node_iter():
				if child.taxon.label != existing_array.id:
					other_child = child.taxon.label
			# Add new node to existing tree as a child of the previous parent of the existing array.
			tree_child_dict[current_parent.id].set_child_nodes([tree_child_dict[ancestor.id], tree_child_dict[other_child]])
			


		else:
			tree.seed_node.set_child_nodes([tree_child_dict[ancestor.id]])
			
		return tree, array_dict, tree_child_dict


def build_tree_single(arrays, tree_namespace, score, all_arrays, node_ids,
	event_costs, branch_support=False):
	"""
	Search treespace for most parsimonious tree using single process.
	Args:
		arrays (list): Ordered list of Array class instances of the arrays to analyse. Will be added to the tree in the provided order.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
		score (int): The score to beat. If at any point during tree construction the tree has total branch lengths above this score construction will be aborted
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): The names for internal nodes in the tree to be assigned.
		event_costs (dict): Dict to look up event types and their parsimony costs.
		branch_support (bool): Should tree always be returned to be included in branch support calculation.
	
	Returns:
		(tuple) Returns array_dict and dendropy.tree object if the tree beats or equals the provide score (or if branch_support = True). Else retuns tuple of (False, False).
	"""
	initial_arrays = copy.deepcopy(arrays)
	tree = dendropy.Tree(taxon_namespace=tree_namespace)
	array_dict = {}
	tree_child_dict = {}
	node_count = 0 # Keep track of which internal node ID should be used for each node
	# Remove the arrays being compared from checks to see what spacers are in other arrays.
	# That way we only worry about ancestral states accounting for deletion of spacers in arrays not yet added to the tree.
	all_arrays = [a for a in all_arrays if a not in [arrays[0].spacers, arrays[1].spacers]]
	results = array_parsimony.resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)

	while results == "No_ID":
		arrays.append(arrays[1])
		del arrays[1]
		results = array_parsimony.resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays,  array_dict, node_ids, node_count, tree, event_costs)
	node_count += 1
	array1, array2, ancestor = results
	for a in [array1, array2, ancestor]:
		# Create tree nodes
		tree_child_dict[a.id] = NodeBS(edge_length=a.distance)
		tree_child_dict[a.id].taxon = tree_namespace.get_taxon(a.id)
		#Store arrays for further comparisons
		array_dict[a.id] = a


	# Add initial relationships to the tree.
	tree.seed_node.add_child(tree_child_dict[ancestor.id])
	tree_child_dict[ancestor.id].add_child(tree_child_dict[array1.id])
	tree_child_dict[ancestor.id].add_child(tree_child_dict[array2.id])

	if len(arrays) != 2:
		for i in range(2, len(arrays)): # Already added the first two so now add the rest 1 by 1
			a = arrays[i]
			all_arrays = [arr for arr in all_arrays if arr != a.spacers]
			seed = False # To check if we are modifying the child of the seed node
			# Find the most similar array already in the tree (measured in parsimony score)
			best_match = find_closest_array(a, array_dict, tree, event_costs)
			while best_match == "No_ID":
				arrays.append(arrays[i])
				del arrays[i]
				a = arrays[i]
				best_match = find_closest_array(a, array_dict, tree, event_costs)
			if best_match.extant: # If the closest match is a array then just join them and replace the existing array with the new node
				current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
				(
					tree, array_dict, tree_child_dict
				) = replace_existing_array(
					best_match,
					a,
					current_parent,
					tree,
					all_arrays,
					node_ids,
					node_count,
					array_dict,
					tree_child_dict,
					seed,
					event_costs,
					tree_namespace
					)
				node_count += 1
			else:
				try:
					current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
				except: # Fails if the best match is a child of the seed node
					seed = True
					current_parent = None
				(
					tree, array_dict, tree_child_dict
				) = replace_existing_array(
					best_match,
					a,
					current_parent,
					tree,
					all_arrays,
					node_ids,
					node_count,
					array_dict,
					tree_child_dict,
					seed,
					event_costs,
					tree_namespace
					)
				node_count += 1

			# Recheck child - ancestor branch length to find indels that would have to occur multiple times
			for node in tree.postorder_node_iter():
				if node.level() != 0:
					if node.parent_node.level() != 0:
						node_array = array_dict[node.taxon.label]
						parent_array = array_dict[node.parent_node.taxon.label]
						node_array.reset()
						parent_array.reset()
						node_array = array_parsimony.count_parsimony_events(node_array, parent_array, array_dict, tree, True)
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							node_array.distance += node_array.events[k] * v
						node.edge_length = node_array.distance

		
			brlen = tree.length()
			if brlen > score and not branch_support:
				return (False, False, False)

					
			if a == arrays[-1]:
				tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works
				return (array_dict, tree, initial_arrays)


	else:
		return (array_dict, tree, initial_arrays)


def build_tree_multi(arrays, tree_namespace, all_arrays, node_ids,
	event_costs):
	"""
	Search treespace for most parsimonious tree using multiple processes.
	Args:
		arrays (list): Ordered list of Array class instances of the arrays to analyse. Will be added to the tree in the provided order.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): The names for internal nodes in the tree to be assigned.
		event_costs (dict): Dict to look up event types and their parsimony costs.
	
	Returns:
		(tuple) Returns array_dict and dendropy.tree object.
	"""

	# Remove the arrays being compared from checks to see what spacers are in other arrays.
	# That way we only worry about ancestral states accounting for deletion of spacers in arrays not yet added to the tree.
	all_arrays = [a for a in all_arrays if a not in [arrays[0].spacers, arrays[1].spacers]]
	initial_arrays = copy.deepcopy(arrays)
	tree = dendropy.Tree(taxon_namespace=tree_namespace)

	array_dict = {}
	tree_child_dict = {}
	node_count = 0 # Keep track of which internal node ID should be used for each node

	results = array_parsimony.resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	while results == "No_ID":

		arrays.append(arrays[1])
		del arrays[1]
		results = array_parsimony.resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	node_count += 1
	array1, array2, ancestor = results

	for a in [array1, array2, ancestor]:
		# Create tree nodes
		tree_child_dict[a.id] = NodeBS(edge_length=a.distance)
		tree_child_dict[a.id].taxon = tree_namespace.get_taxon(a.id)
		#Store arrays for further comparisons
		array_dict[a.id] = a

	# Add initial relationships to the tree.
	tree.seed_node.add_child(tree_child_dict[ancestor.id])
	tree_child_dict[ancestor.id].add_child(tree_child_dict[array1.id])
	tree_child_dict[ancestor.id].add_child(tree_child_dict[array2.id])
	for i in range(2, len(arrays)): # Already added the first two so now add the rest 1 by 1
		a = arrays[i]
		all_arrays = [arr for arr in all_arrays if arr != a.spacers]
		seed = False # To check if we are modifying the child of the seed node

		# Find the most similar array already in the tree (measured in parsimony score)
		best_match = find_closest_array(a, array_dict, tree, event_costs)
		while best_match == "No_ID":
			arrays.append(arrays[i])
			del arrays[i]
			a = arrays[i]
			best_match = find_closest_array(a, array_dict, tree, event_costs)
		if best_match.extant: # If the closest match is a array then just join them and replace the existing array with the new node
			current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
			results = replace_existing_array(
				best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed, event_costs, tree_namespace
				)

			tree, array_dict, tree_child_dict = results
			node_count += 1
		else:
			try:
				current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
			except: # Fails if the best match is a child of the seed node
				seed = True
				current_parent = None
			results = replace_existing_array(
				best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed, event_costs, tree_namespace
				)

			tree, array_dict, tree_child_dict = results
			node_count += 1
		
		if a == arrays[-1]:
			# Recheck child - ancestor branch length to find indels that would have to occur multiple times
			for node in tree.postorder_node_iter():
				if node.level() != 0:
					if node.parent_node.level() != 0:
						node_array = array_dict[node.taxon.label]
						parent_array = array_dict[node.parent_node.taxon.label]
						node_array.reset()
						parent_array.reset()
						node_array = array_parsimony.count_parsimony_events(node_array, parent_array, array_dict, tree, True)
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							node_array.distance += node_array.events[k] * v
						node.edge_length = node_array.distance
			# Need to reroot at the seed so that RF distance works
			tree.reroot_at_node(tree.seed_node, update_bipartitions=False)
			return array_dict, tree, initial_arrays


def reset_anc_mods(tree, array_dict):
	"""Ensures modules identified in ancestor during comparison with
	child aren't plotted"""

	root_array_name = tree.seed_node.taxon.label
	array = array_dict[root_array_name]

	array.reset()
	array.aligned = array.spacers

	new_module = array_parsimony.SpacerModule()
	new_module.type = "shared"
	new_module.spacers = array.spacers
	new_module.indices = [i for i in range(len(new_module.spacers))]

	array.modules.append(new_module)
	array.module_lookup = {i: new_module for i in new_module.indices}

	array_dict[root_array_name] = array

	return array_dict


def build_parser(parser):
	parser.add_argument(
		"-a",
		dest="array_file",
		required=True,
		help="Specify array representatives file."
		)

	output_params = parser.add_argument_group('Output control', 
		"Set which of the optional outputs you want.")
	output_params.add_argument(
		"-o", "--out-file",
		metavar="",
		required=False,
		help="Specify filename for the graphical representation of your tree with hypothetical intermediate arrays."
		)
	output_params.add_argument(
		"--output-arrays",
		metavar="",
		required=False,
		help="Specify filename for the details of you final arrays with hypothetical intermediate arrays in the same format as your input array_file. If there are multiple best trees, one file will be created per tree numbered in the order they are described in the stdout output."
		)
	output_params.add_argument(
		"--print-tree",
		action='store_true',  
		help="Print a graphical representation of the tree using ascii characters."
		)

	run_params = parser.add_argument_group('Running parameters', 
		"Control run behaviour.")
	run_params.add_argument(
		"-x", "--fix-order",
		action='store_true',  
		help="Only build one tree using the provided order of arrays. Good for recreating previously found trees."
		)
	run_params.add_argument(
		"-r", "--replicates",
		metavar="",
		type=int,
		default=100,
		help="Specify number of replicates of tree building to perform. The more replicates, the greater the chance that a better tree will be found. Default: 100"
		)
	run_params.add_argument(
		"--acquisition",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
		)
	run_params.add_argument(
		"--deletion",
		metavar="",
		type=int,
		default=10,
		help="Specify the parsimony cost of a deletion event involving one or more spacers. Default: 10"
		)
	run_params.add_argument(
		"--insertion",
		metavar="",
		type=int,
		default=30,
		help="Specify the parsimony cost of an insertion event involving one or more spacers. Default: 30"
		)
	run_params.add_argument(
		"--rep-indel",
		metavar="",
		type=int,
		default=50,
		help="Specify the parsimony cost of an indel event involving one or more spacers that is independently acquired in multiple arrays. Default: 50"
		)
	run_params.add_argument(
		"--duplication",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of a duplication event involving one or more spacers. Default: 1"
		)
	run_params.add_argument(
		"--trailer-loss",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of the loss of a spacer from the trailer end of the array. Default: 1"
		)
	run_params.add_argument(
		"--no-ident",
		metavar="",
		type=int,
		default=100,
		help="Specify the parsimony cost of a tree in which an array is predicted to have descended from another array with which it shares no spacers. Default: 100"
		)
	
	run_params.add_argument(
		"-t", "--num-threads",
		metavar="",
		type=int,
		default=1,
		help="Specify number of threads to use for building trees. Using multiple threads will speed up the search for trees when performing many replicates with the -r option. Default: 1"
		)
	run_params.add_argument(
		"--seed",
		metavar="",
		type=int,
		required=False,
		default=2,
		help="The order of outline and fill colours assigned to spacers is semi-random. Change it by providing a number here to change which colours are assigned to each spacer."
		)

	cs_files = parser.add_argument_group('Colour scheme files', 
		"Set inputs and outputs for optional colour scheme files.")
	cs_files.add_argument(
		"--colour-file",
		metavar="",
		required=False, 
		help="Specify file with custom colour list (Optional). \
			Colours must be hex codes. One colour per line with no header \
			line in file. e.g. #fd5925."
		)
	cs_files.add_argument(
		"--colour-scheme-outfile",
		metavar="",
		required=False, 
		help="Specify output file to store json format dictionary of the colour schemes used for spacers in this run."
		)
	cs_files.add_argument(
		"--colour-scheme-infile",
		metavar="",
		required=False, 
		help="Specify input file containing json format dictionary of the colour scheme to be used for spacers in this run. Any spacers not in the input file will be coloured according to the normal process."
		)
	
	plot_params = parser.add_argument_group('Plotting parameters', 
		"Control elements of the produced plot.")
	plot_params.add_argument(
		"--no-emphasize-diffs",
		dest="emphasize_diffs",
		action='store_false',  
		help="When plotting a representation of the tree with cartooned arrays, don't emphasize locations where arrays differ from their hypothetical ancestor."
		)
	plot_params.add_argument(
		"-b",
		dest="branch_lengths",
		action='store_true',
		help="When plotting a representation of the tree Include branch lengths at midpoint on branches. N.B. This does not control inclusion of branch lengths in newick format tree output, which are always included."
		)
	plot_params.add_argument(
		"--brlen-scale",
		type=float,
		required=False,
		metavar='',
		default=1,
		help="Factor to scale branch length."
		)
	plot_params.add_argument(
		"--branch-support",
		action='store_true',
		help="Show support at nodes"
		)
	plot_params.add_argument(
		"--no-align-cartoons",
		action='store_true',
		default=False,
		help="When plotting a representation of the tree with cartooned arrays, this option controls whether those cartoons are drawn next to the corresponding node in the tree. By default cartoons will be aligned at the trailer end."
		)
	plot_params.add_argument(
		"--no-align-labels",
		action='store_true',
		default=False,
		help="Should node labels be placed at the corresponding internal node or leaf. By default they are all aligned. N.B. For this setting to work, cartoons must also placed at nodes using the -g option."
		)
	plot_params.add_argument(
		"--dpi",
		type=int,
		required=False,
		metavar='',
		default=600,
		help="The desired resolution of the output image."
		)
	plot_params.add_argument(
		"--no-fade-anc",
		action='store_true', 
		help="Do not de-emphasize ancestral array cartoons using transparency. This option helps to make it clear which are the inferred ancestral arrays and which are the arrays being analyzed. However, this option reduces colour contrast between spacers."
		)
	plot_params.add_argument(
		"--plot-width",
		type=float,
		default=3,
		metavar="",
		help="Width of plot in inches. Default = 3"
		)
	plot_params.add_argument(
		"--plot-height",
		type=float,
		default=3,
		metavar="",
		help="Height of plot in inches. Default = 3"
		)
	plot_params.add_argument(
		"--font-override-labels",
		type=float,
		metavar="",
		help="If you want to force a label font size to be used rather than \
			using scaling to determine font size, provide it here"
		)
	plot_params.add_argument(
		"--font-override-annotations",
		type=float,
		metavar="",
		help="If you want to force a specific font size (pts) to be used for \
			annotations such as branch length labels, rather than using \
			scaling to determine font size, provide it here"
		)

	parser.add_argument(
		"arrays_to_join",
		nargs="*",  
		help="Specify the IDs of the arrays you want to join. \
			If none provided, joins all arrays in the provided array \
			representatives file. **If given, must come at the end of your \
			command after all other arguments.**"
		)

	return parser


def main(args):

	start_time = time.time()

	if args.no_align_cartoons and not args.no_align_labels:
		args.no_align_labels = True

	sys.stderr.write("Running with the following command:\n{}\n".format(
		" ".join(sys.argv)))

	event_costs = { 
					"acquisition": args.acquisition,
					"deletion": args.deletion,
					"insertion": args.insertion,
					"repeated_indel": args.rep_indel,
					"duplication": args.duplication,
					"trailer_loss": args.trailer_loss,
					"no_ident": args.no_ident
					}

	array_spacers_dict = file_handling.read_array_file(args.array_file)

	if args.arrays_to_join:
		arrays = [array_parsimony.Array(
			i, array_spacers_dict[i]) for i in args.arrays_to_join]
		labels = args.arrays_to_join
	else:
		arrays = [array_parsimony.Array(
			i, array_spacers_dict[i]) for i in array_spacers_dict.keys()]
		labels = [l for l in array_spacers_dict.keys()]

	all_arrays = [array.spacers for array in arrays]

	# Generate strings to assign as internal node_IDs (This makes 702)

	node_ids = tree_operations.create_internal_node_ids(
		len(labels), "Anc ")

	# Check if any arrays share no spacers with any others. If so, exit with error message.
	array_parsimony.check_for_no_ident(arrays)
		
	random.seed(args.seed)
	if len(labels) < 9:
		array_choices = [copy.deepcopy(list(i)) for i in permutations(arrays, len(arrays))]
		random.shuffle(array_choices)

		if len(array_choices) > args.replicates:
			sys.stderr.write("\nThere are {} possible trees to check. If you "
				"want to check every possible tree then set -r {}\n".format(
					len(array_choices), len(array_choices)))
			array_choices = [array_choices[i] for i in range(args.replicates)]

		elif len(array_choices) < args.replicates:
			sys.stderr.write("\nThere are only {} possible trees to check. "
				"You specified a greater number of replicates than there are "
				"possible trees. All possible trees will be checked.\n".format(
					len(array_choices)))

		else:
			sys.stderr.write("\nYou specified a number of replicates equal to "
				"the number of possible trees. All possible trees will be "
				"checked.\n")
	else:
		array_choices = [
		random.sample(copy.deepcopy(arrays), len(arrays)) for i in range(
			args.replicates)]

	taxon_namespace = dendropy.TaxonNamespace(labels + node_ids)

	best_score = 99999999

	# Keep track of how many times no common spacers are found.
	# If no common spacers are ever found then print an error message saying so

	no_id_count = 0 

	num_threads = args.num_threads if not args.fix_order else 1

	if len(array_choices) < num_threads:
		num_threads = 1

	if args.branch_support:
		tree_list = []
	if num_threads == 1:
		for i in range(min([args.replicates, len(array_choices)])):
			if args.fix_order:
				if not args.arrays_to_join:
					sys.exit("You must provide the order you want to fix when "
						"using the fixed order option!\n\nABORTING.")
				addition_order = [
				array_parsimony.Array(
					i, array_spacers_dict[i]) for i in args.arrays_to_join]
			else:
				addition_order = array_choices[i]
			try:
				array_dict, tree, addition_order = build_tree_single(
					addition_order,
					taxon_namespace,
					best_score,
					all_arrays,
					node_ids,
					event_costs,
					args.branch_support)
			except Exception as e:
				import os
				exc_type, exc_obj, exc_tb = sys.exc_info()
				fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
				print(sys.exc_info())
				print(exc_type, fname, exc_tb.tb_lineno)
				no_id_count += 1
				continue

			# Returns false if tree exeeds best score while building
			if not tree:
				continue
			if args.branch_support:
				# tuple to be consistent with multithread structure
				tree_list.append((array_dict,tree,addition_order)) 
			score = tree.length()
			if score < best_score:
				# Keep inferred ancestral states and information
				best_arrays = copy.deepcopy(array_dict)
				best_score = copy.deepcopy(score)
				best_addition_order = copy.deepcopy(addition_order)
				# Keep one copy for comparisons as copy.deepcopy makes a new taxon namespace which breaks comparisons.
				best_tree_comparator = dendropy.Tree(tree)
				best_tree = copy.deepcopy(tree)
			elif score == best_score:
				if isinstance(best_tree, list):
					# Check this tree isn't identical to one that's already been found
					if not any(tree_operations.compare_to_trees(
							tree, best_tree_comparator)):
						best_tree_comparator.append(dendropy.Tree(tree))
						best_tree.append(copy.deepcopy(tree))
						best_arrays.append(copy.deepcopy(array_dict))
						best_addition_order.append(copy.deepcopy(addition_order))
				else:
					# Check this tree isn't identical to the one that's already been found
					if tree_operations.compare_to_trees(
							tree, best_tree_comparator):
						best_tree_comparator = [
							best_tree_comparator, dendropy.Tree(tree)]
						best_tree = [best_tree, copy.deepcopy(tree)]
						best_arrays = [best_arrays, copy.deepcopy(array_dict)]
						best_addition_order = [
							best_addition_order, copy.deepcopy(addition_order)]
			if args.fix_order:
				break

				


	else:
		pool = multiprocessing.Pool(processes=num_threads)
		chunksize = args.replicates//num_threads
		options = [(array_list, taxon_namespace, all_arrays, node_ids, event_costs) for array_list in array_choices]
		tree_list = pool.starmap(build_tree_multi, options, chunksize)
		pool.close()
		pool.join()

		all_tree_comparators = dendropy.TreeList()
		for array_dict, tree, addition_order in tree_list:
			if array_dict:
				all_tree_comparators.read(
					data=tree.as_string("newick"), schema="newick")
				score = tree.length()
				if score < best_score:
					# Keep inferred ancestral states and information
					best_arrays = copy.deepcopy(array_dict) 
					best_score = copy.deepcopy(score)
					best_addition_order = copy.deepcopy(addition_order)
					# Keep one copy for comparisons as copy.deepcopy makes 
					# a new taxon namespace which breaks comparisons.
					best_tree_comparator = all_tree_comparators[-1:]
					best_tree = copy.deepcopy(tree)
				elif score == best_score:
					if isinstance(best_tree, list):
						# Check this tree isn't identical to one that's 
						# already been found
						if not any(
							tree_operations.compare_to_trees(
								all_tree_comparators[-1],
								best_tree_comparator)):
							best_tree_comparator.append(
								all_tree_comparators[-1])
							best_tree.append(copy.deepcopy(tree))
							best_arrays.append(copy.deepcopy(array_dict))
							best_addition_order.append(
								copy.deepcopy(addition_order))
					else:
						# Check this tree isn't identical to the one that's already been found
						if tree_operations.compare_to_trees(
								all_tree_comparators[-1],
								best_tree_comparator[0]) != 0.:
							best_tree_comparator.append(
								all_tree_comparators[-1])
							best_tree = [best_tree, copy.deepcopy(tree)]
							best_arrays = [
								best_arrays, copy.deepcopy(array_dict)]
							best_addition_order = [
								best_addition_order,
								copy.deepcopy(addition_order)]


	if no_id_count == min([args.replicates, len(array_choices)]):
		sys.exit("\nERROR:\n\nUnable to construct any trees as at least one "
			"specified array shares no spacers with any other array already "
			"in the tree. \n If you have arrays sharing very few spacers and "
			"did not use enough replicates to explore all tree space, then "
			"consider retrying with more replicates.\n Otherwise, check that "
			"you didn't mistype the array IDs to align.")

	sys.stderr.write("\nScore of best tree is: {}\n".format(best_score))

	if args.branch_support:
		leaf_bits_dict = {k:v for k,v in zip(
			labels, [i for i in range(len(labels))]
			)}
		tree_list_sub = [t[1] for t in tree_list]
		tree_list_polytomy = [
			tree_operations.resolve_polytomies(t) for t in copy.deepcopy(tree_list_sub)
			]

		clade_bins = [
			tree_operations.get_binary_nodes(
					t, leaf_bits_dict
				) for t in tree_list_polytomy
			]

		if isinstance(best_tree, list):
			for good_tree in best_tree:
				tree_operations.resolve_polytomies(good_tree)
				tree_operations.calculate_branch_support(
					good_tree,
					clade_bins,
					leaf_bits_dict)
		else:
			tree_operations.resolve_polytomies(best_tree)
			tree_operations.calculate_branch_support(
				best_tree,
				clade_bins,
				leaf_bits_dict)

	# First check how many spacers will need to be coloured

	occurrences = defaultdict(int)
	for array in array_choices[0]:
		for spacer in set(array.spacers):
			occurrences[spacer]+=1
	non_singleton_spacers = [
		spacer for spacer, count in occurrences.items() if count >1]
		
	# Process colour file related args
	spacer_cols_dict = colour_schemes.process_colour_args(
		args, non_singleton_spacers)

	try:
		if isinstance(best_tree, list):
			sys.stderr.write(
				"{} equivalantly parsimonious trees were identified.".format(
					len(best_tree)))
			for n, good_tree in enumerate(best_tree):
				order = [i.id for i in best_addition_order[n]]
				sys.stderr.write(
					"\nThe addition order to make the following tree was: "
					"{}\n\n".format(" ".join(order)))
				if args.print_tree:
					print(good_tree.as_ascii_plot(
						show_internal_node_labels=True))
					print('\n')
				if args.branch_support:
					out = dendropy.utility.textprocessing.StringIO()
					print(good_tree.seed_node.write_newick_bs(out))
				else:
					print(good_tree.as_string("newick"))
				if args.out_file:
					# Reset events in root array
					best_arrays[n] = reset_anc_mods(good_tree, best_arrays[n])
					filename = "{}_{}.png".format(args.out_file[:-4], n+1)
					sys.stderr.write(
						"Saving image of tree with array diagrams to "
						"{}\n".format(filename))
					plotting.plot_tree(tree=good_tree,
						array_dict=best_arrays[n], filename=filename,
						spacer_cols_dict=spacer_cols_dict,
						branch_lengths=args.branch_lengths,
						emphasize_diffs=args.emphasize_diffs, dpi=args.dpi,
						brlen_scale=args.brlen_scale,
						no_align_cartoons=args.no_align_cartoons,
						no_align_labels=args.no_align_labels,
						no_fade_ancestral=args.no_fade_anc,
						fig_h=args.plot_height, fig_w=args.plot_width,
						label_text_size=args.font_override_labels,
						annot_text_size=args.font_override_annotations,
						branch_support=args.branch_support)

				if args.output_arrays:
					filename = "{}_{}.txt".format(args.output_arrays[:-4], n+1)
					print("Saving details of arrays to {}\n".format(filename))
					with open(filename, 'w') as fout:
						fout.write('\n'.join(["{}\t{}".format(k," ".join(v.spacers)) for k,v in best_arrays[n].items()]))
		else:
			order = [i.id for i in best_addition_order]
			sys.stderr.write(
				"\nThe addition order to make the best tree was: "
				"{}\n\n".format(" ".join(order)))
			if args.print_tree:
				print(best_tree.as_ascii_plot(show_internal_node_labels=True))
				print("\n")
			if args.branch_support:
				out = dendropy.utility.textprocessing.StringIO()
				print(best_tree.seed_node.write_newick_bs(out))
			else:
				print(best_tree.as_string("newick"))
			if args.out_file:
				# Reset events in root array
				best_arrays = reset_anc_mods(best_tree, best_arrays)
				sys.stderr.write(
					"Saving image of tree with array diagrams to "
					"{}\n".format(args.out_file))
				plotting.plot_tree(tree=best_tree,
					array_dict=best_arrays, filename=args.out_file,
					spacer_cols_dict=spacer_cols_dict,
					branch_lengths=args.branch_lengths,
					emphasize_diffs=args.emphasize_diffs, dpi=args.dpi,
					brlen_scale=args.brlen_scale,
					no_align_cartoons=args.no_align_cartoons,
					no_align_labels=args.no_align_labels,
					no_fade_ancestral=args.no_fade_anc,
					fig_h=args.plot_height, fig_w=args.plot_width,
					label_text_size=args.font_override_labels,
					annot_text_size=False,
					branch_support=args.branch_support)
			if args.output_arrays:
				sys.stderr.write(
					"Saving details of arrays to "
					"{}\n".format(args.output_arrays))
				file_handling.write_array_file(best_arrays, args.output_arrays)
	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		exc_tb.print_exception()

	end_time = time.time()
	time_taken = timedelta(seconds=end_time-start_time)

	sys.stderr.write("\nTotal run time: {}\n".format(time_taken))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their relationship."
		)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()
	main(sys.argv[1:])

