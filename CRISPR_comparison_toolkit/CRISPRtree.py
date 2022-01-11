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

from cctk import array_parsimony, sequence_operations, tree_operations


def infer_ancestor(array1, array2, all_arrays, node_ids, node_count, existing_ancestor):
	"""
	Based on the modules found in two aligned arrays, construct a hypothethetical ancestral state of the two arrays
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): A list of names to be used to name internal nodes.
		node_count (int): The index of the name to use in the node_ids list to name this internal node.
		existing_ancestor(Array class instance): The existing ancestor array if there is one
	
	Returns:
		(str or list) A hypothesis of the ancestral state of the provided sequences.
	"""
	
	ancestor = array_parsimony.Array(node_ids[node_count], extant=False)

	array1, array2 = array_parsimony.find_modules(array1, array2)	

	# Process modules to build hypothetical ancestor
	idx = 0
	while idx <= max(array1.module_lookup.keys()):
		mod1 = array1.module_lookup[idx]
		mod2 = array2.module_lookup[idx]
		if mod1.type == "acquisition" or mod2.type == "acquisition": 
			# If either has acquired spacers not in the other then those weren't in the ancestor.
			idx = max([mod1.indices[-1], mod2.indices[-1]]) + 1 # Skip the rest of this module.
			continue
		else:
			if mod1.type == "shared":
				# If spacers are shared in these they are assumed to have been in the ancestor.
				ancestor.modules.append(mod1)
				idx = mod1.indices[-1] + 1
				continue
			elif mod1.type == "duplication":
				# If one has duplicated spacers but the other doesn't then he ancestral state likely didn't have them
				idx = max([mod1.indices[-1], mod2.indices[-1]]) + 1 # Skip the rest of this module.
				continue
			elif mod1.type == "trailer_loss":
				# If one lost a trailer end spacer then the ancestor had it.
				if not '-' in mod1.spacers: # Figure out which module has actual spacers.
					ancestor.modules.append(mod1)
				else:
					ancestor.modules.append(mod2)
				idx += 1
				continue
			elif mod1.type == "indel_gap" or mod1.type == "indel_mm":
				# If one array has just spacers and the other just gaps in the indel then pick the spacers as ancestral.
				if all([i == '-' for i in mod1.spacers]) or all([i == '-' for i in mod2.spacers]):
					# Figure out which one is all gaps and store the other one if it isn't singletons.
					if all([i == '-' for i in mod1.spacers]):
						# Check first if it's a duplication. If so, remove from ancestor.
						count_list = []
						for sp in mod2.spacers:
							if sp != '-':
								count = 0
								for spacer in mod2.spacers:
									if sp == spacer:
										count += 1
								count_list.append(count)
						if all([i>2 for i in count_list]): # If all of the spacers are present twice then duplication has occured in this region.
							idx = mod1.indices[-1] + 1 # Skip this module
							continue
						else:
							ancestor.modules.append(mod2)
					else:
						# Check first if it's a duplication. If so, remove from ancestor.
						count_list = []
						for sp in mod1.spacers:
							if sp != '-':
								count = 0
								for spacer in mod1.spacers:
									if sp == spacer:
										count += 1
								count_list.append(count)
						if all([i>2 for i in count_list]): # If all of the spacers are present twice then duplication has occured in this region.
							idx = mod1.indices[-1] + 1 # Skip this module
							continue
						else:
							ancestor.modules.append(mod1)
						
					idx = mod1.indices[-1] + 1
				else:
					# If both modules contain spacers then this may be two deletions from a larger ancestor or a recombination event.
					# If there is an existing ancestor to consider, use it to resolve this event
					if existing_ancestor:
						# First check if both modules exist in the existing ancestral array. If they do then keep the portion of the ancestor that they align with.
						if all([s in existing_ancestor.spacers for s in mod1.spacers if s != '-']) and all([s in existing_ancestor.spacers for s in mod2.spacers if s != '-']):
							existing_ancestor_indices = sorted([existing_ancestor.spacers.index(s) for s in mod1.spacers + mod2.spacers if s != '-'])
							new_module = array_parsimony.SpacerModule()
							new_module.spacers = existing_ancestor.spacers[existing_ancestor_indices[0]:existing_ancestor_indices[-1]+1]
							ancestor.modules.append(new_module)
							idx = mod1.indices[-1] + 1 
							continue

						#Otherwise check if either module exist in the existing ancestral array. If one does then keep that one.
						elif all([s in existing_ancestor.spacers for s in mod1.spacers if s != '-']):
							new_module = array_parsimony.SpacerModule()
							new_module.spacers = [s for s in mod1.spacers if s != '-']
							ancestor.modules.append(new_module)
							idx = mod1.indices[-1] + 1
							continue

						elif all([s in existing_ancestor.spacers for s in mod2.spacers if s != '-']):
							new_module = array_parsimony.SpacerModule()
							new_module.spacers = [s for s in mod2.spacers if s != '-']
							ancestor.modules.append(new_module)
							idx = mod2.indices[-1] + 1
							continue

						else: # If either has any spacers in the ancestor pick the one with the most
							if any([s in existing_ancestor.spacers for s in mod1.spacers+mod2.spacers if s != '-']):
								n1 = len([s in existing_ancestor.spacers for s in mod1.spacers if s != '-'])
								n2 = len([s in existing_ancestor.spacers for s in mod2.spacers if s != '-'])
								if n1 > n2:
									new_module = array_parsimony.SpacerModule()
									new_module.spacers = [s for s in mod1.spacers if s != '-']
									ancestor.modules.append(new_module)
									idx = mod1.indices[-1] + 1
									continue
								else:
									new_module = array_parsimony.SpacerModule()
									new_module.spacers = [s for s in mod2.spacers if s != '-']
									ancestor.modules.append(new_module)
									idx = mod2.indices[-1] + 1
									continue

					# Next check if all these spacers exits in another array
					spacers_to_check = mod1.spacers + mod2.spacers
					found_array = False # If we find all the spacers in an array keep it to determine order
					for array in all_arrays:
						count = 0
						for spacer in spacers_to_check:
							if spacer in array:
								count += 1
						if count == len(spacers_to_check):
							found_array = array
							continue
					if found_array: # All the spacers were found in 1 array. That looks like 2 deletions from larger array.
						# Were the spacers consecutive in another array
						if " ".join(mod1.spacers + mod2.spacers) in " ".join(found_array):
							new_mod = array_parsimony.SpacerModule()
							new_mod.spacers = mod1.spacers + mod2.spacers
							ancestor.modules.append(new_mod)
							idx = mod1.indices[-1] + 1
						elif " ".join(mod2.spacers + mod1.spacers) in " ".join(found_array):
							new_mod = array_parsimony.SpacerModule()
							new_mod.spacers = mod2.spacers + mod1.spacers
							ancestor.modules.append(new_mod)
							idx = mod1.indices[-1] + 1
						else:
							# Something else has happened. For now don't put these spacers in the ancestral array. May revisit.
							idx = mod1.indices[-1] + 1
					else: # Not all the spacers were found in a single array. Is either set found in another array?
						spacers1_found = False
						spacers2_found = False
						for array in all_arrays:
							count = 0
							for spacer in mod1.spacers:
								if spacer in array:
									count += 1
							if count == len(mod1.spacers):
								spacers1_found = True
								continue
						for array in all_arrays:
							count = 0
							for spacer in mod2.spacers:
								if spacer in array:
									count += 1
							if count == len(mod2.spacers):
								spacers2_found = True
								continue
						if spacers1_found and spacers2_found: 
							# Both are also in another array. Each has parsimony cost. Pick one.
							ancestor.modules.append(mod1)
							idx = mod1.indices[-1] + 1
						elif spacers1_found:
							ancestor.modules.append(mod1)
							idx = mod1.indices[-1] + 1
						elif spacers2_found:
							ancestor.modules.append(mod2)
							idx = mod2.indices[-1] + 1
						else:
							# Neither was found.
							idx = mod2.indices[-1] + 1
							
			



	ancestor.spacers = [i.spacers for i in ancestor.modules]
	if isinstance(ancestor.spacers[0], list):
		ancestor.spacers = [spacer for sublist in ancestor.spacers for spacer in sublist]

	return ancestor


def resolve_pairwise_parsimony(array1, array2, all_arrays, array_dict, node_ids, node_count, tree, event_costs):
	"""
	Given two arrays, make a hypothetical ancestral state and calculate parsimony distance of each input array to that ancestor. 
	Can only build an ancestor state if there are shared spacers so throws an error if none are found.
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		array_dict (dict): Dict of Array class instances with information about the nodes of your tree.
		node_ids (list): A list of names to be used to name internal nodes.
		node_count (int): The index of the name to use in the node_ids list to name this internal node.
		tree (Deondropy Tree class instance): The tree in which the arrays are located.
		event_costs (dict): Dict to look up event types and their parsimony costs.

	Returns:
		(tuple of Array class instances) The input Array class instances with module and distance info added, the ancestral array Array class instance. Tuple order is (array1, array2, ancestor)

	Raises:
		No ID: Raises an exception if the two arrays share no spacers.
	"""

	if len(list(set(array1.spacers) & set(array2.spacers))) > 0:

		array1, array2 = array_parsimony.find_modules(array1, array2)

		if len(tree) > 1:
			if tree.find_node_with_taxon_label(array1.id):
				try:
					existing_ancestor = array_dict[tree.find_node_with_taxon_label(array1.id).parent_node.taxon.label]
				except: #root is parent
					existing_ancestor = False
			elif tree.find_node_with_taxon_label(array2.id):
				try:
					existing_ancestor = array_dict[tree.find_node_with_taxon_label(array2.id).parent_node.taxon.label]
				except: #root is parent
					existing_ancestor = False
			else:
				existing_ancestor = False
		else:
			existing_ancestor = False


		array1.reset() # Make sure the distance and events haven't carried over from a previous usage
		array2.reset()
		
		ancestor = infer_ancestor(array1, array2, all_arrays, node_ids, node_count, existing_ancestor)

		array1 = array_parsimony.count_parsimony_events(array1, ancestor, array_dict, tree, True)
		array2 = array_parsimony.count_parsimony_events(array2, ancestor, array_dict, tree, True)

		for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
			array1.distance += array1.events[k] * v
			array2.distance += array2.events[k] * v

		return array1, array2, ancestor
	

	else:
		return "No_ID"


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
		comparator_array = array_parsimony.count_parsimony_events(comparator_array, array, array_dict, tree, False)
		array.reset()
		array = array_parsimony.count_parsimony_events(array, comparator_array, array_dict, tree, True)
		for k,v in event_costs.items():
			comparator_array.distance += comparator_array.events[k] * v
			array.distance += array.events[k] * v

		if any([i.type == "shared" for i in comparator_array.modules]):
			# if comparator_array.distance < best_score:
			if array.distance < best_score:
				# best_score = comparator_array.distance
				best_score = array.distance
				best_match = comparator_array
			elif comparator_array.distance == best_score:
				if best_match.extant and not comparator_array.extant:
					# Prefer to join arrays to existing ancestors rather than to leaves if they are the same distance.
					best_match = comparator_array
	if not best_match:
		return "No_ID"


	if not any([i.type == "shared" for i in best_match.modules]):
		return "No_ID"


	else:
		return best_match


def replace_existing_array(existing_array, new_array, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed, event_costs, tree_namespace):
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
	results = resolve_pairwise_parsimony(existing_array, new_array, all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	
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
			tree_child_dict[a.id] = dendropy.Node(edge_length=a.distance)
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


def build_tree_single(arrays, tree_namespace, score, all_arrays, node_ids, event_costs):
	"""
	Search treespace for most parsimonious tree using single process.
	Args:
		arrays (list): Ordered list of Array class instances of the arrays to analyse. Will be added to the tree in the provided order.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
		score (int): The score to beat. If at any point during tree construction the tree has total branch lengths above this score construction will be aborted
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): The names for internal nodes in the tree to be assigned.
		event_costs (dict): Dict to look up event types and their parsimony costs.
	
	Returns:
		(tuple) Returns array_dict and dendropy.tree object if the tree beats or equals the provide score. Else retuns tuple of (False, False).
	"""
	initial_arrays = copy.deepcopy(arrays)
	tree = dendropy.Tree(taxon_namespace=tree_namespace)
	array_dict = {}
	tree_child_dict = {}
	node_count = 0 # Keep track of which internal node ID should be used for each node
	# Remove the arrays being compared from checks to see what spacers are in other arrays.
	# That way we only worry about ancestral states accounting for deletion of spacers in arrays not yet added to the tree.
	all_arrays = [a for a in all_arrays if a not in [arrays[0].spacers, arrays[1].spacers]]
	results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)

	while results == "No_ID":
		arrays.append(arrays[1])
		del arrays[1]
		results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays,  array_dict, node_ids, node_count, tree, event_costs)
	node_count += 1
	array1, array2, ancestor = results
	for a in [array1, array2, ancestor]:
		# Create tree nodes
		tree_child_dict[a.id] = dendropy.Node(edge_length=a.distance)
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


			# Recheck child - ancestor branch length to find indels that would have to occur multiple times
			for node in tree:
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
			if brlen > score:
				return (False, False, False)

					
			if a == arrays[-1]:
				tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works
				return (array_dict, tree, initial_arrays)


	else:
		return (array_dict, tree, initial_arrays)


def build_tree_multi(arrays, tree_namespace, all_arrays, node_ids, event_costs):
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

	results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	while results == "No_ID":

		arrays.append(arrays[1])
		del arrays[1]
		results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
	node_count += 1
	array1, array2, ancestor = results

	for a in [array1, array2, ancestor]:
		# Create tree nodes
		tree_child_dict[a.id] = dendropy.Node(edge_length=a.distance)
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
			for node in tree:
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

			tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works
			return array_dict, tree, initial_arrays


def parse_args():
	parser = argparse.ArgumentParser(
		description="Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship."
		)
	parser.add_argument(
		"-a", dest="array_file", required = True,
		help="Specify array representatives file."
		)
	parser.add_argument(
		"-p", dest="print_tree", action='store_true',  
		help="Print a graphical representation of the tree using ascii characters."
		)
	parser.add_argument(
		"-x", dest="fix_order", action='store_true',  
		help="Only build one tree using the provided order of arrays. Good for recreating previously found trees."
		)
	parser.add_argument(
		"-e", dest="emphasize_diffs", action='store_true',  
		help="When plotting a representation of the tree with cartooned arrays, emphasize locations where arrays differ from their hypothetical ancestor."
		)
	parser.add_argument(
		"-b", dest="branch_lengths", action='store_true', 
		help="When plotting a representation of the tree Include branch lengths at midpoint on branches. N.B. This does not control inclusion of branch lengths in newick format tree output, which are always included."
		)
	parser.add_argument(
		"-g", dest="no_align_cartoons", action='store_true', default=False,
		help="When plotting a representation of the tree with cartooned arrays, this option controls whether those cartoons are drawn next to the corresponding node in the tree. By default cartoons will be aligned at the trailer end."
		)
	parser.add_argument(
		"-c", dest="no_align_labels", action='store_true', default=False,
		help="Should node labels be placed at the corresponding internal node or leaf. By default they are all aligned. N.B. For this setting to work, cartoons must also placed at nodes using the -g option."
		)
	parser.add_argument(
		"-r",  dest="replicates", type=int, nargs="?", default = 1,
			help="Specify number of replicates of tree building to perform. The more replicates, the greater the chance that a better tree will be found. Default: 1"
		)
	parser.add_argument(
		"-q",  dest="acquisition", type=int, nargs="?", default = 1,
			help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
		)
	parser.add_argument(
		"-i",  dest="indel", type=int, nargs="?", default = 10,
			help="Specify the parsimony cost of an indel event involving one or more spacers. Default: 10"
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
		"-n",  dest="no_ident", type=int, nargs="?", default = 100,
			help="Specify the parsimony cost of a tree in which an array is predicted to have descended from another array with which it shares no spacers. Default: 100"
		)
	parser.add_argument(
		"-o", dest="output_tree", required = False,
		help="Specify filename for the graphical representation of your tree with hypothetical intermediate arrays as a png."
		)
	parser.add_argument(
		"-t",  dest="num_threads", type=int, nargs="?", default = 1,
			help="Specify number of threads to use for building trees. Using multiple threads will speed up the search for trees when performing many replicates with the -r option. Default: 1"
		)
	parser.add_argument(
		"-y", dest="output_arrays", required = False,
		help="Specify filename for the details of you final arrays with hypothetical intermediate arrays in the same format as your input array_file. If there are multiple best trees, one file will be created per tree numbered in the order they are described in the stdout output."
		)
	parser.add_argument(
		"-m", dest="colour_scheme_outfile", required = False, 
		help="Specify output file to store json format dictionary of the colour schemes used for spacers in this run."
		)
	parser.add_argument(
		"-s", dest="colour_scheme_infile", required = False, 
		help="Specify input file containing json format dictionary of the colour scheme to be used for spacers in this run. Any spacers not in the input file will be coloured according to the normal process."
		)
	parser.add_argument(
		"--seed", dest="seed", type=int, required = False, default = 2,
		help="The order of outline and fill colours assigned to spacers is semi-random. Change it by providing a number here to change which colours are assigned to each spacer."
		)
	parser.add_argument(
		"--dpi", dest="dpi", type=int, required = False, default = 600,
		help="The desired resolution of the output image."
		)
	parser.add_argument(
		"--fade_anc", dest="fade_ancestral", action='store_true', 
		help="De-emphasize ancestral array cartoons using transparency. This option helps to make it clear which are the inferred ancestral arrays and which are the arrays being analyzed. However, this option reduces colour contrast between spacers."
		)
	parser.add_argument(
		"arrays_to_join", nargs="*",  
		help="Specify the IDs of the arrays you want to join. If none provided, joins all arrays in the provided array representatives file. **If given, must come at the end of your command after all other arguments.**"
		)

	args = parser.parse_args()

	if args.no_align_labels and not args.no_align_cartoons:
		print("\n\nYou chose settings that would align node labels but not array cartoons in the output tree image. Labels cannot be aligned without also aligning array cartoons so that setting has been overwritten. Both will be aligned.\n\n")
		args.no_align_cartoons = True

	return args


def main():

	start_time = time.time()

	args = parse_args()

	print("Running with the following command:\n{}\n".format(" ".join(sys.argv)))

	event_costs = { 
					"acquisition" : args.acquisition,
					"indel" : args.indel,
					"repeated_indel" : args.rep_indel,
					"duplication": args.duplication,
					"trailer_loss": args.trailer_loss,
					"no_ident": args.no_ident
					}

	# hex values from this website http://phrogz.net/css/distinct-colors.html

	Cols_hex_27 = ["#fd5925", "#dbc58e", "#008d40", "#304865", "#934270", "#f7b8a2", "#907500", "#45deb2", "#1f4195", "#d67381", "#8e7166", "#afb200", "#005746", "#a598ff", "#8f0f1b", "#b96000", "#667f42", "#00c7ce", "#9650f0", "#614017", "#59c300", "#1a8298", "#b5a6bd", "#ea9b00", "#bbcbb3", "#00b0ff", "#cd6ec6"]


	Cols_tol = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"]

	Cols_hex_12 = ["#07001c", "#ff6f8d", "#4c62ff", "#92ffa9", "#810087", "#bcffe6", "#490046", "#00c8ee", "#b53900", "#ff8cf7", "#5b5800", "#14d625"]

	# Generate strings to assign as internal node_IDs (This makes 702)

	node_ids = ["Anc " + i for i in ascii_uppercase]
	if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
		node_ids += ["Anc " + "".join(i) for i in product(ascii_uppercase, repeat=(ceil(log(len(args.arrays_to_join), 26))))]


	array_spacers_dict = {}
	with open(args.array_file, 'r') as fin:
		for line in fin.readlines():
			bits = line.split()
			array_spacers_dict[bits[0]] = bits[2:]

	if args.arrays_to_join:
		arrays = [array_parsimony.Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
		labels = args.arrays_to_join
	else:
		arrays = [array_parsimony.Array(i, array_spacers_dict[i]) for i in array_spacers_dict.keys()]
		labels = list(array_spacers_dict.keys())

	all_arrays = [array.spacers for array in arrays]

	# Check if any arrays share no spacers with any others. If so, exit with error message.
	all_spacers = [spacer for array in all_arrays for spacer in set(array)]
	all_non_singleton_spacers = set([spacer for spacer, count in Counter(all_spacers).items() if count >1])
	no_id = False
	no_id_arrays = []
	for array in arrays:
		if len(set(array.spacers).intersection(all_non_singleton_spacers)) == 0:
			no_id = True
			no_id_arrays.append(array.id)
	if no_id:
		print("No shared spacers were found between array(s) {} and the other arrays in this dataset. Please ensure that all arrays share at least 1 spacer with another array. If there are no shared spacers between 1 or a group of arrays and another group of arrays then this script will fail".format(", ".join(no_id_arrays)))
		sys.exit()
		
	if len(labels) < 9:
		array_choices = [list(i) for i in permutations(arrays, len(arrays))]
		random.shuffle(array_choices)



		if len(array_choices) > args.replicates:
			print("\nThere are {} possible trees to check. If you want to check every possible tree then set -r {}\n".format(len(array_choices), len(array_choices)))
			array_choices = [array_choices[i] for i in range(args.replicates)]

		elif len(array_choices) < args.replicates:
			print("\nThere are only {} possible trees to check. You specified a greater number of replicates than there are possible trees. All possible trees will be checked.\n".format(len(array_choices)))

		else:
			print("\nYou specified a number of replicates equal to the number of possible trees. All possible trees will be checked.\n")
	else:
		array_choices = [random.sample(arrays, len(arrays)) for i in range(args.replicates)]

	taxon_namespace = dendropy.TaxonNamespace(labels + node_ids)

	best_score = 99999999

	# Keep track of how many times no common spacers are found.
	# If no common spacers are ever found then print an error message saying so

	no_id_count = 0 

	num_threads = args.num_threads if not args.fix_order else 1

	if len(array_choices) < num_threads:
		num_threads = 1

	if num_threads == 1:
		for i in range(min([args.replicates, len(array_choices)])):
			if args.fix_order:
				if not args.arrays_to_join:
					print("You must provide the order you want to fix when using the fixed order option!\n\nABORTING.")
					sys.exit()
				addition_order = [array_parsimony.Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
			else:
				addition_order = array_choices[i]

			if args.fix_order:
				if not args.arrays_to_join:
					print("You must provide the order you want to fix when using the fixed order option!\n\nABORTING.")
					sys.exit()
				addition_order = [array_parsimony.Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
			else:
				addition_order = array_choices[i]
			try:
				array_dict, tree, addition_order = build_tree_single(addition_order, taxon_namespace, best_score, all_arrays, node_ids, event_costs)
			except Exception as e:
				print("Failed while trying to build the tree with the following array order:\n{}\n\nError:\n{}".format(" ".join([i.id for i in addition_order]), e))
				exc_type, exc_obj, exc_tb = sys.exc_info()
				exc_tb.print_exception()
				sys.exit()

			if array_dict:
				score = tree.length()
				if score < best_score:
					best_arrays = copy.deepcopy(array_dict) # Keep inferred ancestral states and information
					best_score = copy.deepcopy(score)
					best_addition_order = copy.deepcopy(addition_order)
					# Keep one copy for comparisons as copy.deepcopy makes a new taxon namespace which breaks comparisons.
					best_tree_comparator = dendropy.Tree(tree)
					best_tree = copy.deepcopy(tree)
				elif score == best_score:
					if isinstance(best_tree, list):
						# Check this tree isn't identical to one that's already been found
						if not any([
							dendropy.calculate.treecompare.weighted_robinson_foulds_distance(good_tree, tree) == 0. for good_tree in best_tree_comparator
							]):
							best_tree_comparator.append(dendropy.Tree(tree))
							best_tree.append(copy.deepcopy(tree))
							best_arrays.append(copy.deepcopy(array_dict))
							best_addition_order.append(copy.deepcopy(addition_order))
					else:
						# Check this tree isn't identical to the one that's already been found
						if dendropy.calculate.treecompare.weighted_robinson_foulds_distance(best_tree_comparator, tree) != 0.:
							best_tree_comparator = [best_tree_comparator, dendropy.Tree(tree)]
							best_tree = [best_tree, copy.deepcopy(tree)]
							best_arrays = [best_arrays, copy.deepcopy(array_dict)]
							best_addition_order = [best_addition_order, copy.deepcopy(addition_order)]
				if args.fix_order:
					break
			else:
				no_id_count += 1


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
				all_tree_comparators.read(data=tree.as_string("newick"), schema="newick")
				score = tree.length()
				if score < best_score:
					best_arrays = copy.deepcopy(array_dict) # Keep inferred ancestral states and information
					best_score = copy.deepcopy(score)
					best_addition_order = copy.deepcopy(addition_order)
					# Keep one copy for comparisons as copy.deepcopy makes a new taxon namespace which breaks comparisons.
					best_tree_comparator = all_tree_comparators[-1:]
					best_tree = copy.deepcopy(tree)
				elif score == best_score:
					if isinstance(best_tree, list):
						# Check this tree isn't identical to one that's already been found
						if not any([
							dendropy.calculate.treecompare.weighted_robinson_foulds_distance(good_tree, all_tree_comparators[-1]) == 0. for good_tree in best_tree_comparator
							]):
							best_tree_comparator.append(all_tree_comparators[-1])
							best_tree.append(copy.deepcopy(tree))
							best_arrays.append(copy.deepcopy(array_dict))
							best_addition_order.append(copy.deepcopy(addition_order))
					else:
						# Check this tree isn't identical to the one that's already been found
						if dendropy.calculate.treecompare.weighted_robinson_foulds_distance(best_tree_comparator[0], all_tree_comparators[-1]) != 0.:
							best_tree_comparator.append(all_tree_comparators[-1])
							best_tree = [best_tree, copy.deepcopy(tree)]
							best_arrays = [best_arrays, copy.deepcopy(array_dict)]
							best_addition_order = [best_addition_order, copy.deepcopy(addition_order)]


	if no_id_count == min([args.replicates, len(array_choices)]):
		print("\nERROR:\n\nUnable to construct any trees as at least one specified array\
	shares no spacers with any other array already in the tree. \n\
	If you have arrays sharing very few spacers and did not use enough replicates\
	to explore all tree space, then consider retrying with more replicates.\n\
	Otherwise, check that you didn't mistype the array IDs to align.\n\
	If not one of these then that may indicate a problem with the program. \
	Please send the data that led to this error to Alan so he can try to figure it out.\n\n")
		sys.exit()


	print("\nScore of best tree is: {}\n".format(best_score))


	# First check how many spacers will need to be coloured

	all_spacers = []
	for array in array_choices[0]:
		all_spacers += array.spacers
	non_singleton_spacers = [spacer for spacer, count in Counter(all_spacers).items() if count >1]
	if len(non_singleton_spacers) > 8:
		if len(non_singleton_spacers) > 12: 
			if len(non_singleton_spacers) > 27:
				print("{} spacers found in multiple arrays. Using fill and outline colour combinations to distinguish spacers.".format(len(non_singleton_spacers)))
				if len(non_singleton_spacers) < 65:
					col_scheme = Cols_tol
				elif len(non_singleton_spacers) < 145:
					col_scheme = Cols_hex_12
				else:
					col_scheme = Cols_hex_27
				
				random.seed(args.seed)
				combos = [i for i in permutations(col_scheme, 2)]
				combos += [(i,i) for i in col_scheme]
				colours = random.sample(combos, len(combos))
				random.seed(None)
				# colours = []
				# for i in range((len(non_singleton_spacers)+len(col_scheme)-1)//len(col_scheme)): # Repeat the same colour scheme.
				# 	for j in col_scheme:
				# 		colours += [(j, col_scheme[i])]
			else:
				colours = [(i, "#000000") for i in Cols_hex_27]
		else:
			colours = [(i, "#000000") for i in Cols_hex_12]
	else:
		colours = [(i, "#000000") for i in Cols_tol]
	# if provided, read in colour scheme.
	spacer_cols_dict = {}
	if args.colour_scheme_infile:
		with open(args.colour_scheme_infile, 'r') as fin:
			user_spacer_cols_dict = json.load(fin)
		for s in non_singleton_spacers:
			if s in user_spacer_cols_dict.keys():
				spacer_cols_dict[s] = user_spacer_cols_dict[s]
			else:
				colour_idx = 0
				while list(colours[colour_idx]) in spacer_cols_dict.values():
					colour_idx += 1
				spacer_cols_dict[s] = colours[colour_idx]
				colour_idx += 1
	else:
		# build a dictionary with colours assigned to each spacer.
		for i, spacer in enumerate(sorted(non_singleton_spacers)):
			spacer_cols_dict[spacer] = colours[i]


	if args.colour_scheme_outfile:
		with open(args.colour_scheme_outfile, 'w', encoding='utf-8') as fout:
			json.dump(spacer_cols_dict, fout, ensure_ascii=False, indent=4)

	try:

		if isinstance(best_tree, list):
			print("{} equivalantly parsimonious trees were identified.".format(len(best_tree)))
			for n, good_tree in enumerate(best_tree):
				order = [i.id for i in best_addition_order[n]]
				print("\nThe addition order to make the following tree was: {}\n".format(" ".join(order)))
				if args.print_tree:
					print(good_tree.as_ascii_plot(show_internal_node_labels=True))
					print('\n\n')
				print(good_tree.as_string("newick"))
				print('\n\n')
				if args.output_tree:
					filename = "{}_{}.png".format(args.output_tree[:-4], n+1)
					print("Saving image of tree with array diagrams to {}\n".format(filename))
					tree_operations.plot_tree(tree=good_tree,
						array_dict=best_arrays[n], filename=filename,
						spacer_cols_dict=spacer_cols_dict,
						branch_lengths=args.branch_lengths,
						emphasize_diffs=args.emphasize_diffs, dpi=args.dpi,
						no_align_cartoons=args.no_align_cartoons,
						no_align_labels=args.no_align_labels,
						fade_ancestral=args.fade_ancestral, fig_h=6, fig_w=6,
						label_text_size=10, annot_text_size=5
						# brlen_scale=0.5, font_scale=1, fig_h=1, fig_w=1,
						# label_text_size=False, annot_text_size=False
						)

				if args.output_arrays:
					filename = "{}_{}.txt".format(args.output_arrays[:-4], n+1)
					print("Saving details of arrays to {}\n".format(filename))
					with open(filename, 'w') as fout:
						fout.write('\n'.join(["{}\t{}".format(k," ".join(v.spacers)) for k,v in best_arrays[n].items()]))
		else:
			order = [i.id for i in best_addition_order]
			print("\nThe addition order to make the best tree was: {}\n\n".format(" ".join(order)))
			if args.print_tree:
				print(best_tree.as_ascii_plot(show_internal_node_labels=True))
			print(best_tree.as_string("newick"))
			if args.output_tree:
				tree_operations.resolve_polytomies(best_tree)
				print("Saving image of tree with array diagrams to {}\n".format(args.output_tree))
				tree_operations.plot_tree(tree=best_tree,
					array_dict=best_arrays, filename=args.output_tree,
					spacer_cols_dict=spacer_cols_dict,
					branch_lengths=args.branch_lengths,
					emphasize_diffs=args.emphasize_diffs, dpi=args.dpi,
					no_align_cartoons=args.no_align_cartoons,
					no_align_labels=args.no_align_labels,
					fade_ancestral=args.fade_ancestral, fig_h=6, fig_w=6,
					label_text_size=10, annot_text_size=5)
				# plot_tree_old(best_tree, best_arrays, args.output_tree, spacer_cols_dict, 
				# 	args.branch_lengths, args.emphasize_diffs, args.dpi, 
				# 	args.no_align_cartoons, args.no_align_labels, args.fade_ancestral)
			if args.output_arrays:
				print("Saving details of arrays to {}\n".format(args.output_arrays))
				with open(args.output_arrays, 'w') as fout:
					fout.write('\n'.join(["{}\t{}".format(k," ".join(v.spacers)) for k,v in best_arrays.items()]))
	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		exc_tb.print_exception()

	end_time = time.time()
	time_taken = timedelta(seconds=end_time-start_time)

	print("Total run time: {}".format(time_taken))

if __name__ == '__main__':
	main()

