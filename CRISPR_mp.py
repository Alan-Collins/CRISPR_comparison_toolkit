#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-7-7
# DESCRIPTION :  Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship.

import sys
import argparse
import numpy as np
from itertools import product
from string import ascii_lowercase
import random
from collections import Counter
import dendropy
import copy


parser = argparse.ArgumentParser(
	description="Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship."
	)
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"-p", dest="print_tree", action='store_true',  
	help="Print a graphical representation of the tree using ascii characters (required ete3 to be installed)."
	)
parser.add_argument(
	"arrays_to_join", nargs="+",  
	help="Specify the IDs of the arrays you want to join. **Must come at the end of your command after all other arguments.**"
	)


args = parser.parse_args(sys.argv[1:])

class Array():
	"""
	Class to store information about extant and inferred ancestral CRISPR arrays to aid in their comparisons.
	
	Attributes:
		id (str): Identifier for this array.
		extant (bool): A boolean indicating if the array is extant in our dataset or if it was inferred as a hypothetical ancestral state.
		modules (list): A list of the contiguous blocks of spacers with common features (e.g. consecutive spacers that are absent in aligned array).
		spacers (list): A list of the spacers in this array.
		aligned (list): A list of spacers in aligned format relative to another array
		module_lookup (dict): A dict with indices as keys and Spacer_Module instances as values where a given Spacer_module instance will be pointed to by all the indices at which it is located.
		distance (int): Parsimony distance from the hypothetical ancestral state of this array.
		events (dict): Record what each event is for later addition of event parsimony weights.
		"""
	def __init__(self, ID, spacers=[], extant=True):
		self.id = ID
		self.extant = extant
		self.modules = []
		self.spacers = spacers
		self.aligned = []
		self.module_lookup = {}
		self.distance = 0
		self.events = { 
						"acquisition" : 0,
						"indel" : 0,
						"duplication": 0,
						}

	def sort_modules(self):
		self.modules.sort(key=lambda x: int(x.indices[0]))

	def reset(self):
		self.distance = 0
		self.events = { 
						"acquisition" : 0,
						"indel" : 0,
						"duplication": 0,
						}


class Spacer_Module():
	"""
	Class to store information about spacers in CRISPR arrays.
	
	Attributes:
		type (str): A string indicating the nature of this module. Possible types:
			N.B. leader region ends once first identical spacer is found.
			- aqcuisition: Spacers at the leader end of one array where no (or different) spacers are found in the other array. 
			- no_acquisition: Gap at leader end of array when other array has spacers not found in this array.
			- indel: non-leader region with spacers present in one array but a gap or different spacers in the other array
			- shared: Region where both arrays have the same spacers.
		spacers (list): A list of the spacer IDs in this module.
		indices (list): A list of the indices in the respective array where this spacer module is located.
	"""
	def __init__(self):
		self.type = ""
		self.spacers = []
		self.indices = []
		

def needle(seq1, seq2, match = 100, mismatch = -1, gap = -2):
	"""
	Perform Needleman-Wunsch pairwise alignment of two sequences.
	Args:
		seq1 (str or list): First sequence of items to align.
		seq2 (str or list): Second sequence of items to align
		match (int): Score for match at a position in alignment.
		mismatch(int): Penalty for mismatch at a position in alignment.
		gap (int): Penalty for a gap at a position in alignment.
	
	Returns:
		(tuple of str lists) Returns a tuple containing the input seq1 and seq2 aligned with '-' added as gaps.
		If strings were given then strings are returned. If lists were given then lists are returned.

	"""

	# Make a list of lists of 0s with dimensions x by y: list containing x lists of y 0s each.
	grid = np.zeros((len(seq2)+1, len(seq1)+1))

	# Fill in grid with scores for all possible alignments
	# First score for no alignment (i.e. all gaps)
	for i in range(len(seq1)+1):
		grid[0][i] = gap*i
	for i in range(len(seq2)+1):
		grid[i][0] = gap*i

	# Then score for each cell if you came to it from the nearest best cell/
	for i in range(len(seq1)):
		for j in range(len(seq2)):
			if seq1[i] == seq2[j]:
				score = match
			else:
				score = mismatch
			grid[j+1][i+1] = max([grid[j][i]+score, grid[j+1][i]+gap, grid[j][i+1]+gap])

	i = len(seq2)
	j = len(seq1)

	# Read back through the grid along the best path to create the best alignment
	align1, align2 = [], []
	while i > 0 and j > 0: # end when it reaches the top or the left edge
		score_current = grid[i][j]
		score_diagonal = grid[i-1][j-1]
		score_up = grid[i][j-1]
		score_left = grid[i-1][j]
		if seq1[j-1] == seq2[i-1]:
			score = match
		else:
			score = mismatch
		# Check to figure out which cell the current score was calculated from,
		# then update i and j to correspond to that cell.
		if score_current == score_diagonal + score:
			align1.append(seq1[j-1])
			align2.append(seq2[i-1])
			i -= 1
			j -= 1
		elif score_current == score_up + gap:
			align1.append(seq1[j-1])
			align2.append('-')
			j -= 1
		elif score_current == score_left + gap:
			align1.append('-')
			align2.append(seq2[i-1])
			i -= 1

	# Finish tracing up to the top left cell
	while j > 0:
		align1.append(seq1[j-1])
		align2.append('-')
		j -= 1
	while i > 0:
		align1.append('-')
		align2.append(seq2[i-1])
		i -= 1
	
	# Since we traversed the score matrix backwards, need to reverse alignments.
	align1 = align1[::-1]
	align2 = align2[::-1]

	if isinstance(seq1, str) and isinstance(seq2, str):
		align1 = ''.join(align1)
		align2 = ''.join(align2)
	
	return align1, align2


def find_modules(array1, array2):
	"""
	Identify contiguous stretches of spacers that share an evolutionary property in terms of how they differ from the comparator array (e.g. indel region). Assigns types to each module according to the types defined in the Spacer_Module class docstring.
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
	
	Returns:
		(tuple of Array class instances) The provided arrays with module information added to the .module attribute.
	"""

	array1.aligned, array2.aligned = needle(array1.spacers, array2.spacers)

	# If this isn't the first comparison these arrays have been part of
	# then need to reset module list and lookup dict.
	array1.modules = []
	array2.modules = []

	array1.module_lookup = {}
	array2.module_lookup = {}

	leader = True
	gap = False
	mismatch = False

	module1 = Spacer_Module()
	module2 = Spacer_Module()

	# Identify modules in aligned arrays

	for n, (a,b) in enumerate(zip(array1.aligned, array2.aligned)):

		# Leader end processing

		if leader:
			if a == '-' or b == '-': # Indicates one array acquired spacers that the other didn't
				if a == '-':

					# Module1 processing for no_acqusition module
					if module1.type != "no_acquisition" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = Spacer_Module()
					module1.type = "no_acquisition"
					module1.indices.append(n)
					module1.spacers.append(a)

					# Module2 processing for acquisition module
					# No need to check module type as gap penalty means needle func will never return stretch like the following:
					# array1	A B C - - -
					# array2	- - - D E F
					module2.type = "acquisition"
					module2.indices.append(n)
					module2.spacers.append(b)

				else:
					if module2.type != "no_acquisition" and module2.type != "":
						array2.modules.append(module2)
						for k in module2.indices:
							array2.module_lookup[k] = module2
						module2 = Spacer_Module()
					module2.type = "no_acquisition"
					module2.indices.append(n)
					module2.spacers.append(b)

					module1.type = "acquisition"
					module1.indices.append(n)
					module1.spacers.append(a)

			elif a != b: # Mismatch at leader end probably means they've both acquired spacers since ancestor
				# Indel also possible. Maybe implement parsimony evaluation of that when comparing ancestor to other arrays?

				if module1.type != "acquisition" and module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = Spacer_Module()
				module1.type = "acquisition"
				module1.indices.append(n)
				module1.spacers.append(a)

				if module2.type != "acquisition" and module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = Spacer_Module()
				module2.type = "acquisition"
				module2.indices.append(n)
				module2.spacers.append(b)

			else: # Other alternative is identical spacers and end of leader region
				leader = False
				if module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = Spacer_Module()
				module1.type = "shared"
				module1.indices.append(n)
				module1.spacers.append(a)
				if module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = Spacer_Module()
				module2.type = "shared"
				module2.indices.append(n)
				module2.spacers.append(b)

		else:
			if a == b:
				if module1.type != "shared" and module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = Spacer_Module()
				module1.type = "shared"
				module1.indices.append(n)
				module1.spacers.append(a)

				if module2.type != "shared" and module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = Spacer_Module()
				module2.type = "shared"
				module2.indices.append(n)
				module2.spacers.append(b)
			else:
				# Module1 processing for indel module
				if module1.type != "indel" and module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = Spacer_Module()
				module1.type = "indel"
				module1.indices.append(n)
				module1.spacers.append(a)
				# Module2 processing for indel module
				if module2.type != "indel" and module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = Spacer_Module()
				module2.type = "indel"
				module2.indices.append(n)
				module2.spacers.append(b)


	if module1.type != "":
		array1.modules.append(module1)
		for k in module1.indices:
			array1.module_lookup[k] = module1
	if module2.type != "":
		array2.modules.append(module2)
		for k in module2.indices:
			array2.module_lookup[k] = module2
	array1.sort_modules()
	array2.sort_modules()

	return array1, array2


def infer_ancestor(array1, array2, all_arrays, node_ids, node_count):
	"""
	Based on the modules found in two aligned arrays, construct a hypothethetical ancestral state of the two arrays
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): A list of names to be used to name internal nodes.
		node_count (int): The index of the name to use in the node_ids list to name this internal node.
	
	Returns:
		(str or list) A hypothesis of the ancestral state of the provided sequences.
	"""

	ancestor = Array(node_ids[node_count], extant=False)

	array1, array2 = find_modules(array1, array2)	

	# Process modules to build hypothetical ancestor
	idx = 0
	while idx < max(array1.module_lookup.keys()):
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
			elif mod1.type == "indel":
				# If one array has just spacers and the other just gaps in the indel then pick the spacers as ancestral unless those spacers are singletons.
				if all([i == '-' for i in mod1.spacers]) or all([i == '-' for i in mod2.spacers]):
					# Figure out which one is all gaps and store the other one if it isn't singletons.
					if all([i == '-' for i in mod1.spacers]):
						# Check first if it's a duplication. If so, remove from ancestor.
						# If not a duplication check if these spacers are in another array. If so keep.
						duplication = False # Start with the assumption this isn't a duplication
						for sp in mod2.spacers:
							count = 0
							for spacer in mod2.spacers:
								if sp == spacer:
									count += 1
								if count == 2: # If the same spacer is present twice then duplication has occured in this region.
									duplication = True
									break
						if duplication == False:
							singleton = True # Start with the assumption that these spacers aren't found in another array.
							for sp in mod2.spacers:
								count = 0
								for array in all_arrays:
									if sp in array:
										count += 1
									if count == 2:
										singleton = False
										continue
							if singleton == False:
								ancestor.modules.append(mod2)
					else:
						duplication = False # Start with the assumption this isn't a duplication
						for sp in mod2.spacers:
							count = 0
							for spacer in mod2.spacers:
								if sp == spacer:
									count += 1
								if count == 2: # If the same spacer is present twice then duplication has occured in this region.
									duplication = True
									break
						if duplication == False:
							singleton = True # Start with the assumption that these spacers aren't found in another array.
							for sp in mod1.spacers:
								count = 0
								for array in all_arrays:
									if sp in array:
										count += 1
									if count == 2:
										singleton = False
										continue
							if singleton == False:
								ancestor.modules.append(mod1)
					idx = mod1.indices[-1] + 1
				else:
					# If both modules contain spacers then this may be two deletions from a larger ancestor or a recombination event.
					# First check if all these spacers exits in another array
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
							new_mod = Spacer_Module()
							new_mod.spacers = mod1.spacers + mod2.spacers
							ancestor.modules.append(new_mod)
						elif " ".join(mod2.spacers + mod1.spacers) in " ".join(found_array):
							new_mod = Spacer_Module()
							new_mod.spacers = mod2.spacers + mod1.spacers
							ancestor.modules.append(new_mod)
							idx = mod1.indices[-1] + 1
						else:
							# Something else has happened. For now don't put these spacers in the ancestral array. May revisit.
							pass
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


def find_dupes(child, ancestor):
	"""
	Look for spacers in two copies in the child that exist in one copy in the ancestor
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
	
	Returns:
		(list) list of indices of duplicated spacers if found. Empty list if none found
	"""

	spacer_counts = Counter(child.aligned)

	dupe_indices = []

	for k,v in spacer_counts.items():
		if k != '-':
			if v > 1:
				indices = [idx for idx, el in enumerate(child.aligned) if el == k]
				dupe_indices.append(indices)

	dupe_groups = []
	if len(dupe_indices) > 0: # If there were no duplications, return the empty list.
		# Identify groups of consecutive spacers to seperate duplicated modules in case there is more than 1 event.

		total_spacers = sum([len(i) for i in dupe_indices])
		included_spacers = [[],[]] # First list stores which spacer has been included. Second stores which copy of that spacer has been used.
		checked_spacers = []
		group = [min([x for y in dupe_indices for x in y])] # Start with the first duplicated spacer

		x = 0
		# while sum([len(i) for i in dupe_groups]) < total_spacers:
		while all([_ != [] for _ in dupe_indices]):
			for n, spacer in enumerate(dupe_indices):
				if spacer != [] and n not in included_spacers[0] and n not in checked_spacers:
					for m, idx in enumerate(spacer):
						if idx == group[-1]:
							included_spacers[0].append(n)
							included_spacers[1].append(m)
							break
						elif idx == group[-1] + 1: # If this spacer follows the last spacer added, add it.
							group.append(dupe_indices[n][m])
							checked_spacers = [] # If we change the latest spacer we need to recheck ones that haven't added yet
							included_spacers[0].append(n)
							included_spacers[1].append(m)
							break
						elif idx != group[-1] + 1 and idx == spacer[-1]: 
							checked_spacers.append(n) # If this spacer doesn't have an index following the last spacer then we don't need to check it again unless we find another spacer elsewhere.
			if len(checked_spacers) + len(included_spacers[0]) == len(dupe_indices):
				dupe_groups.append(group)
				checked_spacers = []
				for x, y in zip(reversed(included_spacers[0]), reversed(included_spacers[1])):
					del dupe_indices[x][y]
				try:
					group = [min([x for y in dupe_indices for x in y])] # Restart the process with the next earliest spacer
				except: # If all the dupe_indices have been deleted this will fail and we'll be finished
					break
				included_spacers = [[],[]]


	return dupe_groups


def count_parsimony_events(child, ancestor):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
	
	Returns:
		(Array class instance) The child Array class instance with parsimony distance added to the .distance attribute.
	"""

	child, ancestor = find_modules(child, ancestor)


	# Process modules to count events since ancestor

	# First look for duplications

	dupes = find_dupes(child, ancestor)
	dupe_indices = [x for y in dupes for x in y] # Flatten list to check against later

	child.events['duplication'] = int(len(dupes) / 2)

	idx = 0
	while idx < max(child.module_lookup.keys()):
		if idx in dupe_indices:
			idx += 1
			continue
		mod = child.module_lookup[idx]
		if mod.type == 'acquisition' or mod.type == 'no_acquisition': 
		# Checking no acquisition allows comparison of child-child as well as child-ancestor arrays.
			child.events['acquisition'] += len(mod.indices)
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
			continue
		if mod.type == 'indel':
			child.events['indel'] += 1
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
		else:
			idx += 1

	return child


def resolve_pairwise_parsimony(array1, array2, all_arrays, node_ids, node_count):
	"""
	Given two arrays, make a hypothetical ancestral state and calculate parsimony distance of each input array to that ancestor. 
	Can only build an ancestor state if there are shared spacers so throws an error if none are found.
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
		node_ids (list): A list of names to be used to name internal nodes.
		node_count (int): The index of the name to use in the node_ids list to name this internal node.

	Returns:
		(tuple of Array class instances) The input Array class instances with module and distance info added, the ancestral array Array class instance. Tuple order is (array1, array2, ancestor)

	Raises:
		No ID: Raises an exception if the two arrays share no spacers.
	"""

	if len(list(set(array1.spacers) & set(array2.spacers))) > 0:

		array1, array2 = find_modules(array1, array2)

		array1.reset() # Make sure the distance and events haven't carried over from a previous usage
		array2.reset()

		ancestor = infer_ancestor(array1, array2, all_arrays, node_ids, node_count)

		array1 = count_parsimony_events(array1, ancestor)

		array2 = count_parsimony_events(array2, ancestor)
		

		for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
			array1.distance += array1.events[k] * v
			array2.distance += array2.events[k] * v

		return array1, array2, ancestor
	

	else:
		return "No_ID"


def find_closest_array(array, array_dict):
	"""
	Args:
		array (Array class instance): The array you want to find the closest match for.
		array_dict (dict): The dictionary with values of Array class instances of arrays already in your tree
	
	Returns:
		(Array class instance) The array already in your tree that is the most parsimonious match to the query array.
	"""

	best_score = 9999999999

	for array_id, comparator_array in array_dict.items():
		x = count_parsimony_events(comparator_array, array)
		for k,v in event_costs.items():
			x.distance += x.events[k] * v
		if x.distance < best_score:
			best_score = x.distance
			best_match = array_dict[x.id]

	return best_match


def replace_existing_array(existing_array, new_array, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed):
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
	Returns:
		(tuple) of the following:
			(dendropy Tree instance) Tree with the existing array replaced with the newly formed node.
			(dict) Dict of arrays in the tree so far.
			(dict) Dict of dendopy tree node class instances in the tree.
	"""

	# Make a hypothetical ancestor for the pair of arrays and calculate the distances
	results = resolve_pairwise_parsimony(existing_array, new_array, all_arrays, node_ids, node_count)
	if results == "No_ID":
		return results


	else:
		existing_array, new_array, ancestor = results
		tree_child_dict[existing_array.id].edge_length = existing_array.distance
		# If the current parent is not seed then do this
		if not seed:
			# Calculate distance from this hypothetical ancestor to its new parent in the tree
			ancestor.reset()
			ancestor = count_parsimony_events(ancestor, current_parent)
			for k,v in event_costs.items():
				ancestor.distance += ancestor.events[k] * v

		

		for a in [new_array, ancestor]:
			# Create tree nodes
			tree_child_dict[a.id] = dendropy.Node(edge_length=a.distance)
			tree_child_dict[a.id].taxon = taxon_namespace.get_taxon(a.id)
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


event_costs = { 
				"acquisition" : 1,
				"indel" : 1,
				"duplication": 1,
				}


# Generate strings to assign as internal node_IDs (This makes 702)

node_ids = ["Internal_" + i for i in ascii_lowercase]
if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
	node_ids += ["Internal_" + "".join(i) for i in product(ascii_lowercase, repeat = 2)]


array_spacers_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_spacers_dict[bits[0]] = bits[2:]

arrays = [Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
all_arrays = [array.spacers for array in arrays]
labels = args.arrays_to_join
 

taxon_namespace = dendropy.TaxonNamespace(args.arrays_to_join + node_ids)

best_score = 99999999

for i in range(100):

	addition_order = random.sample(arrays, len(arrays)) # Shuffle array order to build tree.
	# addition_order = [Array(i, array_spacers_dict[i]) for i in ['355', '433', '761', '1685', '1254', '159', '146', '487']]
	try:
		tree = dendropy.Tree(taxon_namespace=taxon_namespace)

		array_dict = {}
		tree_child_dict = {}
		node_count = 0 # Keep track of which internal node ID should be used for each node

		results = resolve_pairwise_parsimony(addition_order[0], addition_order[1], all_arrays, node_ids, node_count)
		if results != "No_ID":
			node_count += 1
			array1, array2, ancestor = results
		else:
			continue


		for a in [array1, array2, ancestor]:
			# Create tree nodes
			tree_child_dict[a.id] = dendropy.Node(edge_length=a.distance)
			tree_child_dict[a.id].taxon = taxon_namespace.get_taxon(a.id)
			#Store arrays for further comparisons
			array_dict[a.id] = a

		# Add initial relationships to the tree.
		tree.seed_node.add_child(tree_child_dict[ancestor.id])
		tree_child_dict[ancestor.id].add_child(tree_child_dict[array1.id])
		tree_child_dict[ancestor.id].add_child(tree_child_dict[array2.id])

		for a in addition_order[2:]: # Already added the first two so now add the rest 1 by 1
			seed = False # To check if we are modifying the child of the seed node

			# Find the most similar array already in the tree (measured in parsimony score)
			best_match = find_closest_array(a, array_dict)
			if best_match.extant: # If the closest match is a array then just join them and replace the existing array with the new node
				current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
				results = replace_existing_array(
					best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
					)
				if results == "No_ID":
					break
				else:
					tree, array_dict, tree_child_dict = results
				node_count += 1
			else:
				try:
					current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
				except: # Fails if the best match is a child of the seed node
					seed = True
					current_parent = None
				results = replace_existing_array(
					best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
					)
				if results == "No_ID":
					break
				else:
					tree, array_dict, tree_child_dict = results
				node_count += 1
			if a != addition_order[-1]:
				score = sum([v.distance for v in array_dict.values()])
				if score > best_score:
					break
		else:
			score = sum([v.distance for v in array_dict.values()])
			if score < best_score:
				best_score = copy.deepcopy(score)
				best_tree = copy.deepcopy(tree)
	except:
		print('Something went wrong when running with the following array order:')
		print([i.id for i in addition_order])


print('\n\n\n\n')
print(best_score)
print(best_tree.as_ascii_plot(show_internal_node_labels=True))
print(best_tree.as_string("newick"))#, suppress_leaf_node_labels=False, suppress_annotations=False))
