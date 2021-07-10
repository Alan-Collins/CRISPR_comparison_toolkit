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
		"""
	def __init__(self, ID, spacers=[], extant=True):
		self.id = ID
		self.extant = extant
		self.modules = []
		self.spacers = spacers
		self.aligned = []
		self.module_lookup = {}
		self.distance = 0

	def sort_modules(self):
		self.modules.sort(key=lambda x: int(x.indices[0]))

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
		


def needle(seq1, seq2, match = 20, mismatch = -1, gap = -2):
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
	node_count += 1

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
						elif spacers1_found:
							ancestor.modules.append(mod1)
						elif spacers2_found:
							ancestor.modules.append(mod2)
						else:
							# Neither was found.
							pass



	ancestor.spacers = [i.spacers for i in ancestor.modules]
	if isinstance(ancestor.spacers[0], list):
		ancestor.spacers = [spacer for sublist in ancestor.spacers for spacer in sublist]

	return ancestor


def count_parsimony_events(child, ancestor):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
	
	Returns:
		(Array class instance) The child Array class instance with parsimony distance added to the .distance attribute.
	"""

	child, ancestor = find_modules(child, ancestor)

	print(child.aligned)
	print(["{} : {}".format(i.type, i.spacers) for i in child.modules])


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

	array1, array2 = find_modules(array1, array2)

	ancestor = infer_ancestor(array1, array2, all_arrays, node_ids, node_count)

	array2 = count_parsimony_events(array2, ancestor)


# Generate strings to assign as internal node_IDs (This makes 702)

node_ids = ["Internal_" + i for i in ascii_lowercase]
if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
	node_ids += ["Internal_" + "".join(i) for i in product(ascii_lowercase, repeat = 2)]
node_count = 0 # Keep track of which internal node ID should be used for each node

array_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

arrays = [Array(i, array_dict[i]) for i in args.arrays_to_join]
labels = args.arrays_to_join


ancestor = resolve_pairwise_parsimony(arrays[0], arrays[1], [array.spacers for array in arrays], node_ids, node_count)

# print("{} : {}".format(ancestor.id, ancestor.spacers))

if args.print_tree:
	from ete3 import Tree
	#print(Tree(tree))
