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
from collections import defaultdict
import dendropy
import copy
from itertools import permutations
import matplotlib.pyplot as plt
from matplotlib import rcParams
import multiprocessing


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
	"-r",  dest="replicates", type=int, nargs="?", default = 1,
		help="Specify number of replicates of tree building to perform. The more replicates, the greater the chance that a better tree will be found. Default: 1"
	)
parser.add_argument(
	"-q",  dest="acquisition", type=int, nargs="?", default = 1,
		help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
	)
parser.add_argument(
	"-i",  dest="indel", type=int, nargs="?", default = 1,
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
	"-t",  dest="num_threads", type=int, nargs="?", default = 1,
		help="Specify number of threads to use for building trees. Using multiple threads will speed up the search for trees when performing many replicates with the -r option. Default: 1"
	)
parser.add_argument(
	"-y", dest="output_arrays", required = False,
	help="Specify filename for the details of you final arrays with hypothetical intermediate arrays in the same format as your input array_file. If there are multiple best trees, one file will be created per tree numbered in the order they are described in the stdout output."
	)
parser.add_argument(
	"arrays_to_join", nargs="*",  
	help="Specify the IDs of the arrays you want to join. If none provided, joins all arrays in the provided array representatives file. **If given, must come at the end of your command after all other arguments.**"
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
						"repeated_indel" : 0,
						"duplication": 0,
						"trailer_loss": 0
						}

	def sort_modules(self):
		self.modules.sort(key=lambda x: int(x.indices[0]))

	def reset(self):
		self.distance = 0
		self.events = { 
						"acquisition" : 0,
						"indel" : 0,
						"repeated_indel" : 0,
						"duplication": 0,
						"trailer_loss": 0
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
		partner (str): If the type of this module is "repeated_indel", then this stores the ID of the array that was identified as the other instance of this indel.
	"""
	def __init__(self):
		self.type = ""
		self.spacers = []
		self.indices = []
		self.partner = []
		

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
				# Check if duplication
				if a != b and Counter(array2.aligned)[b] > 1 and b != '-':
					# If this spacer is in multiple copies and aligns with a gap it means the comparator array has no or fewer copies of this spacer.
					if module2.type != "duplication" and module2.type != "":
						array2.modules.append(module2)
						for k in module2.indices:
							array2.module_lookup[k] = module2
						module2 = Spacer_Module()
					module2.type = "duplication"
					module2.indices.append(n)
					module2.spacers.append(b)

					if module1.type != "duplication" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = Spacer_Module()
					module1.type = "duplication"
					module1.indices.append(n)
					module1.spacers.append(a)
				elif a != b and Counter(array1.aligned)[a] > 1  and a != '-':
					if module1.type != "duplication" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = Spacer_Module()
					module1.type = "duplication"
					module1.indices.append(n)
					module1.spacers.append(a)

					if module2.type != "duplication" and module2.type != "":
						array2.modules.append(module2)
						for k in module2.indices:
							array2.module_lookup[k] = module2
						module2 = Spacer_Module()
					module2.type = "duplication"
					module2.indices.append(n)
					module2.spacers.append(b)

				else:
					if a != '-' and b != '-':
						# Module1 processing for indel module where mismatched
						if module1.type != "indel_mm" and module1.type != "":
							array1.modules.append(module1)
							for k in module1.indices:
								array1.module_lookup[k] = module1
							module1 = Spacer_Module()
						module1.type = "indel_mm"
						module1.indices.append(n)
						module1.spacers.append(a)
						# Module2 processing for indel module where mismatched
						if module2.type != "indel_mm" and module2.type != "":
							array2.modules.append(module2)
							for k in module2.indices:
								array2.module_lookup[k] = module2
							module2 = Spacer_Module()
						module2.type = "indel_mm"
						module2.indices.append(n)
						module2.spacers.append(b)
					else:
						if n == len(array1.aligned)-1: # If this is the last spacer, call it trailer loss
							array1.modules.append(module1)
							for k in module1.indices:
								array1.module_lookup[k] = module1
							module1 = Spacer_Module()
							module1.type = "trailer_loss"
							module1.indices.append(n)
							module1.spacers.append(a)
							array2.modules.append(module2)
							for k in module2.indices:
								array2.module_lookup[k] = module2
							module2 = Spacer_Module()
							module2.type = "trailer_loss"
							module2.indices.append(n)
							module2.spacers.append(b)
						else:
							# Module1 processing for indel module where one is gap
							if module1.type != "indel_gap" and module1.type != "":
								array1.modules.append(module1)
								for k in module1.indices:
									array1.module_lookup[k] = module1
								module1 = Spacer_Module()
							module1.type = "indel_gap"
							module1.indices.append(n)
							module1.spacers.append(a)
							# Module2 processing for indel module where one is gap
							if module2.type != "indel_gap" and module2.type != "":
								array2.modules.append(module2)
								for k in module2.indices:
									array2.module_lookup[k] = module2
								module2 = Spacer_Module()
							module2.type = "indel_gap"
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
				# If one array has just spacers and the other just gaps in the indel then pick the spacers as ancestral unless those spacers are singletons.
				if all([i == '-' for i in mod1.spacers]) or all([i == '-' for i in mod2.spacers]):
					# Figure out which one is all gaps and store the other one if it isn't singletons.
					if all([i == '-' for i in mod1.spacers]):
						# Check first if it's a duplication. If so, remove from ancestor.
						# If not a duplication check if these spacers are in another array. If so keep.
						duplication = False # Start with the assumption this isn't a duplication
						for sp in mod2.spacers:
							if sp != '-':
								count = 0
								for spacer in mod2.spacers:
									if sp == spacer:
										count += 1
									if count == 2: # If the same spacer is present twice then duplication has occured in this region.
										duplication = True
										break
						if duplication == False:
							singleton = [True for _ in mod2.spacers] # Start with the assumption that these spacers aren't found in another array.
							for n, sp in enumerate(mod2.spacers):
								count = 0
								for array in all_arrays:
									if sp in array:
										count += 1
									if count == 1: # If the spacer is found in another array than the two being considered here then it's not a singleton for our purposes.
										singleton[n] = False
										continue

							if all([_ == False for _ in singleton]):
								ancestor.modules.append(mod2)

							elif not singleton[0] and not singleton[-1]: 
								# If the outtermost are shared with other arrays, perhaps this is an insertion in one and a deletion in the other. That may be more likely than independent acquisition of two stretches of spacers twice.
								ancestor.modules.append(mod2)
							elif not all(singleton):
								# If spacers found in another array are on the edge of this indel it may be more parsimonious to have them in the ancestor than not.
								if not singleton[0]: # Start looking from the leader end
									spacers = []
									for n, tf in enumerate(singleton):
										if not tf:
											spacers.append(mod2.spacers[n])
										else:
											x = Spacer_Module()
											x.spacers = spacers
											ancestor.modules.append(x)
											break
								elif not singleton[-1]: # Start looking from the trailer end
									spacers = []
									for n, tf in enumerate(reversed(singleton)):
										if not tf:
											spacers.append(mod2.spacers[::-1][n])
										else:
											
											spacers.reverse() # We were working from trailer to leader so need to reverse now.
											x = Spacer_Module()
											x.spacers = spacers
											ancestor.modules.append(x)
											break

					else:
						duplication = False # Start with the assumption this isn't a duplication
						for sp in mod2.spacers:
							if sp != '-':
								count = 0
								for spacer in mod2.spacers:
									if sp == spacer:
										count += 1
									if count == 1: # If the same spacer is present twice then duplication has occured in this region.
										duplication = True
										break
						if duplication == False:
							singleton = [True for _ in mod1.spacers] # Start with the assumption that these spacers aren't found in another array.
							for n, sp in enumerate(mod1.spacers):
								count = 0
								for array in all_arrays:
									if sp in array:
										count += 1
									if count == 1:
										singleton[n] = False
										continue
							if all([_ == False for _ in singleton]):
								ancestor.modules.append(mod1)
							elif not singleton[0] and not singleton[-1]: 
								# If the outtermost are shared with other arrays, perhaps this is an insertion in one and a deletion in the other. That may be more likely than independent acquisition of two stretches of spacers twice.
								ancestor.modules.append(mod1)
							elif not all(singleton):
								# If spacers found in another array are on the edge of this indel it may be more parsimonious to have them in the ancestor than not.
								if not singleton[0]: # Start looking from the leader end
									spacers = []
									for n, tf in enumerate(singleton):
										if not tf:
											spacers.append(mod1.spacers[n])
										else:
											x = Spacer_Module()
											x.spacers = spacers
											ancestor.modules.append(x)
											break
								elif not singleton[-1]: # Start looking from the trailer end
									spacers = []
									for n, tf in enumerate(reversed(singleton)):
										if not tf:
											spacers.append(mod1.spacers[::-1][n])
										else:
											spacers.reverse() # We were working from trailer to leader so need to reverse now.
											x = Spacer_Module()
											x.spacers = spacers
											ancestor.modules.append(x)
											break
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
							idx = mod1.indices[-1] + 1
						elif " ".join(mod2.spacers + mod1.spacers) in " ".join(found_array):
							new_mod = Spacer_Module()
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


def count_parsimony_events(child, ancestor, array_dict, tree, parent_comparison):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
		array_dict (dict): Dict of Array class instances with information about the nodes of your tree.
		tree (Deondropy Tree class instance): The tree in which the arrays are located.	
		parent_comparison (bool): Are the arrays being compared parent-child nodes in the tree? E.g. if you are trying to find the best matching array in the tree this is False. If you are counting events between an array and its ancestor array this is True

	Returns:
		(Array class instance) The child Array class instance with parsimony events added to the .events dict attribute.
	"""

	child, ancestor = find_modules(child, ancestor)

	# Process modules to count events since ancestor

	# First look for duplications

	dupes = find_dupes(child, ancestor)
	dupe_indices = [x for y in dupes for x in y] # Flatten list to check against later

	child.events['duplication'] = int(len(dupes) / 2)

	idx = 0
	while idx <= max(child.module_lookup.keys()):
		if idx in dupe_indices:
			idx += 1
			continue
		mod = child.module_lookup[idx]
		if mod.type == 'acquisition' or mod.type == 'no_acquisition': 
		# Checking no acquisition allows comparison of child-child as well as child-ancestor arrays.
			if parent_comparison: 
				# If comparing to a nodes ancestor, then deletions of spacers from the leader end may look like acquisition. This must be checked to make sure deletions are not being miscalled as acquisitions.
				child = identify_repeat_indels(child, ancestor, array_dict, mod, ancestor.module_lookup[idx], tree)
			else:
				child.events['acquisition'] += len(mod.indices)
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
			continue
		if mod.type == 'indel_gap' or mod.type == "indel_mm":
			if parent_comparison:
				child = identify_repeat_indels(child, ancestor, array_dict, mod, ancestor.module_lookup[idx], tree)
			else:
				child.events['indel'] += 1
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
		if mod.type == 'duplication':
			child.events['duplication'] += 1
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
		if mod.type == "trailer_loss":
			child.events['trailer_loss'] += 1
			idx = mod.indices[-1] + 1
			pass
		else:
			idx += 1

	return child


def identify_repeat_indels(child, ancestor, array_dict, module, ancestor_module, tree):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
		array_dict (dict): Dict of Array class instances with information about the nodes of your tree.
		module (Spacer_Module class instance): The indel module to be processed.
		ancestor_module (Spacer_Module class instance): The module of the parent of the array to be processed.
		tree (Deondropy Tree class instance): The tree in which the arrays are located.
	
	Returns:
		(Array class instance) The child array, modified to contain newly identified Spacer_Module instances if any are found.
	"""

	indels = repeated_indels = acquisition = 0 # Start counters to keep track of events identified
	new_modules = [] # Initialize a list to contain new modules if any are found
	if len(tree) > 1:
		ancestor_node = tree.find_node_with_taxon_label(ancestor.id)
		if ancestor_node:
			ancestor_children = ancestor_node.leaf_nodes() + [node for node in ancestor_node.levelorder_iter(lambda x: x.is_internal())]
			ancestor_children_ids = [node.taxon.label for node in ancestor_children]

			non_ancestor_children = []
			for node in tree:
				if node.taxon:
					if node.taxon.label not in ancestor_children_ids:
						non_ancestor_children.append(node)

			non_ancestor_children_ids = [node.taxon.label for node in non_ancestor_children]

			non_ancestor_children_spacers = [array_dict[array].spacers for array in non_ancestor_children_ids]

			child_node = tree.find_node_with_taxon_label(child.id)
			other_child_children = [] # Store a list of the children that are not the child array being assessed. If the child is already set as the child of this ancestor then this will just be that child's sibling. Otherwise this is a comparison used to place the "child" array which is not yet in the tree and so all the children of the ancestor should be included.
			if child_node:
				if child_node.parent_node == ancestor_node:
					other_child = child_node.sibling_nodes()[0]
					other_child_children = other_child.leaf_nodes() + [node for node in other_child.levelorder_iter(lambda x: x.is_internal())]
				else:
					for n in ancestor_node.child_nodes():
						other_child_children += n.leaf_nodes() + [node for node in n.levelorder_iter(lambda x: x.is_internal())]
			else:
				for n in ancestor_node.child_nodes():
					other_child_children += n.leaf_nodes() + [node for node in n.levelorder_iter(lambda x: x.is_internal())]
			other_child_children_ids = [node.taxon.label for node in other_child_children]
			other_child_children_spacers = [array_dict[array].spacers for array in other_child_children_ids]



			spacer_hits = defaultdict(list)
			repeated_indel_instances = []
			repeat_search = True
			original_spacers_to_check = spacers_to_check = module.spacers # Store which spacers to look for and an original copy to maintain index information
			array_IDs_of_concern = non_ancestor_children_ids + other_child_children_ids
			arrays_of_concern = non_ancestor_children_spacers + other_child_children_spacers
			while repeat_search:
				longest_match = 0
				longest_indices = []
				arrays_to_check = []
				# If lost spacers are found in either of these two groups it indicates independent gain of the same spacers in different lineages.
				for spacer in spacers_to_check:
					for a, array in enumerate(arrays_of_concern):
						if spacer in array:
							arrays_to_check.append(a)		
				if len(arrays_to_check) == 0:
					repeat_search = False 
				else:
					for a in arrays_to_check:
						x = Array("x", spacers_to_check)
						y = Array("y", arrays_of_concern[a])
						x.aligned, y.aligned = needle(x.spacers, y.spacers) # Find where the match is by aligning
						x,y = find_modules(x,y) # Then identify any shared regions
						for m in x.modules: # Pull out the shared modules
							if m.type == "shared":
								if len(m.indices) > longest_match:
									longest_match = len(m.indices)
									longest_indices = [original_spacers_to_check.index(s) for s in m.spacers if s in original_spacers_to_check]
									longest_spacers = [s for s in spacers_to_check if s in m.spacers]
									partner = [array_IDs_of_concern[a]]
									partner_extant = array_dict[array_IDs_of_concern[a]].extant
								elif len(m.indices) == longest_match and (m.spacers) == set(longest_spacers):
									# If the modules are equally long, prefer to keep extant arrays.
									if not partner_extant and array_dict[array_IDs_of_concern[a]].extant:
										longest_match = len(m.indices)
										longest_indices = [original_spacers_to_check.index(s) for s in m.spacers if s in original_spacers_to_check]
										longest_spacers = [s for s in spacers_to_check if s in m.spacers]
										partner = [array_IDs_of_concern[a]]
										partner_extant = array_dict[array_IDs_of_concern[a]].extant
									else:
										if array_dict[array_IDs_of_concern[a]].extant:
											if array_IDs_of_concern[a] not in partner:
												partner.append(array_IDs_of_concern[a])
					new_rep_indel_mod = Spacer_Module()
					new_rep_indel_mod.spacers = longest_spacers
					new_rep_indel_mod.indices = [module.indices[i] for i in longest_indices]
					new_rep_indel_mod.type = "repeated_indel"
					new_rep_indel_mod.partner = partner
					new_modules.append(new_rep_indel_mod)
					repeated_indels += 1
					# Remove spacers that have been found from list and look again
					spacers_to_check = [s for s in spacers_to_check if s not in longest_spacers]
			if len(spacers_to_check) > 0:
				# identify consecutive runs of spacers that are unique to this array
				spacer_indices = [n for n, s in enumerate(module.spacers) if s in spacers_to_check]
				if module.type == "acquisition":
					expected_index = 0 # If this is an acquisition then the first spacer will still be present.
					new_ac_mod = Spacer_Module()
					new_ac_mod.type = "acquisition"
					while spacer_indices[0] == expected_index:
						new_ac_mod.indices.append(module.indices[spacer_indices[0]])
						new_ac_mod.spacers.append(module.spacers[spacer_indices[0]])
						acquisition += 1 # If true then this spacer was acquired since ancestor.
						expected_index = spacer_indices.pop(0)+1 # Update expected index to be the one higher. Remove the index from the list so it isn't later counted as an indel
						if len(spacer_indices) == 0: # If we finish the list then break the while loop
							break
					if len(new_ac_mod.spacers) > 0: # If any acquisition was found then add the new acquisition module to the list of new modules
						new_modules.append(new_ac_mod)
				if spacers_to_check != module.indices: # If we haven't found anything, no need to update the existing indel. This check is passed if we found something
					if len(spacer_indices) > 0: # Then find consecutive runs of spacers. Each must be an indel
						new_indel_mod = Spacer_Module() # Make a new Spacer_Module to store the indel.
						new_indel_mod.type = "indel" # Call it just an idel as gap or mm relative to comparator is not assessed here.
						last_n = False # Keep track of what the last index was to know if this is a consecutive run.
						for n in spacer_indices:
							if n == spacer_indices[-1]: # If the list is over then add the last indel and break
								new_indel_mod.indices.append(module.indices[n])
								new_indel_mod.spacers.append(module.spacers[n])
								new_modules.append(new_indel_mod)
								indels += 1
								break
							if last_n:
								if n != last_n + 1: # If this number is not one more than the last then consecutive run is over
									new_modules.append(new_indel_mod)
									new_indel_mod = Spacer_Module() # Make a new Spacer_Module to store the next indel.
									new_indel_mod.type = "indel"
									indels += 1
								last_n = n
								new_indel_mod.indices.append(module.indices[n])
								new_indel_mod.spacers.append(module.spacers[n])
							else: # Otherwise move to the next number
								last_n = n 
								new_indel_mod.indices.append(module.indices[n])
								new_indel_mod.spacers.append(module.spacers[n])
			if module.type == 'no_acquisition':
				# None of the above will have run if this module is just gaps.
				# In that case we need to check if the absence of those spacers in the child array represents simply not having acquired them or if it would require the loss and regain of the same spacers at different points in the tree.
				# If the lost spacers are present in the children of this array then there has to have been a deletion and regain. That probably means this arrangement is wrong so I will assign a high parsimony cost in order to favour other topologies.
				ancestor_children_spacers = [array_dict[array].spacers for array in ancestor_children_ids]
				# Do a simple check to see if all the lost spacers are in any of the children
				found = False
				for array in ancestor_children_spacers:
					if all([spacer in array for spacer in ancestor_module.spacers]):
						found = True
						break
				if found: # If we found all the lost spacers in a child array then this is at least one repeated deletion. Add one event for now. May implement more detailed characterization in future if trees don't look right.
					repeated_indels += 1
				else: # If not then just add the cost of the acquisitions that must have occured between the child and ancestor (even though the direction of that change is opposite to normal comparisons.)
					acquisition += len(module.indices)
		else:
			indels = 1
	else:
		if module.type == "acquisition":
			acquisition += len(module.indices)
		else:
			indels = 1

	child.events['acquisition'] += acquisition
	child.events['indel'] += indels
	child.events['repeated_indel'] += repeated_indels

	if len(new_modules) > 0: # If modules have been found within the query module then it should be replaced in the .modules list with the newly found modules and the .module_lookup dict updated.
		new_modules.sort(key=lambda x: int(x.indices[0])) # Sort the modules in order of indices
		idx = child.modules.index(module) # find the index of the module being replaced in the .modules list
		del child.modules[idx] # Remove the old module from the list.
		for new_mod in new_modules:
			# Add the new module to the .modules and .lookup_modules
			child.modules.append(new_mod)
			for idx in new_mod.indices:
				child.module_lookup[idx] = new_mod
		child.sort_modules() # Sort the modules so the newly added ones are in the correct order.

			
	return child


def resolve_pairwise_parsimony(array1, array2, all_arrays, array_dict, node_ids, node_count, tree):
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

		array1 = count_parsimony_events(array1, ancestor, array_dict, tree, True)
		array2 = count_parsimony_events(array2, ancestor, array_dict, tree, True)

		for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
			array1.distance += array1.events[k] * v
			array2.distance += array2.events[k] * v

		return array1, array2, ancestor
	

	else:
		return "No_ID"


def find_closest_array(array, array_dict, tree):
	"""
	Args:
		array (Array class instance): The array you want to find the closest match for.
		array_dict (dict): The dictionary with values of Array class instances of arrays already in your tree
		tree (Deondropy Tree class instance): The tree in which the arrays are located.

	Returns:
		(Array class instance) The array already in your tree that is the most parsimonious match to the query array.
	"""

	best_score = 9999999999
	best_match = False
	for array_id in array_dict.keys():
		comparator_array = copy.deepcopy(array_dict[array_id])
		comparator_array.reset()
		comparator_array = count_parsimony_events(comparator_array, array, array_dict, tree, False)
		array.reset()
		array = count_parsimony_events(array, comparator_array, array_dict, tree, True)
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
	results = resolve_pairwise_parsimony(existing_array, new_array, all_arrays, array_dict, node_ids, node_count, tree)
	
	if results == "No_ID":
		return results

	else:
		existing_array, new_array, ancestor = results
		tree_child_dict[existing_array.id].edge_length = existing_array.distance
		# If the current parent is not seed then do this
		if not seed:
			# Calculate distance from this hypothetical ancestor to its new parent in the tree
			ancestor.reset()
			ancestor = count_parsimony_events(ancestor, current_parent, array_dict, tree, True)
			for k,v in event_costs.items():
				ancestor.distance += ancestor.events[k] * v

		# modify the edge length for the existing array and update its distance in the array_dict
		tree_child_dict[existing_array.id].edge_length = existing_array.distance
		array_dict[existing_array.id] = existing_array
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


def plot_tree(tree, array_dict, filename):
	"""
	Args:
		tree (dendopy Tree class instance): The tree you want to plot.
		array_dict (dict): Dict of Array class instances with information about the nodes of your tree.
		filename (str): Path to the file you want created with this plot.
	"""

	# Find tree dimensions
	tree_width = max(tree.calc_node_root_distances(return_leaf_distances_only=True))
	tree_height = len(tree.nodes())

	#Find deepest node so plotting can start at the most distant leaf
	for node in tree.postorder_node_iter():
		if node.root_distance == tree_width:
			start_node = node

	node_locs = {} # Store where each node is located to draw lines to it.

	start_position = [0.5,0.5]
	node_locs[start_node.taxon.label] = start_position
	dim_x, dim_y = 10, 10

	hscale = (dim_x+1)/(tree_width + max([len(array.spacers) for array in array_dict.values()])) # Factor to scale all branch lengths to fit them in the plot
	vscale = (dim_y+1)/tree_height

	max_depth = 0.5+tree_width*hscale

	fig, ax = plt.subplots()

	fig.set_size_inches(dim_x, dim_y)

	node = start_node
	highest_y = 0.5

	nodes_to_revisit = {} # Store nodes with subtrees and the y value to start at for them

	while True: # Keep going until reaching the root triggers a break.

		first_node = node
		if not first_node.taxon.label in node_locs.keys():
			node_locs[first_node.taxon.label] = position
		if len(node.sibling_nodes())==1:
			if node.taxon.label not in nodes_to_revisit.keys():
				second_node = node.sibling_nodes()[0]
			else:
				del nodes_to_revisit[node.taxon.label] # Remove that from the list as we've finished its tree
				if len(nodes_to_revisit) == 0:
					break
				else:
					node = tree.find_node_with_taxon_label(list(nodes_to_revisit.keys())[0]) # start the next node to revisit
					first_node = node.leaf_nodes()[0]# start drawing from a random leaf in the subtree
					second_node = first_node.sibling_nodes()[0]
					highest_y = y2 = nodes_to_revisit[node.taxon.label] # Set the y position reserved for this subtree
					second_node = first_node.sibling_nodes()[0]
					node_locs[first_node.taxon.label] = [max_depth-first_node.root_distance*hscale,y2]
				
				
				


		else: # We've reached the root. Check if any subtrees need to be figured out.
			if len(nodes_to_revisit) == 0: # No subtrees. Exit the loop
				break
			else:
				node = tree.find_node_with_taxon_label(list(nodes_to_revisit.keys())[0]) # start the first node to revisit
				first_node = node.child_nodes()[0].leaf_nodes()[0] # start drawing from a random leaf in the subtree
				second_node = first_node.sibling_nodes()[0]
				highest_y = y2 = nodes_to_revisit[node.taxon.label] # Set the y position reserved for this subtree
				node_locs[first_node.taxon.label] = [max_depth-first_node.root_distance*hscale,y2]
		

		# figure out first branch location
		
		x1 = node_locs[first_node.taxon.label][0]
		x2 = node_locs[first_node.taxon.label][0] + first_node.edge_length * hscale
		y1 = node_locs[first_node.taxon.label][1]
		y2 = node_locs[first_node.taxon.label][1]

		highest_y = max([highest_y,y2])
		
		num_leaves = len(second_node.leaf_nodes()) # Figure out how much space is needed based on the number of leaves below this node
		num_internal = len([i for i in second_node.levelorder_iter(lambda x: x.is_internal())])

		if num_internal > 0:  
			num_internal += num_leaves # - 1 # Counts self so need to subtract 1.

			node_locs[second_node.parent_node.taxon.label] = [
			max_depth-second_node.parent_node.root_distance*hscale,
			highest_y+1.5*vscale
			]

			y2 = highest_y+num_internal*1.5*vscale

			# Leave space for subtree
			highest_y = highest_y+(num_internal+1)*1.5*vscale


			position = [max_depth-second_node.root_distance*hscale ,y2]
			if second_node.taxon.label not in node_locs.keys():
				node_locs[second_node.taxon.label] = position

			# figure out second branch location

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length * hscale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			nodes_to_revisit[second_node.taxon.label] = y2-((num_internal-1)*1.5)*vscale+1.5*vscale # store name of subtree parent and position to start drawing subtree


		else:
			y2 = highest_y+3*vscale

			highest_y = y2

			position = [max_depth-second_node.root_distance*hscale ,y2]
			if second_node.taxon.label not in node_locs.keys():
				node_locs[second_node.taxon.label] = position

			# figure out second branch location

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length * hscale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]


			y1 = node_locs[first_node.taxon.label][1]
			position = [max_depth-second_node.parent_node.root_distance*hscale ,y2-1.5*vscale]
		node = second_node.parent_node

	# Add cartoon arrays to show hypothetical ancestral states
	# plot each array using the coordinates of the array label on the plotted tree.
	rep_indel_report_count = 1 # If repeat_indel found with multiple arrays as possible partners, annotate one on the tree and print the rest to stdout if user wants emphasis.
	rep_indel_message_printed = False # Print a help message the first time a list of possible rep indel partners are found.

	for array, location in node_locs.items():
		# Add label first
		x ,y = location
		ax.text(x-0.25*hscale, y-0.4*hscale, array, ha='right', fontsize=50*hscale)
		# then add branches
		first_node = tree.find_node_with_taxon_label(array)

		# First add branch lengths if user desires

		if args.branch_lengths:
			if first_node.edge_length != 0:
				ax.text(x+(first_node.edge_length/2)*hscale, y-1.5*hscale, first_node.edge_length, ha='center', fontsize=30*hscale)
		
		# Draw first branch
		
		x1 = node_locs[first_node.taxon.label][0]
		x2 = node_locs[first_node.taxon.label][0] + first_node.edge_length * hscale
		y1 = node_locs[first_node.taxon.label][1]
		y2 = node_locs[first_node.taxon.label][1]

		ax.plot([x1, x2], [y1, y2], color='black', linewidth = 1*vscale, solid_capstyle="butt")

		if len(first_node.sibling_nodes()) == 1:
			second_node = first_node.sibling_nodes()[0]
			# draw second branch

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length * hscale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			ax.plot([x1, x2], [y1, y2], color='black', linewidth = 1*vscale, solid_capstyle="butt")

			# draw line between branches

			x1 = x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length * hscale
			y1 = node_locs[first_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			ax.plot([x2, x2], [y1, y2], color='black', linewidth = 1*vscale, solid_capstyle="butt")

			# Then add spacers and highligh differences

			ancestor = 	array_dict[first_node.parent_node.taxon.label]
			child = array_dict[array]
			child = count_parsimony_events(child, ancestor, array_dict, tree, True)

			if args.emphasize_diffs:
				start_pos_x = location[0]-5*hscale # Start a bit to the left to leave room for the label
				start_pos_y = location[1]
				spacer_count = 0 # How many spacers have been plotted?
				reshift_loc = 1000
				for n, diff_type in reversed(child.module_lookup.items()):
					spacer = child.aligned[n]

					# Add change info
					if n == reshift_loc:
						start_pos_x-=0.5*hscale # Shift future spacers to make space for line

					if n == diff_type.indices[-1]:
						if diff_type.type != 'shared':
							if diff_type.type == "repeated_indel":
								# Plot a red box around repeated indels
								nspacers = len([child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])
								# First bar
								ax.fill_between([start_pos_x-2*spacer_count*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale],start_pos_y+0.5*vscale, start_pos_y-0.5*vscale, color="#cc3300", edgecolor='none')
								# Second bar
								ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale], start_pos_y+0.5*vscale, start_pos_y-0.5*vscale, color="#cc3300", edgecolor='none')
								# Top bar
								ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale], start_pos_y+0.3*vscale, start_pos_y+0.5*vscale, color="#cc3300", edgecolor='none')
								# Bottom bar
								ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale], start_pos_y-0.3*vscale, start_pos_y-0.5*vscale, color="#cc3300", edgecolor='none')
								
								# Add Array ID of the array in which the spacers of this predicted repeated_indel can be found
								if len(diff_type.partner) > 2:
									ax.text(start_pos_x-2*(spacer_count+nspacers/2)*hscale-0.5*hscale, start_pos_y-1*vscale, "\n".join([diff_type.partner[0], "event {}".format(rep_indel_report_count)]), color="#cc3300", ha='center', fontsize=40*hscale)
									if not rep_indel_message_printed:
										print("Repeated indels were identified with multiple possible partners. Those cases will be annotated in the tree png file with the one of the arrays identified as a partner followed by an event number corresponding to one of the lists of partner arrays below:\n\n")
										rep_indel_message_printed = True
									print("Event {}: {}\n\n".format(rep_indel_report_count, " ".join(diff_type.partner)))
									rep_indel_report_count += 1

								else:
									if len(diff_type.partner) == 2:
										ax.text(start_pos_x-2*(spacer_count+nspacers/2)*hscale-0.5*hscale, start_pos_y-1*vscale, "\n".join(diff_type.partner), color="#cc3300", ha='center', fontsize=40*hscale)
									else:
										ax.text(start_pos_x-2*(spacer_count+nspacers/2)*hscale-0.5*hscale, start_pos_y-0.8*vscale, diff_type.partner[0], color="#cc3300", ha='center', fontsize=40*hscale)

								start_pos_x-=0.5*hscale # Shift future spacers a bit to make spacer for this line.
								# Shift again after the indel region
								reshift_loc = diff_type.indices[0]-1
								
							if diff_type.type == 'indel_gap' or diff_type.type == 'indel_mm' or diff_type.type == 'indel': 
								if spacer == '-':
									ax.plot([start_pos_x-2*spacer_count*hscale, start_pos_x-2*spacer_count*hscale-2*hscale],[start_pos_y+0.3*vscale, start_pos_y-0.3*vscale], color="#666666", linewidth=3*vscale, solid_capstyle="butt")
									ax.plot([start_pos_x-2*spacer_count*hscale, start_pos_x-2*spacer_count*hscale-2*hscale],[start_pos_y-0.3*vscale, start_pos_y+0.3*vscale], color="#666666", linewidth=3*vscale, solid_capstyle="butt")
									spacer_count+=1 # Shift future spacers a bit to make spacer for this line.
								else:
									nspacers = len([child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])
									# First bar
									ax.fill_between([start_pos_x-2*spacer_count*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale],start_pos_y+0.5*vscale, start_pos_y-0.5*vscale, color="#666666", edgecolor='none')
									# Second bar
									ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale], start_pos_y+0.5*vscale, start_pos_y-0.5*vscale, color="#666666", edgecolor='none')
									# Top bar
									ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale], start_pos_y+0.3*vscale, start_pos_y+0.5*vscale, color="#666666", edgecolor='none')
									# Bottom bar
									ax.fill_between([start_pos_x-2*(spacer_count+nspacers)*hscale-0.5*hscale-0.5*hscale, start_pos_x-2*spacer_count*hscale], start_pos_y-0.3*vscale, start_pos_y-0.5*vscale, color="#666666", edgecolor='none')

									start_pos_x-=0.5*hscale # Shift future spacers a bit to make spacer for this line.
									# Shift again after the indel region
									reshift_loc = diff_type.indices[0]-1

							elif diff_type.type == "acquisition":
								nspacers = len(diff_type.indices)

								rcParams['path.sketch'] = (50*vscale, 100*vscale, 1)
								ax.plot(np.linspace(start_pos_x-2*(spacer_count+nspacers)*hscale,start_pos_x-2*spacer_count*hscale,3),[start_pos_y-0.5*vscale]*3,color="#666666", linewidth=3*vscale, solid_capstyle="butt")
								rcParams['path.sketch'] = (0, 0, 0)

							elif diff_type.type == "trailer_loss":
								if spacer == '-': # Draw a single sloped line
									ax.plot([start_pos_x-2*spacer_count*hscale, start_pos_x-2*spacer_count*hscale-2*hscale],[start_pos_y+0.3*vscale, start_pos_y-0.3*vscale], color="#666666", linewidth=3*vscale, solid_capstyle="butt")
									spacer_count+=1 # Shift future spacers a bit to make spacer for this line.


					# Plot spacer cartoon
					if spacer != '-':
						if spacer in spacer_cols_dict.keys():
							spcolour = spacer_cols_dict[spacer]
							line_width = 0.1
						else:
							spcolour = ("#000000", "#000000") #black
							line_width = 0.04
						# ax.plot([start_pos_x-2*spacer_count*hscale, start_pos_x-2*spacer_count*hscale-2*hscale],[start_pos_y, start_pos_y], color=spcolour, linewidth=line_width, solid_capstyle="butt")
						ax.fill_between([start_pos_x-(2*spacer_count+0.03)*hscale, start_pos_x-(2*spacer_count+0.03)*hscale-1.97*hscale], start_pos_y-line_width*vscale, start_pos_y+line_width*vscale, color=spcolour[0], edgecolor=spcolour[1], linewidth=4*hscale)
						spacer_count+=1




		else: # Draw root ancestral array
			if args.emphasize_diffs:
			# Then add spacers
				spacers = array_dict[array].spacers
				start_pos_x = location[0]-5*hscale # Start a bit to the left to leave room for the label
				start_pos_y = location[1] 
				for n, spacer in enumerate(reversed(spacers)): # work backwards through the array plotting from right to left
					if spacer in spacer_cols_dict.keys():
						spcolour = spacer_cols_dict[spacer]
						line_width = 0.1
					else:
						spcolour = ("#000000", "#000000") #black
						line_width = 0.04
					# ax.plot([start_pos_x-2*n*hscale, start_pos_x-2*n*hscale-2*hscale],[start_pos_y, start_pos_y], color=spcolour, linewidth=line_width, solid_capstyle="butt")
					ax.fill_between([start_pos_x-(2*n+0.03)*hscale, start_pos_x-(2*n+0.03)*hscale-1.97*hscale], start_pos_y-line_width*vscale, start_pos_y+line_width*vscale, color=spcolour[0], edgecolor=spcolour[1], linewidth=4*hscale)


		if not args.emphasize_diffs:
			# Then add spacers
			spacers = array_dict[array].spacers
			start_pos_x = location[0]-5*hscale # Start a bit to the left to leave room for the label
			start_pos_y = location[1] 
			for n, spacer in enumerate(reversed(spacers)): # work backwards through the array plotting from right to left
				if spacer in spacer_cols_dict.keys():
					spcolour = spacer_cols_dict[spacer]
					line_width = 0.1 # 10*vscale
				else:
					spcolour = ("#000000", "#000000") #black
					line_width = 0.04 # 4*vscale
				# ax.plot([start_pos_x-2*n*hscale, start_pos_x-2*n*hscale-2*hscale],[start_pos_y, start_pos_y], color=spcolour, linewidth=line_width, solid_capstyle="butt")
				ax.fill_between([start_pos_x-(2*n+0.03)*hscale, start_pos_x-(2*n+0.03)*hscale-1.97*hscale], start_pos_y-line_width*vscale, start_pos_y+line_width*vscale, color=spcolour[0], edgecolor=spcolour[1], linewidth=4*hscale)

	plt.axis('off')
	plt.tight_layout()
	plt.savefig(filename, dpi=600)


def build_tree_single(arrays, tree_namespace, score, all_arrays):
	"""
	Search treespace for most parsimonious tree using single process.
	Args:
		arrays (list): Ordered list of Array class instances of the arrays to analyse. Will be added to the tree in the provided order.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
		score (int): The score to beat. If at any point during tree construction the tree has total branch lengths above this score construction will be aborted
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
	
	Returns:
		(tuple) Returns array_dict and dendropy.tree object if the tree beats or equals the provide score. Else retuns tuple of (False, False).
	"""

	tree = dendropy.Tree(taxon_namespace=taxon_namespace)

	array_dict = {}
	tree_child_dict = {}
	node_count = 0 # Keep track of which internal node ID should be used for each node
	# Remove the arrays being compared from checks to see what spacers are in other arrays.
	# That way we only worry about ancestral states accounting for deletion of spacers in arrays not yet added to the tree.
	all_arrays = [a for a in all_arrays if a not in [arrays[0].spacers, arrays[1].spacers]]
	results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree)

	while results == "No_ID":
		arrays.append(arrays[1])
		del arrays[1]
		results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays,  array_dict, node_ids, node_count, tree)
	node_count += 1
	array1, array2, ancestor = results
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

	if len(arrays) != 2:
		for i in range(2, len(arrays)): # Already added the first two so now add the rest 1 by 1
			a = arrays[i]
			all_arrays = [arr for arr in all_arrays if arr != a.spacers]
			seed = False # To check if we are modifying the child of the seed node
			# Find the most similar array already in the tree (measured in parsimony score)
			best_match = find_closest_array(a, array_dict, tree)
			while best_match == "No_ID":

				arrays.append(arrays[i])
				del arrays[i]
				a = arrays[i]
				best_match = find_closest_array(a, array_dict, tree)
			if best_match.extant: # If the closest match is a array then just join them and replace the existing array with the new node
				current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
				results = replace_existing_array(
					best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
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
					best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
					)

				tree, array_dict, tree_child_dict = results
				node_count += 1

			# Recheck child - ancestor branch length to find indels that would have to occur multiple times
			count = 0
			for node in tree:
				if node.level() != 0:
					if node.parent_node.level() != 0:
						count+=1
						node_array = array_dict[node.taxon.label]
						parent_array = array_dict[node.parent_node.taxon.label]
						node_array.reset()
						parent_array.reset()
						node_array = count_parsimony_events(node_array, parent_array, array_dict, tree, True)
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							node_array.distance += node_array.events[k] * v
						node.edge_length = node_array.distance

		
			brlen = tree.length()
			if brlen > score:
				return (False, False)

					
			if a == arrays[-1]:
				tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works
				return (array_dict, tree)


	else:
		return (array_dict, tree)


def build_tree_multi(arrays, tree_namespace, all_arrays):
	"""
	Search treespace for most parsimonious tree using multiple processes.
	Args:
		arrays (list): Ordered list of Array class instances of the arrays to analyse. Will be added to the tree in the provided order.
		tree_namespace (dendropy.TaxonNamespace): Namespace for taxa to add to the tree.
		all_arrays (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
	
	Returns:
		(tuple) Returns array_dict and dendropy.tree object.
	"""

	# Remove the arrays being compared from checks to see what spacers are in other arrays.
	# That way we only worry about ancestral states accounting for deletion of spacers in arrays not yet added to the tree.
	all_arrays = [a for a in all_arrays if a not in [arrays[0].spacers, arrays[1].spacers]]

	tree = dendropy.Tree(taxon_namespace=taxon_namespace)

	array_dict = {}
	tree_child_dict = {}
	node_count = 0 # Keep track of which internal node ID should be used for each node

	results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree)
	while results == "No_ID":

		arrays.append(arrays[1])
		del arrays[1]
		results = resolve_pairwise_parsimony(arrays[0], arrays[1], all_arrays, array_dict, node_ids, node_count, tree)
	node_count += 1
	array1, array2, ancestor = results

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
	for i in range(2, len(arrays)): # Already added the first two so now add the rest 1 by 1
		a = arrays[i]
		all_arrays = [arr for arr in all_arrays if arr != a.spacers]
		seed = False # To check if we are modifying the child of the seed node

		# Find the most similar array already in the tree (measured in parsimony score)
		best_match = find_closest_array(a, array_dict, tree)
		while best_match == "No_ID":
			arrays.append(arrays[i])
			del arrays[i]
			a = arrays[i]
			best_match = find_closest_array(a, array_dict, tree)
		if best_match.extant: # If the closest match is a array then just join them and replace the existing array with the new node
			current_parent = array_dict[tree_child_dict[best_match.id].parent_node.taxon.label]
			results = replace_existing_array(
				best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
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
				best_match, a, current_parent, tree, all_arrays, node_ids, node_count, array_dict, tree_child_dict, seed
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
						node_array = count_parsimony_events(node_array, parent_array, array_dict, tree, True)
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							node_array.distance += node_array.events[k] * v
						node.edge_length = node_array.distance

			tree.reroot_at_node(tree.seed_node, update_bipartitions=False) # Need to reroot at the seed so that RF distance works
			return array_dict, tree, arrays

	
print("Running with the following command:\n{}\n".format(" ".join(sys.argv)))

event_costs = { 
				"acquisition" : args.acquisition,
				"indel" : args.indel,
				"repeated_indel" : args.rep_indel,
				"duplication": args.duplication,
				"trailer_loss": args.trailer_loss
				}

# hex values from this website http://phrogz.net/css/distinct-colors.html

Cols_hex_27 = ['#fd5925', '#dbc58e', '#008d40', '#304865', '#934270', '#f7b8a2', '#907500', '#45deb2', '#1f4195', '#d67381', '#8e7166', '#afb200', '#005746', '#a598ff', '#8f0f1b', '#b96000', '#667f42', '#00c7ce', '#9650f0', '#614017', '#59c300', '#1a8298', '#b5a6bd', '#ea9b00', '#bbcbb3', '#00b0ff', '#cd6ec6']

#hex values from https://mokole.com/palette.html

Cols_hex_40 = ["#696969","#556b2f","#a0522d","#800000","#006400","#808000","#483d8b","#3cb371","#008080","#bdb76b","#4682b4","#000080","#9acd32","#32cd32","#daa520","#7f007f","#ff4500","#00ced1","#ff8c00","#c71585","#0000cd","#00ff00","#9400d3","#dc143c","#00bfff","#f4a460","#adff2f","#da70d6","#ff00ff","#1e90ff","#db7093","#fa8072","#ffff54","#dda0dd","#7b68ee","#afeeee","#98fb98","#7fffd4","#ffe4c4","#ffc0cb"]

Cols_tol = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"]

Cols_hex_12 = ["#07001c", "#ff6f8d", "#4c62ff", "#92ffa9", "#810087", "#bcffe6", "#490046", "#00c8ee", "#b53900", "#ff8cf7", "#5b5800", "#14d625"]

# Generate strings to assign as internal node_IDs (This makes 702)

node_ids = ["Int_" + i for i in ascii_lowercase]
if len(args.arrays_to_join) > 27: # Maximum internal nodes in tree is n-2 so only need more than 26 if n >= 28
	node_ids += ["Int_" + "".join(i) for i in product(ascii_lowercase, repeat=(len(args.arrays_to_join)//26)+1)]


array_spacers_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_spacers_dict[bits[0]] = bits[2:]

if args.arrays_to_join:
	arrays = [Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
	labels = args.arrays_to_join
else:
	arrays = [Array(i, array_spacers_dict[i]) for i in array_spacers_dict.keys()]
	labels = list(array_spacers_dict.keys())

all_arrays = [array.spacers for array in arrays]


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
			addition_order = [Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
		else:
			addition_order = array_choices[i]

		if args.fix_order:
			if not args.arrays_to_join:
				print("You must provide the order you want to fix when using the fixed order option!\n\nABORTING.")
				sys.exit()
			addition_order = [Array(i, array_spacers_dict[i]) for i in args.arrays_to_join]
		else:
			addition_order = array_choices[i]
		try:
			array_dict, tree = build_tree_single(addition_order, taxon_namespace, best_score, all_arrays)
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
	options = [(array_list, taxon_namespace, all_arrays) for array_list in array_choices]
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
				plot_tree(good_tree, best_arrays[n], filename)
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
			print("Saving image of tree with array diagrams to {}\n".format(args.output_tree))
			plot_tree(best_tree, best_arrays, args.output_tree)
		if args.output_arrays:
			print("Saving details of arrays to {}\n".format(args.output_arrays))
			with open(args.output_arrays, 'w') as fout:
				fout.write('\n'.join(["{}\t{}".format(k," ".join(v.spacers)) for k,v in best_arrays.items()]))
except Exception as e:
	exc_type, exc_obj, exc_tb = sys.exc_info()
	exc_tb.print_exception()

