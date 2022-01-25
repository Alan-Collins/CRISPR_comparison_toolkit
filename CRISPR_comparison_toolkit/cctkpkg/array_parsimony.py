from collections import Counter
import copy

import sys

from . import sequence_operations

class Array():
	"""Store information about arrays and aid in their comparison.
	
	Attributes:
		id (str): Identifier for this array.
		extant (bool): A boolean indicating if the array is extant in 
			our dataset or if it was inferred as a hypothetical 
			ancestral state.
		modules (list): A list of the contiguous blocks of spacers with 
			common features 
			(e.g. consecutive spacers that are absent in aligned array).
		spacers (list): A list of the spacers in this array.
		aligned (list): A list of spacers in aligned format relative to 
			another array
		module_lookup (dict): A dict with indices as keys and 
			SpacerModule instances as values where a given 
			SpacerModule instance will be pointed to by all the 
			indices at which it is located.
		distance (int): Parsimony distance from the hypothetical 
			ancestral state of this array.
		events (dict): Record what each event is for later addition of 
			event parsimony weights.
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
						"trailer_loss": 0,
						"no_ident": 0,
						}

	def sort_modules(self):
		"""Sort SpacerModules by their indices within the array.
		"""
		self.modules.sort(key=lambda x: int(x.indices[0]))

	def reset(self):
		"""Reset distance, modules, and event counts.
		"""
		self.distance = 0
		self.modules = []
		self.module_lookup = {}
		self.events = { 
						"acquisition" : 0,
						"indel" : 0,
						"repeated_indel" : 0,
						"duplication": 0,
						"trailer_loss": 0,
						"no_ident" : 0,
						}


class SpacerModule():
	"""
	Class to store information about spacers in CRISPR arrays.
	
	Attributes:
		type (str): A string indicating the nature of this module. 
			Possible types:
			N.B. leader region ends once first identical spacer is 
			found.
			- aqcuisition: Spacers at the leader end of one array where 
			no (or different) spacers are found in the other array. 
			- no_acquisition: Gap at leader end of array when other 
			array has spacers not found in this array.
			- indel: non-leader region with spacers present in one 
			array but a gap or different spacers in the other array
			- shared: Region where both arrays have the same spacers.
		spacers (list): A list of the spacer IDs in this module.
		indices (list): A list of the indices in the respective array 
			where this spacer module is located.
		partner (str): If the type of this module is "repeated_indel", 
		then this stores the ID(s) of the array(s) identified as the 
		other instance of this indel(s).
	"""
	def __init__(self):
		self.type = ""
		self.spacers = []
		self.indices = []
		self.partner = []


def find_modules(array1, array2, parent_comparison=False):
	"""
	Identify contiguous stretches of spacers that share an evolutionary
	property in terms of how they differ from the comparator array (e.g.
	indel region). Assigns types to each module according to the types
	defined in the SpacerModule class docstring.
	
	Args:
		array1 (Array class instance): The first array to be compared.
		array2 (Array class instance): The second array to be compared.
		parent_comparison (Bool): Is this a comparison with array1 as
			the child array of array2?
	
	Returns:
		(tuple of Array class instances) The provided arrays with module
		information added to the .module attribute.
	"""

	array1.aligned, array2.aligned = sequence_operations.needle(
		array1.spacers, array2.spacers)

	# If this isn't the first comparison these arrays have been part of
	# then need to reset module list and lookup dict.
	array1.reset()
	array2.reset()

	leader = True
	gap = False
	mismatch = False

	module1 = SpacerModule()
	module2 = SpacerModule()

	# If the two arrays being compared have no shared spacers, a
	# high-cost event should be assigned.

	if not any([a==b for a,b in zip(array1.aligned, array2.aligned)]):
		module1.type = "no_ident"
		module1.indices = [n for n in range(len(array1.aligned))]
		module1.spacers = [a for a in array1.aligned]
		for k in module1.indices:
			array1.module_lookup[k] = module1
		
		module2.type = "no_ident"
		module2.indices = [n for n in range(len(array2.aligned))]
		module2.spacers = [b for b in array2.aligned]
		for k in module2.indices:
			array2.module_lookup[k] = module1

		return array1, array2


	if parent_comparison:
		array1, array2 = modify_parent_child_alignment(array1, array2)

	# Identify modules in aligned arrays

	for n, (a,b) in enumerate(zip(array1.aligned, array2.aligned)):

		# Leader end processing

		if leader:
			if a == '-' or b == '-': 
				# If this is a child-parent comparison, use specialized
				# function for this module

				if parent_comparison:
					module1, array1, module2, array2 = process_parent_leader(
					module1, array1, module2, array2, n, a, b)
					continue

			# Indicates one array acquired spacers that the other didn't
				if a == '-':

					# Module1 processing for no_acqusition module
					if module1.type != "no_acquisition" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = SpacerModule()
					module1.type = "no_acquisition"
					module1.indices.append(n)
					module1.spacers.append(a)

					# Module2 processing for acquisition module
					# No need to check module type as gap penalty means
					# needle func will never return stretch like the
					# following:
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
						module2 = SpacerModule()
					module2.type = "no_acquisition"
					module2.indices.append(n)
					module2.spacers.append(b)

					module1.type = "acquisition"
					module1.indices.append(n)
					module1.spacers.append(a)

			elif a != b: # Mismatch at leader end probably means they've
				# both acquired spacers since ancestor
				# Indel also possible. Maybe implement parsimony
				# evaluation of that when comparing ancestor to other
				# arrays?

				if module1.type != "acquisition" and module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = SpacerModule()
				module1.type = "acquisition"
				module1.indices.append(n)
				module1.spacers.append(a)

				if module2.type != "acquisition" and module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = SpacerModule()
				module2.type = "acquisition"
				module2.indices.append(n)
				module2.spacers.append(b)

			else: # Other alternative is identical spacers and end of
				# leader region
				leader = False
				if module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = SpacerModule()
				module1.type = "shared"
				module1.indices.append(n)
				module1.spacers.append(a)
				if module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = SpacerModule()
				module2.type = "shared"
				module2.indices.append(n)
				module2.spacers.append(b)

		else:
			if a == b:
				if module1.type != "shared" and module1.type != "":
					array1.modules.append(module1)
					for k in module1.indices:
						array1.module_lookup[k] = module1
					module1 = SpacerModule()
				module1.type = "shared"
				module1.indices.append(n)
				module1.spacers.append(a)

				if module2.type != "shared" and module2.type != "":
					array2.modules.append(module2)
					for k in module2.indices:
						array2.module_lookup[k] = module2
					module2 = SpacerModule()
				module2.type = "shared"
				module2.indices.append(n)
				module2.spacers.append(b)
			else:
				# Check if duplication
				# Only call it a duplication if the comparator has the spacers, otherwise prefer calling this an insertion.
				if a != b and Counter(array2.aligned)[b] > 1 and b != '-' and b in array1.spacers:
					# If this spacer is in multiple copies and aligns with a gap it means the comparator array has no or fewer copies of this spacer.
					if module2.type != "duplication" and module2.type != "":
						array2.modules.append(module2)
						for k in module2.indices:
							array2.module_lookup[k] = module2
						module2 = SpacerModule()
					module2.type = "duplication"
					module2.indices.append(n)
					module2.spacers.append(b)

					if module1.type != "duplication" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = SpacerModule()
					module1.type = "duplication"
					module1.indices.append(n)
					module1.spacers.append(a)
				elif a != b and Counter(array1.aligned)[a] > 1 and a != '-' and a in array2.spacers:
					if module1.type != "duplication" and module1.type != "":
						array1.modules.append(module1)
						for k in module1.indices:
							array1.module_lookup[k] = module1
						module1 = SpacerModule()
					module1.type = "duplication"
					module1.indices.append(n)
					module1.spacers.append(a)

					if module2.type != "duplication" and module2.type != "":
						array2.modules.append(module2)
						for k in module2.indices:
							array2.module_lookup[k] = module2
						module2 = SpacerModule()
					module2.type = "duplication"
					module2.indices.append(n)
					module2.spacers.append(b)

				else:
					if a != '-' and b != '-':
						# Module1 processing for indel module where mismatched
						if module1.type not in ["indel_mm", "indel_gap"] and module1.type != "":
							array1.modules.append(module1)
							for k in module1.indices:
								array1.module_lookup[k] = module1
							module1 = SpacerModule()
						module1.type = "indel_mm"
						module1.indices.append(n)
						module1.spacers.append(a)
						# Module2 processing for indel module where mismatched
						if module2.type not in ["indel_mm", "indel_gap"] and module2.type != "":
							array2.modules.append(module2)
							for k in module2.indices:
								array2.module_lookup[k] = module2
							module2 = SpacerModule()
						module2.type = "indel_mm"
						module2.indices.append(n)
						module2.spacers.append(b)
					else:
						if n == len(array1.aligned)-1 and module1.type not in ["indel_mm", "indel_gap"]: # If this is the last spacer and not a continuation of an indel, call it trailer loss
							array1.modules.append(module1)
							for k in module1.indices:
								array1.module_lookup[k] = module1
							module1 = SpacerModule()
							module1.type = "trailer_loss"
							module1.indices.append(n)
							module1.spacers.append(a)
							array2.modules.append(module2)
							for k in module2.indices:
								array2.module_lookup[k] = module2
							module2 = SpacerModule()
							module2.type = "trailer_loss"
							module2.indices.append(n)
							module2.spacers.append(b)
						else:
							# Module1 processing for indel module where one is gap
							if module1.type not in ["indel_mm", "indel_gap"] and module1.type != "":
								array1.modules.append(module1)
								for k in module1.indices:
									array1.module_lookup[k] = module1
								module1 = SpacerModule()
								module1.type = "indel_gap"
							module1.indices.append(n)
							module1.spacers.append(a)
							# Module2 processing for indel module where one is gap
							if module2.type not in ["indel_mm", "indel_gap"] and module2.type != "":
								array2.modules.append(module2)
								for k in module2.indices:
									array2.module_lookup[k] = module2
								module2 = SpacerModule()
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


def modify_parent_child_alignment(child, parent):
	""" When comparing known child-parent arrays, events may occur that
	CRISPRtree-adapted alignment will not catch. The alignment must be
	modified for that case.
	"""

	# find first shared spacer location

	for n, (a,b) in enumerate(zip(child.aligned, parent.aligned)):
		if a == b:
			shared_idx = n
			break

	# Count number of non-gap spacers in each leader end
	n_spacers_child = 0
	n_spacers_parent = 0
	# and store the non-gap leader spacers
	child_lead_sps = []
	parent_lead_sps = []

	for a,b in zip(child.aligned[:shared_idx], parent.aligned[:shared_idx]):
		if a != '-':
			n_spacers_child += 1
			child_lead_sps.append(a)
		if b != '-':
			n_spacers_parent += 1
			parent_lead_sps.append(b)

	# Align leaders so that differences are scored as gaps.
	# i.e.
	# A B C - - -
	# - - - D E F

	child.aligned = (child_lead_sps # leader spacers in child
		+ ['-']*n_spacers_parent # gap character per leader sp in parent
		+ child.aligned[shared_idx:]) # Non-leader child spacers

	parent.aligned = (['-']*n_spacers_child
		+ parent_lead_sps
		+ parent.aligned[shared_idx:])

	return child, parent


def process_parent_leader(
					module1, array1, module2, array2, n, a, b):
	""" Process leader differences between ancestral and child array.

	If an ancestral array has spacers at its leader end that are missing
	in a child array then these spacers must have been deleted in the
	child array. This should only be seen for arrays produced by
	evolve_array.py as CRISPRtree should never predict this event.


	"""

	if a != '-':
		# Child leader acquisition which will come first in alignment
		module1.type = "acquisition"
		module1.indices.append(n)
		module1.spacers.append(a)

		module2.type = "no_acquisition"
		module2.indices.append(n)
		module2.spacers.append(b)

		return module1, array1, module2, array2

	elif b != '-':
		# loss of parent leader spacers in child.

		# Module1 processing for indel_gap module
		if module1.type != "indel_gap" and module1.type != "":
			array1.modules.append(module1)
			for k in module1.indices:
				array1.module_lookup[k] = module1
			module1 = SpacerModule()
		module1.type = "indel_gap"
		module1.indices.append(n)
		module1.spacers.append(a)

		if module2.type != "indel_gap" and module2.type != "":
			array2.modules.append(module2)
			for k in module2.indices:
				array2.module_lookup[k] = module2
			module2 = SpacerModule()
		module2.type = "indel_gap"
		module2.indices.append(n)
		module2.spacers.append(b)

		return module1, array1, module2, array2


	else:
		print("Unexpected parent-child alignment encountered while processing\
			 arrays {} and {}. Exiting...".format(child.id, parent.id))
		sys.exit()

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
		if k != '-' and k in ancestor.aligned:
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


def count_parsimony_events(child, ancestor, array_dict, tree,
	parent_comparison=False):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of
			the child array.
		array_dict (dict): Dict of Array class instances with
			information about the nodes of your tree.
		tree (Deondropy Tree class instance): The tree in which the
			arrays are located.	
		parent_comparison (bool): Are the arrays being compared
			parent-child nodes in the tree? E.g. if you are trying to
			find the best matching array in the tree this is False. If
			you are counting events between an array and its ancestor
			array this is True

	Returns:
		(Array class instance) The child Array class instance with
			parsimony events added to the .events dict attribute.
	"""

	child, ancestor = find_modules(child, ancestor, parent_comparison)

	# Process modules to count events since ancestor

	# First look for duplications

	dupes = find_dupes(child, ancestor)
	# Flatten list to check against later
	dupe_indices = [x for y in dupes for x in y]

	child.events['duplication'] = int(len(dupes) / 2)

	idx = 0
	while idx <= max(child.module_lookup.keys()):
		if idx in dupe_indices:
			idx += 1
			continue
		mod = child.module_lookup[idx]
		if mod.type == 'acquisition' or mod.type == 'no_acquisition': 
		# Checking no acquisition allows comparison of child-child as
		# well as child-ancestor arrays.
			if parent_comparison: 
				# If comparing to a nodes ancestor, then deletions of
				# spacers from the leader end may look like acquisition.
				# This must be checked to make sure deletions are not
				# being miscalled as acquisitions.
				child = identify_repeat_indels(
					child, ancestor, array_dict, mod,
					ancestor.module_lookup[idx], tree)
			else:
				child.events['acquisition'] += len(mod.indices)
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
			continue
		elif mod.type == 'indel_gap' or mod.type == "indel_mm":
			if parent_comparison:
				if not all([s == '-' for s in mod.spacers]):
					child = identify_repeat_indels(
						child, ancestor, array_dict, mod,
						ancestor.module_lookup[idx], tree)
				else:
					child.events['indel'] += 1
			else:
				child.events['indel'] += 1
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
		elif mod.type == 'duplication':
			child.events['duplication'] += 1
			idx = mod.indices[-1] + 1 # Skip the rest of this module.
		elif mod.type == "trailer_loss":
			# Only score as trailer loss if the ancestor has a spacer
			# and the child doesn't. If it's mismatch or child has
			# spacer then it's indel
			if mod.spacers[0] == '-':
				child.events['trailer_loss'] += 1
			else:
				mod.type = 'indel_mm'
				child.module_lookup[idx] = mod
				child.events['indel'] += 1
			idx = mod.indices[-1] + 1
		elif mod.type == "no_ident":
			child.events['no_ident'] += 1
			idx = mod.indices[-1] + 1
		else:
			idx += 1
	return child


def identify_repeat_indels(child, ancestor, array_dict, module,
	ancestor_module, tree):
	"""
	Args:
		child (Array class instance): The child array.
		ancestor (Array class instance): The hypothetical ancestor of the child array.
		array_dict (dict): Dict of Array class instances with information about the nodes of your tree.
		module (SpacerModule class instance): The indel module to be processed.
		ancestor_module (SpacerModule class instance): The module of the parent of the array to be processed.
		tree (Deondropy Tree class instance): The tree in which the arrays are located.
	
	Returns:
		(Array class instance) The child array, modified to contain newly identified SpacerModule instances if any are found.
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

			repeat_search = True
			original_spacers_to_check = copy.deepcopy(module.spacers) 
			spacers_to_check = copy.deepcopy(module.spacers) # Store which spacers to look for and an original copy to maintain index information
			array_IDs_of_concern = non_ancestor_children_ids + other_child_children_ids
			arrays_of_concern = non_ancestor_children_spacers + other_child_children_spacers
			while repeat_search:
				longest_match = 0
				longest_indices = []
				arrays_to_check = []
				# If lost spacers are found in either of these two groups it indicates independent gain of the same spacers in different lineages.
				for spacer in spacers_to_check:
					for a, array in enumerate(arrays_of_concern):
						if spacer in set(array):
							arrays_to_check.append(a)
				if len(arrays_to_check) == 0:
					repeat_search = False 
				else:
					for a in set(arrays_to_check):
						x = Array("x", spacers_to_check)
						y = Array("y", arrays_of_concern[a])
						x.aligned, y.aligned = sequence_operations.needle(x.spacers, y.spacers) # Find where the match is by aligning
						x,y = find_modules(x,y) # Then identify any shared regions
						for m in x.modules: # Pull out the shared modules
							if m.type == "shared":
								if len(m.indices) > longest_match:
									longest_match = len(m.indices)
									longest_indices = [original_spacers_to_check.index(s) for s in m.spacers if s in original_spacers_to_check]
									longest_spacers = [s for s in spacers_to_check if s in m.spacers]
									partner = [array_IDs_of_concern[a]]
									partner_extant = array_dict[array_IDs_of_concern[a]].extant
								elif len(m.indices) == longest_match and set(m.spacers) == set(longest_spacers):
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
					new_rep_indel_mod = SpacerModule()
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
					new_ac_mod = SpacerModule()
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
				if spacers_to_check != module.spacers: # If we haven't found anything, no need to update the existing indel. This check is passed if we found something
					if len(spacer_indices) > 0: # Then find consecutive runs of spacers. Each must be an indel
						new_indel_mod = SpacerModule() # Make a new SpacerModule to store the indel.
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
									new_indel_mod = SpacerModule() # Make a new SpacerModule to store the next indel.
									new_indel_mod.type = "indel"
									indels += 1
								last_n = n
								new_indel_mod.indices.append(module.indices[n])
								new_indel_mod.spacers.append(module.spacers[n])
							else: # Otherwise move to the next number
								last_n = n 
								new_indel_mod.indices.append(module.indices[n])
								new_indel_mod.spacers.append(module.spacers[n])
				else: # If we haven't found anything then if this was an indel then set the indel count to 1 and move on. If it was an acquisition then the count has already been added so just move on.
					if module.type in ["indel_gap", "indel_mm"]:
						indels = 1
					elif module.type == "acquisition":
						pass
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
