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

	# If the two arrays being compared have no shared spacers, a high-cost event should be assigned.

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
				# Only call it a duplication if the comparator has the spacers, otherwise prefer calling this an insertion.
				if a != b and Counter(array2.aligned)[b] > 1 and b != '-' and b in array1.spacers:
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
				elif a != b and Counter(array1.aligned)[a] > 1 and a != '-' and a in array2.spacers:
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
						if module1.type not in ["indel_mm", "indel_gap"] and module1.type != "":
							array1.modules.append(module1)
							for k in module1.indices:
								array1.module_lookup[k] = module1
							module1 = Spacer_Module()
						module1.type = "indel_mm"
						module1.indices.append(n)
						module1.spacers.append(a)
						# Module2 processing for indel module where mismatched
						if module2.type not in ["indel_mm", "indel_gap"] and module2.type != "":
							array2.modules.append(module2)
							for k in module2.indices:
								array2.module_lookup[k] = module2
							module2 = Spacer_Module()
						module2.type = "indel_mm"
						module2.indices.append(n)
						module2.spacers.append(b)
					else:
						if n == len(array1.aligned)-1 and module1.type not in ["indel_mm", "indel_gap"]: # If this is the last spacer and not a continuation of an indel, call it trailer loss
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
							if module1.type not in ["indel_mm", "indel_gap"] and module1.type != "":
								array1.modules.append(module1)
								for k in module1.indices:
									array1.module_lookup[k] = module1
								module1 = Spacer_Module()
								module1.type = "indel_gap"
							module1.indices.append(n)
							module1.spacers.append(a)
							# Module2 processing for indel module where one is gap
							if module2.type not in ["indel_mm", "indel_gap"] and module2.type != "":
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
