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
		"""Sort SpacerModules by their indices within the array."""
		self.modules.sort(key=lambda x: int(x.indices[0]))

	def reset(self):
		"""Set distance and event counts to zero."""
		self.distance = 0
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