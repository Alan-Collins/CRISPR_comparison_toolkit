from string import ascii_lowercase
from math import ceil, log
from itertools import product
import dendropy


def create_internal_node_ids(n_leaves, prefix="", chars="letters"):
	"""Make list of internal node IDs for tree.

	given the number of leaf nodes in a tree, generates enough IDs to
	uniquely identify all internal nodes of a bifurcating tree.

	Args:
	  n_leaves (int): 
	  	The number of leaf nodes in the tree.
	  prefix (str, optional): 
	  	Prefix with which all node IDs should begin.
	  chars (str, optional): 
	  	Options "letters", "numbers". The kind of character that should 
	  	be used to identify each node.
	Returns:
	  list of str: 
	  	List of unique node IDs.

	Raises:
	  ValueError: If n_leaves is less than 2.
	  TypeError: If n_leaves is not int.
	  TypeError: If prefix is not str.
	  ValueError: If chars is not either "letters" or "numbers".
	"""
	if n_leaves < 2:
		raise ValueError("n_leaves must be 2 or more.")
	if type(n_leaves) is not int:
		raise TypeError(
			"n_leaves must be int, not {}.".format(type(n_leaves).__name__))
	if type(prefix) is not str:
		raise TypeError(
			"prefix must be str, not {}.".format(type(prefix).__name__))
	if chars not in ["letters", "numbers"]:
		raise ValueError('chars must be either "letters" or "numbers".')

	if chars == "letters":
		node_ids = [prefix + i for i in ascii_lowercase]
		# Number of internal nodes in tree is n-1 
		# so only need more than 26 if n >= 28
		# Numer of characters needed to generate N combinations 
		# is ceil(log26(N))
		if n_leaves > 27: 
			node_ids += [prefix + "".join(i) for i in product(
				ascii_lowercase, repeat=(ceil(log(n_leaves, 26))))]
		# End up with more IDs than needed so slice appropriate number.
		node_ids = node_ids[:n_leaves]
	else:
		# Start counting at 1 for general audience.
		node_ids = [prefix+str(i) for i in range(1, n_leaves+1)]
	return node_ids


def scale_branches(tree, max_len=10):
	"""Scales tree branch lengths linearly up to a user-defined maximum.
	
	Args:
	  tree (dendropy.Tree instance):
		Tree to scale
	  max_len (int or float):
		Value to set the longest branch to

	Returns:
	  tree:
	  	dendropy.Tree instance with scaled branch lengths.

	Raises:
	  TypeError: If tree is not a dendropy Tree class instance
	  ValueError: If max_len is not >0
	  TypeError: If max_len is not an int of float.
	"""
	edge_lens = []
	for node in tree:
		if node.edge_length != None:
			edge_lens.append(node.edge_length)

	max_branch = max(edge_lens)
	scale = max_len/max_branch

	for node in tree:
		if node.edge_length != None:
			new_branch = node.edge_length * scale
			# Handle ints vs floats
			if new_branch%1 == 0:
				node.edge_length = int(new_branch)
			else:
				node.edge_length = new_branch
		else:
			node.edge_length = 0
	
	 
	return tree


def yield_nodes(node):
	""" Work through a tree or subtree, yield the taxon labels of nodes
	"""
	if node.is_internal():
		# decide when to insert this node among its children
		# Should be about half way through its list of child nodes
		label_loc = int(len(node.child_nodes())/2)-1
		for i, n in enumerate(node.child_nodes()):
			for x in yield_nodes(n):
				yield x
			if i == label_loc:
				yield node.taxon.label

	else:
		yield node.taxon.label


def find_node_locs(tree, brlen_scale=0.5, branch_spacing=2.3):

	node_list = [n for n in yield_nodes(tree.seed_node)]

	node_locs = {}

	# Prepare tree for root distance queries
	tree.calc_node_root_distances(
		return_leaf_distances_only=False)

	for i, n in enumerate(node_list):
		node = tree.find_node_with_taxon_label(n)
		x_position = -node.root_distance*brlen_scale
		y_position = i*branch_spacing
		node_locs[n] = (x_position, y_position)

	return node_locs


def resolve_polytomies(tree):
	""" Collapse internal node branches with len == 0
	"""
	for node in tree.seed_node.postorder_internal_node_iter(
		exclude_seed_node=True):
		if node.edge_length == 0:
			# Transfer node children to parent node
			parent_children = [
			i for i in node.parent_node.child_nodes() if i != node]
			node.parent_node.set_child_nodes(
				node.child_nodes() + parent_children)
	
	return tree


def compare_to_trees(tree, comparator):
	if isinstance(comparator, list) or isinstance(
		comparator, dendropy.datamodel.treecollectionmodel.TreeList):
		rf_list = [
			dendropy.calculate.treecompare.weighted_robinson_foulds_distance(
				good_tree, tree) for good_tree in comparator]
		identical_list = [i == 0. for i in rf_list]
		return rf_list
	
	else:
		rf = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(
				comparator, tree)
		ident = True if rf != 0. else False
		return ident


def process_duplicate_arrays_contrain(tree,	genome_array_dict):
	new_genome_array_dict = {}
	# extract taxon namespace to add new entries later
	tns = tree.taxon_namespace

	for current_id, array_id_list in genome_array_dict.items():
		if len(array_id_list) == 1:
			array = array_id_list[0]
			new_id = current_id + ".{}".format(array)
			new_genome_array_dict[new_id] = array

			tax = dendropy.Taxon(new_id)
			tns.add_taxon(tax)

			node_mod = tree.find_node_with_taxon_label(current_id)
			node_mod.taxon = tax
			node_mod.label = new_id
			continue

		# Replace corresponding leaf node with a blank node
		# A new leaf will be added as a child of this node for each
		# array found in this genome
		node_mod = tree.find_node_with_taxon_label(current_id)
		node_mod.taxon = None

		for array in array_id_list:
			new_id = current_id + ".{}".format(array)
			new_genome_array_dict[new_id] = array

			# Add new taxon to namespace
			tax = dendropy.Taxon(new_id)
			tns.add_taxon(tax)
			# Create node using new taxon and add as child of node_mod
			node_mod.add_child(dendropy.Node(taxon=tax, label=new_id))

	return tree, new_genome_array_dict

