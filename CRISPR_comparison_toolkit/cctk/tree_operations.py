from string import ascii_lowercase
from math import ceil, log
from itertools import product


def create_internal_node_ids(n_leaves, prefix="", chars="letters"):
	"""Make list of internal node IDs for tree.

	given the number of leaf nodes in a tree, generates enough IDs to
	uniquely identify all internal nodes of a bifurcating tree.

	Args:
		n_leaves (int): The number of leaf nodes in the tree.
		prefix (str, optional): Prefix with which all node IDs should begin.
		chars (str, optional): Options "letters", "numbers". The kind
			of character that should be used to identify each node.
	Returns:
		list of str: List of unique node IDs.

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


def scale_branches(branch, all_branches, max_len=10):
	max_branch = max(all_branches)
	scale = max_len/max_branch
	branch *= scale
	return branch
