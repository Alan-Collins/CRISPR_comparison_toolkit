from string import ascii_lowercase
from math import ceil, log
from itertools import product


def create_internal_node_ids(n_leaves):
	"""Make list of internal node IDs for tree.

	As the number of internal nodes in a bifurcating tree = n-1,
	if fewer than 28 internal nodes then they can be uniquely
	identified using single letters. For larger numbers of internal 
	nodes multiple character identifiers are needed.
	"""
	node_ids = ["Anc " + i for i in ascii_lowercase]
	# Number of internal nodes in tree is n-1 
	# so only need more than 26 if n >= 28
	if n_leaves > 27: 
		node_ids += ["Anc " + "".join(i) for i in product(
			ascii_lowercase, repeat=(ceil(log(n_leaves, 26))))]
	return node_ids