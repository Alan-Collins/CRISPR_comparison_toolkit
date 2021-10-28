from string import ascii_lowercase
from math import ceil, log
from itertools import product


def create_internal_node_ids(n_leaves):
	"""Make list of internal node IDs for tree."""
	node_ids = ["Anc " + i for i in ascii_lowercase]
	# Number of internal nodes in tree is n-1 
	# so only need more than 26 if n >= 28
	if n_leaves > 27: 
		node_ids += ["Anc " + "".join(i) for i in product(
			ascii_lowercase, repeat=(ceil(log(n_leaves, 26))))]
	return node_ids