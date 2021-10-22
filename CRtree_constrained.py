#! /usr/bin/env python3

import dendropy
import sys
import CRISPR_mp
from string import ascii_lowercase
from itertools import product
from math import ceil, log


intree = sys.argv[1]
array_file = sys.argv[2]
genome_array_file = sys.argv[3]


event_costs = { 
	"acquisition" : 1,
	"indel" : 10,
	"repeated_indel" : 50,
	"duplication": 1,
	"trailer_loss": 1,
	"no_ident": 100
	}

array_spacers_dict = {}
with open(array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_spacers_dict[bits[0]] = bits[2:]

genome_array_dict = {}
with open(genome_array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		if bits[0].lower() == "outgroup":
			outgroup_taxon = bits[1]
			continue
		genome_array_dict[bits[1]] = bits[0]
		


array_dict = {}
for array in genome_array_dict.values():
	array_dict[array] = CRISPR_mp.Array(array, array_spacers_dict[array], extant=True)
all_arrays = [array.spacers for array in array_dict.values()]

tree = dendropy.Tree.get(path=intree, schema="newick")


# If tree not rooted using outgroup node then root it
outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

# Now remove outgroup as it isn't used for array analysis
tree.prune_taxa_with_labels([outgroup_taxon])

print(tree.as_ascii_plot(plot_metric='length'))

node_ids = ["Anc " + i for i in ascii_lowercase]

if len(genome_array_dict) > 27: # Number of internal nodes in tree is n-1 so only need more than 26 if n >= 28
	node_ids += ["Anc " + "".join(i) for i in product(ascii_lowercase, repeat=(ceil(log(len(genome_array_dict), 26))))]

taxon_namespace = dendropy.TaxonNamespace(list(genome_array_dict.values()) + node_ids)


for leaf in tree.leaf_node_iter():
	old_label = leaf.taxon.label
	leaf.taxon = taxon_namespace.get_taxon(genome_array_dict[old_label])

print(tree.as_ascii_plot(plot_metric='length'))


node_count = 0
for node in tree.seed_node.postorder_iter():
	if node.taxon:
		a = node.taxon.label 
		if len(node.sibling_nodes()) != 0:
			sister = node.sibling_nodes()[0]
			parent = node.parent_node
			if sister.taxon:
				b = sister.taxon.label
				if not parent.taxon: # Only need to add the ancestor once.
					results = CRISPR_mp.resolve_pairwise_parsimony(array_dict[a], array_dict[b], all_arrays, array_dict, node_ids, node_count, tree, event_costs)
					# if results == "No_ID":
					# 	Incomplete_tree = True
					# 	break
					array_dict[a], array_dict[b], ancestor = results
					node_count+=1

					for n, n_id in [(node, a), (sister, b)]:
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							array_dict[n_id].distance += array_dict[n_id].events[k] * v
						# n.edge_length = array_dict[n_id].distance

					array_dict[ancestor.id] = ancestor
					node.parent_node.taxon = taxon_namespace.get_taxon(ancestor.id)
					# node.parent_node.edge_length = 0 # Start the ancestor with 0 branch length


print(tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))

