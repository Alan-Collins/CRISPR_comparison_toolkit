#! /usr/bin/env python3

import sys
from math import ceil, log
from collections import Counter

import dendropy

import CRISPRtree
from cctk import (
	colour_schemes,
	file_handling,
	tree_operations
	)


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

array_spacers_dict = file_handling.read_array_file(array_file)


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
	array_dict[array] = CRISPR_mp.Array(
		array, array_spacers_dict[array], extant=True)
all_arrays = [array.spacers for array in array_dict.values()]


all_spacers = []
for array in array_dict.keys():
	all_spacers += array_spacers_dict[array]
non_singleton_spacers = [
	spacer for spacer, count in Counter(all_spacers).items() if count >1]

colours = colour_schemes.choose_col_scheme(len(non_singleton_spacers))
# build a dictionary with colours assigned to each spacer.
spacer_cols_dict  = {}

for i, spacer in enumerate(sorted(non_singleton_spacers)):
	spacer_cols_dict[spacer] = colours[i]


tree = dendropy.Tree.get(path=intree, schema="newick")


# If tree not rooted using outgroup node then root it
outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

# Now remove outgroup as it isn't used for array analysis
tree.prune_taxa_with_labels([outgroup_taxon])


tree = tree_handling.scale_branches(tree, 10)



labels = [g for g in genome_array_dict.values()]
node_ids = tree_operations.create_internal_node_ids(
	len(genome_array_dict))


taxon_namespace = dendropy.TaxonNamespace(labels + node_ids)


for leaf in tree.leaf_node_iter():
	old_label = leaf.taxon.label
	leaf.taxon = taxon_namespace.get_taxon(genome_array_dict[old_label])


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
					results = CRISPR_mp.resolve_pairwise_parsimony(
						array_dict[a], 
						array_dict[b], 
						all_arrays, 
						array_dict, 
						node_ids, 
						node_count, 
						tree, 
						event_costs,
						)
					# if results == "No_ID":
					# 	Incomplete_tree = True
					# 	break
					array_dict[a], array_dict[b], ancestor = results
					node_count+=1

					for n, n_id in [(node, a), (sister, b)]:
						for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
							d = array_dict[n_id].events[k]*v
							array_dict[n_id].distance += d
						# n.edge_length = array_dict[n_id].distance

					array_dict[ancestor.id] = ancestor
					node.parent_node.taxon = taxon_namespace.get_taxon(
						ancestor.id)
					# node.parent_node.edge_length = 0 # Start the ancestor with 0 branch length

# Repeat iteration now that tree is built to add repeat indels.
for node in tree:
	if node.level() != 0:
		parent = node.parent_node				
		array_dict[node.taxon.label].reset()
		array_dict[parent.taxon.label].reset()
		array_dict[node.taxon.label] = CRISPR_mp.count_parsimony_events(
			array_dict[node.taxon.label], 
			array_dict[parent.taxon.label],
			array_dict, 
			tree, 
			True,
			)
		for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
			d = array_dict[node.taxon.label].events[k]*v
			array_dict[node.taxon.label].distance += d
			# node.edge_length = array_dict[node.taxon.label].distance

print(tree.as_string("newick"))

print(tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))

CRISPR_mp.plot_tree(tree, 
	array_dict, 
	"test_out.png", 
	spacer_cols_dict, 
	branch_lengths=False, 
	emphasize_diffs=True,
	)
