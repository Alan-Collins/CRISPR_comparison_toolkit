#! /usr/bin/env python3

import dendropy
import sys
import CRISPR_mp
from string import ascii_lowercase
from itertools import product
from math import ceil, log
from collections import Counter

def scale_branches(branch, all_branches):
	max_branch = max(all_branches)
	branch = int(9*branch/max_branch)+1
	return branch


Cols_hex_27 = ["#fd5925", "#dbc58e", "#008d40", "#304865", "#934270", "#f7b8a2", "#907500", "#45deb2", "#1f4195", "#d67381", "#8e7166", "#afb200", "#005746", "#a598ff", "#8f0f1b", "#b96000", "#667f42", "#00c7ce", "#9650f0", "#614017", "#59c300", "#1a8298", "#b5a6bd", "#ea9b00", "#bbcbb3", "#00b0ff", "#cd6ec6"]


Cols_tol = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255"]

Cols_hex_12 = ["#07001c", "#ff6f8d", "#4c62ff", "#92ffa9", "#810087", "#bcffe6", "#490046", "#00c8ee", "#b53900", "#ff8cf7", "#5b5800", "#14d625"]


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


all_spacers = []
for array in array_dict.keys():
	all_spacers += array_spacers_dict[array]
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


tree = dendropy.Tree.get(path=intree, schema="newick")


# If tree not rooted using outgroup node then root it
outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

# Now remove outgroup as it isn't used for array analysis
tree.prune_taxa_with_labels([outgroup_taxon])

edge_lens = []
for node in tree:
	if node.edge_length != None:
		edge_lens.append(node.edge_length)
	



for node in tree:
	if node.edge_length != None:
		node.edge_length = scale_branches(node.edge_length, edge_lens)
	else:
		node.edge_length = 0


node_ids = ["Anc " + i for i in ascii_lowercase]

if len(genome_array_dict) > 27: # Number of internal nodes in tree is n-1 so only need more than 26 if n >= 28
	node_ids += ["Anc " + "".join(i) for i in product(ascii_lowercase, repeat=(ceil(log(len(genome_array_dict), 26))))]

taxon_namespace = dendropy.TaxonNamespace(list(genome_array_dict.values()) + node_ids)


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

# Repeat iteration now that tree is built to add repeat indels.
for node in tree:
	if node.level() != 0:
		parent = node.parent_node				
		array_dict[node.taxon.label].reset()
		array_dict[parent.taxon.label].reset()
		array_dict[node.taxon.label] = CRISPR_mp.count_parsimony_events(array_dict[node.taxon.label], array_dict[parent.taxon.label], array_dict, tree, True)
		for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
			array_dict[node.taxon.label].distance += array_dict[node.taxon.label].events[k] * v
			# node.edge_length = array_dict[node.taxon.label].distance

print(tree.as_string("newick"))

print(tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))

CRISPR_mp.plot_tree(tree, array_dict, "test_out.png", spacer_cols_dict, branch_lengths=False, emphasize_diffs=True)