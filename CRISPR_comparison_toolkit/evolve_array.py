#!/usr/bin/env python

import sys
import os
import argparse
import random
from copy import deepcopy

import dendropy

from cctk import array_parsimony, tree_operations, colour_schemes


def cmdline_args():

	p = argparse.ArgumentParser(
		description="evolve CRISPR arrays in silico. Starts with a single \
		array and performs a simple in silico evolution process."
		)
	p.add_argument(
		"-s", "--seed", type=int, metavar="", help="Set seed for consistency")
	p.add_argument(
		"-n", "--num-events", type=int, default=20, metavar="",
		help="How many events should be allowed to occur before the \
		simulation ends? Default = 20"
		)
	p.add_argument(
		"-o", "--outdir", type=str, default="./", metavar="",
		help="Directory in which output files should be written. \
		Defalt is your current directory"
		)

	parameters = p.add_argument_group('Evolution parameters', 
		"Specify the relative frequencies with which different events \
		should occur.")
	parameters.add_argument(
		"-i", "--initial-length", type=int, default=5, metavar="",
		help="Length of the starting array. Default = 5")
	parameters.add_argument(
		"-a", "--acquisition", type=int, default=75, metavar="",
		help="Relative frequency of spacer acquisitions. Default = 75")
	parameters.add_argument(
		"-t", "--trailer-loss", type=int, default=15, metavar="",
		help="Relative frequency of trailer spacer decay. Default = 15")
	parameters.add_argument(
		"-d", "--deletion", type=int, default=10, metavar="",
		help="Relative frequency of deletions of one or more spacers.\
		Default = 10")
	parameters.add_argument(
		"-l", "--loss-rate", type=int, default=50, metavar="",
		help="Rate (percent) at which arrays are lost after spawning a \
		descendant. Default = 50")

	plotting = p.add_argument_group('Plotting parameters', 
		"Set parameters for plotting the tree to file")
	plotting.add_argument(
		"-f", "--font-scale", type=float, default=1, metavar="",
		help="Font size is multiplied by this number. Default is to scale \
		based on the plot size")
	plotting.add_argument(
		"--font-override-labels", type=float, metavar="",
		help="If you want to force a label font size to be used rather than \
		using scaling to determine font size, provide it here")
	plotting.add_argument(
		"--font-override-annotations", type=float, metavar="",
		help="If you want to force a specific font size (pts) to be used for \
		annotations such as branch length labels, rather than using scaling \
		to determine font size, provide it here")
	plotting.add_argument(
		"-b", "--brlen-labels", action="store_true",
		help="Add branch length labels to plot.")
	plotting.add_argument(
		"--dpi", type=float, default=300, metavar="",
		help="Pixel density in pixels per inch (only relevant for PNG output.\
		Default = 300")
	plotting.add_argument(
		"--tree-width", type=float, default=3, metavar="",
		help="Width of plot in inches. Default = 3")
	plotting.add_argument(
		"--tree-height", type=float, default=3, metavar="",
		help="Height of plot in inches. Default = 3")
	plotting.add_argument(
		"--branch-weight", type=float, default=1, metavar="",
		help="Thickness of branch lines. Default = 1")
	plotting.add_argument(
		"--branch-spacing", type=float, default=1, metavar="",
		help="The vertical space between branches will be multiplied by this \
		number. Default = 1")
	plotting.add_argument(
		"--brlen-scale", type=float, default=1, metavar="",
		help="Branch lengths will be multiplied by this number. Default = 1")
	plotting.add_argument(
		"--align", action="store_true",
		help="Align array labels and cartoons rather than drawing them at \
		leaf tips and intenal nodes.")


	return p.parse_args()


def check_outdir(outdir):
	"""Check directory is a valid and add a trailing '/' if needed
	"""
	if not os.path.isdir(outdir):
		raise Exception(
			"{} is not a directory. Please make sure you ".format(outdir)
			+ "provide an existing directory as the output dir.")

	# Add / to outdir path if user forgot
	outdir = outdir + '/' if outdir[-1] != '/' else outdir

	return outdir


def tick(active_arrays, tree, tree_namespace, spacer_n, array_name, events_list,
	all_arrays, loss_rate):
	"""Choose an event and return the resulting array.
	
	Picks an existing array to modify and an event at random. Creates a
	modified copy of the chosen array according to the event chosen and
	returns the modified copy of the array

	Args:
	  active_arrays (list of Array):
		List of Arrays that are still available to use as source arrays.
	  tree (Dendropy.Tree):
	    Tree instance containing array histories.
	  tree_namespace (Dendropy.TaxonNamespace):
	    namespace for tree nodes.
	  spacer_n (int):
		Next spacer ID to use.
	  array_name (int):
	    The ID to assign to newly create Array instances.
	  events_list (list):
	    List of events to choose from where events are repeated in list
	    to control chance of their being selected.
	  all_arrays (list):
	    list of all Array instances
	  loss_rate (int):
	    Percent rate at which to remove old arrays from the 
	    active_arrays after they are the source of a new array

	Returns:
	  A modified form of the chosen array as an Array class instance and
	  other updates objects.
	"""

	source_array = random.choice(active_arrays)
	event = random.choice(events_list)

	# Don't allow deletion if array is only 1 spacer long

	if len(source_array.spacers) == 1:
		while event == "Deletion":
			event = events_dict[random.randint(1,len(events_dict))]

	array = array_parsimony.Array(str(array_name), 
		spacers=[i for i in source_array.spacers],
		extant=False)
	array_name+=1

	if event == "Acquisition":
		array, spacer_n = do_acquisition(array, spacer_n)
	elif event == "Trailer_loss":
		array = do_trailer_loss(array)
	else:
		array = do_deletion(array)

	# Add new array to active_arrays and all_arrays lists
	active_arrays.append(array)
	all_arrays.append(array)

	# Check whether to remove source array from active_arrays
	x = random.randint(1,100)
	if x <= loss_rate:
		del_idx = active_arrays.index(source_array)
		del active_arrays[del_idx]

	# Add array relationships to tree

	new_node = dendropy.Node()
	new_node.taxon = tree_namespace.get_taxon(array.id)

	tree.find_node_with_taxon_label(source_array.id).add_child(new_node)


	return active_arrays, tree, spacer_n, array_name, all_arrays


def do_acquisition(array, spacer_n):
	array.spacers.insert(0,spacer_n)
	spacer_n+=1

	return array, spacer_n


def do_trailer_loss(array):
	del array.spacers[-1]

	return array
	

def do_deletion(array):
	# Define parameters for a normal distribution
	# This distribution will be used to select random spacers that
	# tend to be close to the middle of the array
	array_len = len(array.spacers)
	# Define the central point of the array
	mean = array_len/2

	# Set the level of variation from the middle of the array
	stddev = array_len/4

	a = pick_normal_index(array_len, mean, stddev)
	b = pick_normal_index(array_len, mean, stddev)

	# To ensure a deletion takes place, but not the whole array
	while b == a or (a in [0, array_len] and b in [0, array_len]):
		b = pick_normal_index(array_len, mean, stddev)

	if a > b:
		a,b = b,a

	spacers_to_del = [i for i in range(a,b)]

	array.spacers = [
		i for n,i in enumerate(array.spacers) if n not in spacers_to_del]

	return array


def pick_normal_index(max, mean, stddev):

	while True: # Keep going until a valid index is chosen
		idx = int(random.normalvariate(mean, stddev))
		if 0 <= idx <= max:
			return idx


def summarise_tree(tree):
	""" Removes nodes with only a single descendant from the given tree.
	"""

	new_tree = dendropy.Tree(tree)

	for node in new_tree.postorder_node_iter():
		if node == new_tree.seed_node: # Skip seed node
			continue
		if len(node.child_nodes()) == 1:
			current_child = node.child_nodes()[0]
			# Retrive existing children of parent
			parent_children = node.parent_node.child_nodes()
			# Find index of current node in parents child_nodes
			idx = parent_children.index(node)
			parent_children[idx] = current_child
			node.parent_node.set_child_nodes(parent_children)

	return new_tree


def remove_lost_spacers(array_dict, active_spacers):
	"""Remove spacers from array dict that are not present in final arrays"""

	new_dict = deepcopy(array_dict)

	for array in new_dict.values():
		array.spacers = [i for i in array.spacers if i in active_spacers]

	return new_dict


def main(args):

	outdir = check_outdir(args.outdir)

	if args.seed:
		random.seed(args.seed)	

	spacer_n = 1 # Track ID of spacers
	array_name = 0
	active_arrays = []
	events_list = []

	init_array = array_parsimony.Array(str(array_name), extant=False)
	array_name+=1

	# Initialize starting array
	for _ in range(args.initial_length):
		init_array.spacers.insert(0,spacer_n)
		spacer_n+=1
	all_arrays = [init_array]
	active_arrays = [init_array]

	# build list of events where repetitions of event types controls
	# chance that the event will be selected.
	events_list += (['Acquisition']*args.acquisition
		+ ['Trailer_loss']*args.trailer_loss
		+ ['Deletion']*args.deletion)

	# Initialize tree
	tree_namespace = dendropy.TaxonNamespace(
		[str(i) for i in range(args.num_events+1)])
	tree = dendropy.Tree(taxon_namespace=tree_namespace)

	first_node = dendropy.Node()
	first_node.taxon = tree_namespace.get_taxon('0')
	tree.seed_node.add_child(first_node)

	# Run evolution model
	for _ in range(args.num_events):
		active_arrays, tree, spacer_n, array_name, all_arrays = tick(
			active_arrays, tree, tree_namespace, spacer_n, array_name, events_list,
			all_arrays, args.loss_rate)

	# for a in all_arrays:
	# 	print(str(a.id) + '\t' + ' '.join([str(i) for i in a.spacers]))

	# Remove intermediate nodes from tree
	new_tree = summarise_tree(tree)

	# Add arbitrary branch lengths for plotting purposes.
	for node in new_tree.postorder_node_iter():
		# Skip seed node and its child
		if node == new_tree.seed_node:		 
			continue
		elif node == new_tree.seed_node.child_nodes()[0]:
			node.edge_length = 0
		else:
			node.edge_length = 5

	# print(new_tree.as_string(schema='newick'))

	# print(new_tree.as_ascii_plot(show_internal_node_labels=True))

	# Plot tree
	active_spacers = []
	for array in active_arrays:
		active_spacers += array.spacers

	active_spacers = set(active_spacers)

	ncols = len(active_spacers)

	col_scheme = colour_schemes.choose_col_scheme(ncols)

	spacer_colours = {i: j for i,j in zip(active_spacers, col_scheme)}

	array_dict = {a.id: a for a in all_arrays}

	# Only keep spacers that there is evidence for the existence of in
	# final arrays

	final_array_dict = remove_lost_spacers(array_dict, active_spacers)

	# tree_operations.plot_tree(new_tree, final_array_dict, args.outdir+"test.png", spacer_colours)

	fig_h = args.tree_height
	fig_w = args.tree_width
	# Scale font size according to image dimensions and user-input
	font_scale = min(fig_h,fig_w)*args.font_scale
	dpi = args.dpi
	line_scale = args.branch_weight
	label_text_size = args.font_override_labels
	annot_text_size = args.font_override_annotations
	branch_spacing = 2.3*args.branch_spacing
	brlen_scale=0.5*args.brlen_scale


	tree_operations.plot_tree(
		new_tree, final_array_dict, outdir+"test_temp.png",
		spacer_colours, fig_h=fig_h, fig_w=fig_w, font_scale=font_scale,
		dpi=dpi, line_scale=line_scale, branch_lengths=args.brlen_labels,
		branch_spacing=branch_spacing, brlen_scale=brlen_scale,
		label_text_size=label_text_size, annot_text_size=annot_text_size,
		align_labels=args.align)

if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
