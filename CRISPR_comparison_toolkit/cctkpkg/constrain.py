#! /usr/bin/env python3

import sys
from math import ceil, log
from collections import Counter
import argparse

import dendropy

from . import (
	array_parsimony,
	colour_schemes,
	file_handling,
	tree_operations,
	plotting
	)


def build_parser(parser):

	input_params = parser.add_argument_group('Required inputs')
	input_params.add_argument(
		"-a", "--array-file",
		required=True,
		help="Array_IDs.txt or Array_seqs.txt file."
		)
	input_params.add_argument(
		"-t", "--input-tree",
		required=True,
		help="file containing tree in newick format."
		)
	input_params.add_argument(
		"-g", "--genome-array-file",
		required=True,
		help="file specifying which array ID is present in each genome in "
		"your tree."
		)

	output_params = parser.add_argument_group('Output control', 
		"Set which of the optional outputs you want.")
	output_params.add_argument(
		"-o", "--out-plot",
		required=True, 
		help="Specify output plot file name."
		)
	output_params.add_argument(
		"--output-arrays",
		metavar="",
		required=False,
		help="Specify filename for the details of you final arrays with hypothetical intermediate arrays in the same format as your input array_file. If there are multiple best trees, one file will be created per tree numbered in the order they are described in the stdout output."
		)
	output_params.add_argument(
		"--print-tree",
		action='store_true',  
		help="Print a graphical representation of the tree using ascii characters."
		)

	run_params = parser.add_argument_group('Running parameters', 
		"Control run behaviour.")
	run_params.add_argument(
		"--acquisition",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
		)
	run_params.add_argument(
		"--indel",
		metavar="",
		type=int,
		default=10,
		help="Specify the parsimony cost of an indel event involving one or more spacers. Default: 10"
		)
	run_params.add_argument(
		"--rep-indel",
		metavar="",
		type=int,
		default=50,
		help="Specify the parsimony cost of an indel event involving one or more spacers that is independently acquired in multiple arrays. Default: 50"
		)
	run_params.add_argument(
		"--duplication",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of a duplication event involving one or more spacers. Default: 1"
		)
	run_params.add_argument(
		"--trailer-loss",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of the loss of a spacer from the trailer end of the array. Default: 1"
		)
	run_params.add_argument(
		"--no-ident",
		metavar="",
		type=int,
		default=100,
		help="Specify the parsimony cost of a tree in which an array is predicted to have descended from another array with which it shares no spacers. Default: 100"
		)
	run_params.add_argument(
		"--seed",
		metavar="",
		type=int,
		required=False,
		default=2,
		help="The order of outline and fill colours assigned to spacers is semi-random. Change it by providing a number here to change which colours are assigned to each spacer."
		)

	cs_files = parser.add_argument_group('Colour scheme files', 
		"Set inputs and outputs for optional colour scheme files.")
	cs_files.add_argument(
		"--colour-file",
		metavar=' ',
		required=False, 
		help="Specify file with custom colour list (Optional). \
		Colours must be hex codes. One colour per line with no header \
		line in file. e.g. #fd5925."
		)
	cs_files.add_argument(
		"--colour-scheme-outfile",
		metavar=' ',
		required=False, 
		help="Specify output file to store json format dictionary of \
		the colour schemes used for spacers in this run."
		)
	cs_files.add_argument(
		"--colour-scheme-infile",
		metavar=' ',
		required=False,
		help="Specify input file containing json format dictionary of \
		the colour scheme to be used for spacers in this run. Any \
		spacers not in the input file will be coloured according to \
		the normal process."
		)

	plot_params = parser.add_argument_group('Plotting parameters', 
		"Control elements of the produced plot.")
	plot_params.add_argument(
		"-e",
		dest="emphasize_diffs",
		action='store_true',  
		help="When plotting a representation of the tree with cartooned arrays, emphasize locations where arrays differ from their hypothetical ancestor."
		)
	plot_params.add_argument(
		"-b",
		dest="branch_lengths",
		action='store_true', 
		help="When plotting a representation of the tree Include branch lengths at midpoint on branches. N.B. This does not control inclusion of branch lengths in newick format tree output, which are always included."
		)
	plot_params.add_argument(
		"--brlen-scale",
		type=float,
		required=False,
		metavar='',
		default=1,
		help="Factor to scale branch length."
		)
	plot_params.add_argument(
		"--no-align-cartoons",
		action='store_true',
		default=False,
		help="When plotting a representation of the tree with cartooned arrays, this option controls whether those cartoons are drawn next to the corresponding node in the tree. By default cartoons will be aligned at the trailer end."
		)
	plot_params.add_argument(
		"--no-align-labels",
		action='store_true',
		default=False,
		help="Should node labels be placed at the corresponding internal node or leaf. By default they are all aligned. N.B. For this setting to work, cartoons must also placed at nodes using the -g option."
		)
	plot_params.add_argument(
		"--dpi",
		type=int,
		required=False,
		metavar='',
		default=600,
		help="The desired resolution of the output image."
		)
	plot_params.add_argument(
		"--no-fade-anc",
		action='store_true', 
		help="Do not de-emphasize ancestral array cartoons using transparency. This option helps to make it clear which are the inferred ancestral arrays and which are the arrays being analyzed. However, this option reduces colour contrast between spacers."
		)
	plot_params.add_argument(
		"--plot-width",
		type=float,
		default=3,
		metavar="",
		help="Width of plot in inches. Default = 3"
		)
	plot_params.add_argument(
		"--plot-height",
		type=float,
		default=3,
		metavar="",
		help="Height of plot in inches. Default = 3"
		)
	plot_params.add_argument(
		"--font-override-labels",
		type=float,
		metavar="",
		help="If you want to force a label font size to be used rather than \
			using scaling to determine font size, provide it here"
		)
	plot_params.add_argument(
		"--font-override-annotations",
		type=float,
		metavar="",
		help="If you want to force a specific font size (pts) to be used for \
			annotations such as branch length labels, rather than using \
			scaling to determine font size, provide it here"
		)


	return parser


def main(args):

	if args.no_align_cartoons and not args.no_align_labels:
		args.no_align_labels = True

	event_costs = { 
		"acquisition" : args.acquisition,
		"indel" : args.indel,
		"repeated_indel" : args.rep_indel,
		"duplication": args.duplication,
		"trailer_loss": args.trailer_loss,
		"no_ident": args.no_ident
		}

	array_spacers_dict = file_handling.read_array_file(args.array_file)


	genome_array_dict, outgroup_taxon = file_handling.read_genome_reps_file(
		args.genome_array_file)

	array_dict = {}
	for array in genome_array_dict.values():
		array_dict[array] = array_parsimony.Array(
			array, array_spacers_dict[array], extant=True)
	all_arrays = [array.spacers for array in array_dict.values()]


	all_spacers = []
	for array in array_dict.keys():
		all_spacers += array_spacers_dict[array]
	non_singleton_spacers = [
		spacer for spacer, count in Counter(all_spacers).items() if count >1]

	spacer_cols_dict = colour_schemes.process_colour_args(
		args, non_singleton_spacers)

	tree = dendropy.Tree.get(path=args.input_tree, schema="newick")


	# If tree not rooted using outgroup node then root it
	outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
	tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

	# Now remove outgroup as it isn't used for array analysis
	tree.prune_taxa_with_labels([outgroup_taxon])


	tree = tree_operations.scale_branches(tree, 10)

	labels = [g for g in genome_array_dict.values()]
	node_ids = tree_operations.create_internal_node_ids(
		len(genome_array_dict))


	taxon_namespace = dendropy.TaxonNamespace(labels + node_ids)


	for leaf in tree.leaf_node_iter():
		old_label = leaf.taxon.label
		leaf.taxon = taxon_namespace.get_taxon(genome_array_dict[old_label])


	node_count = 0
	for node in tree.seed_node.postorder_iter():
		if not node.taxon:
			continue

		a = node.taxon.label 
		if len(node.sibling_nodes()) == 0:
			continue

		sister = node.sibling_nodes()[0]
		parent = node.parent_node
		if not sister.taxon:
			continue

		b = sister.taxon.label
		if parent.taxon: # Only need to add the ancestor once.
			continue

		# First make a list of all arrays not in subtree rooted here
		# list of nodes descended from a and b to exclude:
		descendent_array_IDs = [ad.taxon.label for ad in node.postorder_iter() 
			] + [bd.taxon.label for bd in sister.postorder_iter()
			]
		# exclude nodes descended from these. Just include leaves
		extra_subtree_array_IDs = [
			a for a in labels if a not in descendent_array_IDs]
		extra_subtree_array_spacers = [
			array_dict[a].spacers for a in extra_subtree_array_IDs]

		results = array_parsimony.resolve_pairwise_parsimony(
			array_dict[a], 
			array_dict[b], 
			all_arrays, 
			array_dict, 
			node_ids, 
			node_count, 
			tree, 
			event_costs,
			extra_subtree_arrays=extra_subtree_array_spacers
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
	for node in tree.seed_node.postorder_iter():
		if node.level() != 0:
			parent = node.parent_node				
			array_dict[node.taxon.label].reset()
			array_dict[parent.taxon.label].reset()
			array_dict[node.taxon.label] = array_parsimony.count_parsimony_events(
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

	plotting.plot_tree(
		tree=tree,
		array_dict=array_dict,
		filename=args.out_plot,
		spacer_cols_dict=spacer_cols_dict,
		branch_lengths=args.branch_lengths,
		emphasize_diffs=args.emphasize_diffs,
		dpi=args.dpi,
		brlen_scale=args.brlen_scale,
		no_align_cartoons=args.no_align_cartoons,
		no_align_labels=args.no_align_labels,
		no_fade_ancestral=args.no_fade_anc,
		fig_h=args.plot_height,
		fig_w=args.plot_width,
		label_text_size=args.font_override_labels,
		annot_text_size=args.font_override_annotations
		)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="predict array relationships constrained by a tree"
		)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()
	main(sys.argv[1:])
