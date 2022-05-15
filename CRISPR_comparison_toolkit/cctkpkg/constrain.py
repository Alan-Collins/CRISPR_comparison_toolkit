#! /usr/bin/env python3

import sys
from math import ceil, log
from collections import Counter
import argparse
from copy import copy, deepcopy
import textwrap

import dendropy

description = """
usage: cctk constrain [-h] -a -t -g -o [--output-arrays] [--print-tree] [-u] \
[--acquisition] [--deletion] [--insertion] [--rep-indel] [--duplication] \
[--trailer-loss] [--no-ident] [--seed] [--colour-file] [--colour-scheme-outfile] \
[--colour-scheme-infile] [-e] [-b] [--brlen-scale] [--no-align-cartoons] \
[--no-align-labels] [--dpi] [--no-fade-anc] [--plot-width] [--plot-height] \
[--font-override-labels] [--font-override-annotations]

optional arguments:
  -h, --help        show this help message and exit

inputs:
  -a, --array-file  Array_IDs.txt or Array_seqs.txt file
  -t, --input-tree  file containing tree in newick format
  -g, --genome-array-file	
                    file corresponding array ID and genome ID

output control:
  set which of the optional outputs you want.

  -o, --out-plot    output plot file name
  --output-arrays   file to store analyzed arrays and hypothetical ancestors
  --print-tree      print an ascii symbol representation of the tree
  -u, --unrooted    input tree is unrooted

running parameters:
  control run behaviour.

  --acquisition     parsimony cost of a spacer acquisition event. Default: 1
  --deletion        parsimony cost of a deletion event. Default: 10
  --insertion       parsimony cost of an insertion event. Default: 30
  --rep-indel       parsimony cost independently acquiring spacers. Default: 50
  --duplication     parsimony cost of a duplication event. Default: 1
  --trailer-loss    parsimony cost of trailer spacer loss. Default: 1
  --no-ident        parsimony cost of an array having no identity with its \
ancestor. Default: 100
  --seed            set seed for random processes

colour scheme files:
  set inputs and outputs for optional colour scheme files.

  --colour-file     file with custom colour list
  --colour-scheme-outfile
                    output file to store json format colour schemes
  --colour-scheme-infile
                    input file json format colour scheme

plotting parameters:
  control elements of the produced plot.

  -b                include branch lengths in tree plot
  --brlen-scale     factor to scale branch length.
  --no-emphasize-diffs
                    don't emphasize events in each array since its ancestor
  --no-align-cartoons
                    draw array cartoons next to leaf node
  --no-align-labels
                    draw leaf labels next to leaf node
  --replace-brlens  replace input tree branch lengths with array parsimony costs
  --no-fade-anc     do not apply transparency to ancestral array depiction
  --plot-width      width of plot in inches. Default = 3
  --plot-height     height of plot in inches. Default = 3
  --font-override-labels
                    set label font size in pts
  --font-override-annotations
                    set annotation font size in pts
"""

from . import (
	array_parsimony,
	colour_schemes,
	file_handling,
	tree_operations,
	plotting
	)

def find_best_array_no_id(array_a, array_b, array_dict):
	best_array = None
	best_score = 0
	for array in array_dict.values():
		if array.id in [array_a.id, array_b.id]:
			continue
		a_score = len(set(array.spacers).intersection(set(array_a.spacers)))
		b_score = len(set(array.spacers).intersection(set(array_b.spacers)))
		
		if a_score > best_score:
			best_array = array_a
			best_score = copy(a_score)

		if b_score > best_score:
			best_array = array_b
			best_score = copy(b_score)
	return best_array


def build_parser(parser):

	input_params = parser.add_argument_group('Iinputs')
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
	output_params.add_argument(
		"-u", "--unrooted",
		action='store_true',  
		help="Specify that the tree is unrooted. Will be rooted using the "
		"first outgroup taxon in the input genome array file."
		)
	run_params.add_argument(
		"--acquisition",
		metavar="",
		type=int,
		default=1,
		help="Specify the parsimony cost of a spacer acquisition event. Default: 1"
		)
	run_params.add_argument(
		"--deletion",
		metavar="",
		type=int,
		default=10,
		help="Specify the parsimony cost of a deletion event involving one or more spacers. Default: 10"
		)
	run_params.add_argument(
		"--insertion",
		metavar="",
		type=int,
		default=30,
		help="Specify the parsimony cost of an insertion event involving one or more spacers. Default: 30"
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
		"--no-emphasize-diffs",
		dest="emphasize_diffs",
		action='store_false',  
		help="When plotting a representation of the tree with cartooned arrays, don't emphasize locations where arrays differ from their hypothetical ancestor."
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
		"--replace-brlens",
		action='store_true',
		default=False,
		help="Should branch lengths in the tree be replaced by array parsimony costs?"
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
		"acquisition": args.acquisition,
		"deletion": args.deletion,
		"insertion": args.insertion,
		"repeated_indel": args.rep_indel,
		"duplication": args.duplication,
		"trailer_loss": args.trailer_loss,
		"no_ident": args.no_ident
		}

	array_spacers_dict = file_handling.read_array_file(args.array_file)
	array_spacers_dict["none"] = []


	genome_array_dict, outgroup_taxons = file_handling.read_genome_reps_file(
		args.genome_array_file)

	array_dict = {}
	for array_list in genome_array_dict.values():
		for array in array_list:
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

	tree = dendropy.Tree.get(
		path=args.input_tree,
		schema="newick"
		)

	if args.unrooted:
		# If tree not rooted, use outgroup node to root it
		outgroup_node = tree.find_node_with_taxon_label(outgroup_taxons[0])
		tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

	# Now remove outgroup as it isn't used for array analysis
	tree.prune_taxa_with_labels(outgroup_taxons)

	(
		tree,
		genome_array_dict
	) = tree_operations.process_duplicate_arrays_contrain(
		tree,
		genome_array_dict
	)

	tree.reroot_at_node(tree.seed_node, update_bipartitions=True)

	tree = tree_operations.scale_branches(tree, 10)


	labels = [g for g in genome_array_dict.keys()]

	# Add new labels to array_dict
	for label in labels:
		array_id = label.split(".")[-1]
		array_dict[label] = deepcopy(array_dict[array_id])
		# update ID
		array_dict[label].id = label

	# Remove unused arrays
	to_del = []
	for array in array_dict:
		if array not in labels:
			to_del.append(array)
	for a in to_del:
		del array_dict[a]

	node_ids = tree_operations.create_internal_node_ids(
		len(genome_array_dict), "Anc ")

	tns = tree.taxon_namespace
	for new_id in node_ids:
		tax = dendropy.Taxon(new_id)
		tns.add_taxon(tax)


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
		if results == "No_ID":
			# Choose the one that shares the most spacers with another array
			# and make a copy of that the ancestor
			array_a = array_dict[a]
			array_b = array_dict[b]
			best_array = find_best_array_no_id(array_a, array_b, array_dict)
			ancestor = array_parsimony.Array(
				node_ids[node_count],
				spacers=copy(best_array.spacers),
				extant=False)

		else:
			array_dict[a], array_dict[b], ancestor = results
		
		node_count+=1

		for n, n_id in [(node, a), (sister, b)]:
			for k,v in event_costs.items(): # Get weighted distance based on each event's cost.
				d = array_dict[n_id].events[k]*v
				array_dict[n_id].distance += d
				

		array_dict[ancestor.id] = ancestor
		node.parent_node.taxon = tns.get_taxon(ancestor.id)
		if args.replace_brlens:
			node.parent_node.edge_length = 0 # Start the ancestor with 0 branch length

	# Repeat iteration now that tree is built to add repeat indels.
	# and calculate tree score
	score = 0
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
				score += d
				if args.replace_brlens:
					node.edge_length = array_dict[node.taxon.label].distance

	sys.stderr.write("Total tree score is {}\n\n".format(score))

	print(tree.as_string(
		"newick",
		suppress_internal_node_labels=True,
		suppress_rooting=True
		))

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
