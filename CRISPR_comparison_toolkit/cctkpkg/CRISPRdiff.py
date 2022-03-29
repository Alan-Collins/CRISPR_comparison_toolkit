#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
from collections import defaultdict
from itertools import combinations, permutations
from random import sample, randrange, seed
import argparse

from . import (colour_schemes,
	file_handling,
	sequence_operations,
	plotting)

description = """
usage: cctk CRISPRdiff [-h] -a -o [--iterations] [--preordered] \
[--approx-ordered] [--seed] [--colour-file] [--colour-scheme-outfile] \
[--colour-scheme-infile] [--line-width] [--dpi] [--connection-outline] \
[--plot-width] [--plot-height] [--font-size] [arrays_to_align]

positional arguments:
  arrays_to_align    IDs of the arrays you want to analyse. Default: all

optional arguments:
  -h, --help         show this help message and exit
  -a, --array-file   Array_IDs.txt or Array_seqs.txt
  -o, --out-file     output plot file name

running parameters:
  control run behaviour

  --iterations      number of iterations of order determination. Default = 10
  --preordered      array order you provided is the one you want plotted
  --approx-ordered  array order you provided should be optimized slightly
  --seed            set seed for random processes

colour scheme files:
  set inputs and outputs for optional colour scheme files

  --colour-file     file with custom colour list
  --colour-scheme-outfile
                    output file to store json format colour schemes
  --colour-scheme-infile
                    input file json format colour scheme

plotting parameters:
  control elements of the produced plot

  --line-width      scale factor for lines between identical spacers
  --dpi             resolution of output image
  --connection-outline 
                    add outline colour to lines connecting identical spacers
  --plot-width      width of plot in inches. Default = 3
  --plot-height     height of plot in inches. Default = 3
  --font-size       font size. Defualt 10pt.
"""

def get_list_score(arrays_dict, arrays_order):
	"""Returns number of shared spacers between all neighbouring arrays.
	
	Given an ordered list of arrays and a dict of array spacer content.
	For each set of neighbouring arrays in the list, finds the number of
	spacers shared among those two arrays. Sums number of shared spacers
	between all neighbouring arrays in list.
	e.g. 
	Array1 has spacers 1,2,3,4,5
	Array2 has spacers 9,2,3,5,6,7,8,
	Array3 has spacers 1,6,7,8
	Array1 and Array2 share spacers 2,3,5 = 3 spacers
	Array2 and Array3 share spacers 6,7,8 = 3 spacers
	Array1 and Array3 are not neighbours so their shared spacers are not
	counted.
	The score for this order would therefore be 3 + 3 = 6
	
	Args:
	  arrays_dict (dict):
	    dict object containing all arrays in your dataset as keys and a
	    list of the spacers in those arrays as values. 
	    Format: {
	      Array1: [Spacer1, Spacer2, ...], 
	      Array2: [Spacer45, Spacer2, ...], 
	      ...
	      }
	  arrays_order (list):
	    ordered list of arrays to be scored. 
	    Format: [Array1, Array2, Array3, Array4]
	
	returns:
	  Score (int) 
	  The total number of shared spacers among neighbouring arrays in
	  the provided order list of arrays.)
	"""
	overlap_dict = {}
	for i in combinations(list(arrays_dict.keys()), 2):
		overlap_dict[i]=len(
			list(set(arrays_dict[i[0]]).intersection(arrays_dict[i[1]])))
	score = 0
	for i in range(len(arrays_order)-1):
		if (arrays_order[i], arrays_order[i+1]) in overlap_dict.keys():
			score+= overlap_dict[(arrays_order[i], arrays_order[i+1])]
		elif (arrays_order[i+1], arrays_order[i]) in overlap_dict.keys():
			score+= overlap_dict[(arrays_order[i+1], arrays_order[i])]
		else:
			print("can't find combination %s, %s" %(
				arrays_order[i+1], arrays_order[i]))
	return score


def decide_array_order_global_best(arrays_dict):
	"""Find an array order to maximise shared spacers in plot.
	
	Given a dict of arrays and their spacers, where the arrays in that 
	dict are to be ordered such that the neighbours share the most 
	possible spacers, generates a list of all possible orders of spacers
	and calls get_list_score to find the best possible array order.

	Args:
	  arrays_dict (dict): 
	    dict object containing all arrays in your dataset as keys and a
	    list of the spacers in those arrays as values. 
	    Format: {
	      Array1: [Spacer1, Spacer2, ...], 
	      Array2: [Spacer45, Spacer2, ...], 
	      ...
	      }
	
	returns:
	  best_order (list)
	    The order of arrays that results in the highest total number of
	    shared spacers among neighbouring arrays in the list.
	"""
	elements = list(arrays_dict.keys())
	possible_orders = list(permutations(elements, len(elements)))
	score = 0
	best_order = list(possible_orders[0])
	for i in possible_orders:
		new_score = get_list_score(arrays_dict, i)
		if new_score > score:
			score = int(new_score)
			best_order = list(i)
	return best_order


def jiggle_list_to_local_max(arrays_dict, order):
	"""Perform slight changes to array order to improve score.

	Given a dict of arrays and the spacers they contain, and a list of 
	arrays, first calculates the score for the order of arrays using 
	get_list_score function. Then swaps the order to the first and 
	second array in the list and calculates score again. If score 
	improved, keeps the new order. Otherwise discards new order.
	Proceeds through the list, swapping each array with the preceding 
	and then the following array in the list, keeping the new order 
	whenever the score improves. Repeats this process until it goes 
	through the list of arrays from first to last array without any of 
	the rearrangements increasing the score.

	Args:
	  arrays_dict (dict):
	    dict object containing all arrays in your dataset as keys and a
	    list of the spacers in those arrays as values. 
	    Format: {
	    Array1: [Spacer1, Spacer2, ...], 
	    Array2: [Spacer45, Spacer2, ...], 
	    ...
	    }
	  order (list): 
	    ordered list of arrays to be scored. 
	    Format: [Array1, Array2, Array3, Array4]
	
	returns: 
	  tuple of
	    (
	    best_score (int) The score of the best order found,
	    best_order (list) the best order found
	    )

	"""
	repeat = True
	best_order = list(order)
	best_score = get_list_score(arrays_dict, order)
	while repeat:
		repeat = False
		for i in range(len(order)):
			try:
				# switch list element with the element before it 
				order[i], order[i-1] = order[i-1], order[i]
				new_score = get_list_score(arrays_dict, order)
				if new_score > best_score:
					repeat = True
					best_score = int(new_score)
					best_order = list(order)
			except:
				pass
			try:
				# switch list element with the element after it
				order[i], order[i+1] = order[i+1], order[i] 
				new_score = get_list_score(arrays_dict, order)
				if new_score > best_score:
					repeat = True
					best_score = int(new_score)
					best_order = list(order)
			except:
				pass
	return best_order, best_score


def shuffle_random_arrays(order):
	"""Shuffle order of list of arrays.

	Given a list of arrays, picks two random indices in the list and
	swaps those elements of the list.

	Args:
	  order (list):
	    ordered list of arrays to be scored.
	    Format: [Array1, Array2, Array3, Array4]
	
	returns:
	  order (list) new order of arrays with two random indices swapped
	"""
	ran = len(order)
	a = randrange(ran)
	b = randrange(ran)
	while b == a:
		b = randrange(ran)
	order[a], order[b] = order[b], order[a]

	return order


def decide_array_order_local_best(arrays_dict, reps, nits):
	"""Find good order for arrays to maximise shared spacers.

	Given a dict of arrays and the spacers they contain and a list of 
	arrays, creates a random initial order of arrays in a list. 
	Then calls shuffle_random_arrays function to switch the position of
	two random arrays in the list. If this shuffling improves the score 
	calculated by get_list_score, keep the new order. If not discard the
	new order and proceed with the old order. Repeat this shuffling step
	until it has failed to improve the order how ever many times is
	defined by the "reps" argument.	Repeat this whole process with a new
	initial random order however many times is defined by the "nits"
	argument (n iterations). Returns the order and score that were
	highest among all of the iterations

	Args:
	  arrays_dict (dict):
	    dict object containing all arrays in your dataset as keys and a
	    list of the spacers in those arrays as values. 
	    Format: {
	      Array1: [Spacer1, Spacer2, ...],
	      Array2: [Spacer45, Spacer2, ...],
	      ...
	      }
	  reps (int):
	    Number of consecutive shuffling steps that should result in no
	    increase in score before the function stops optimizing this
	    order of arrays
	  nits (int):
	    Number of times to repeat the whole process with a different
	    initial random array order
	  returns:
	    tuple of:
	    (
	      overall_best_order (list)
	        Highest scoring order of arrays found among all the 
	        iterations,
		  overall_best_score (list)
		    Highest score found among all the iterations,
		  )
	"""
	elements = list(arrays_dict.keys())
	overall_best_score = 0
	for i in range(nits):
		best_order = ''
		best_score = 0
		rep = 0	
		order = sample(elements, len(elements))
		while rep < reps:
			rep+=1
			order = shuffle_random_arrays(order) 
			score = get_list_score(arrays_dict, order)
			if score > best_score:
				rep = 0
				best_score = int(score)
				best_order = list(order)
		if best_score > overall_best_score:
			overall_best_score = int(best_score)
			overall_best_order = list(best_order)
	return overall_best_order, overall_best_score


def build_parser(parser):
	parser.add_argument(
		"-a", "--array-file",
		required=True,
		help="Specify array representatives file."
		)
	parser.add_argument(
		"-o", "--out-file",
		required=True, 
		help="Specify output file."
		)

	run_params = parser.add_argument_group('Running parameters', 
		"Control run behaviour.")
	run_params.add_argument(
		"--iterations",
		required=False,
		default=10,
		metavar=' ',
		type=int,
		help="Number of attempts to find the best order of arrays for \
		plotting. Higher number can improve the array order. Default = 10"
		)
	run_params.add_argument(
		"--preordered",
		required=False, 
		action='store_true',  
		help="Declare that the array order you provided is the one you \
		want plotted."
		)
	run_params.add_argument(
		"--approx-ordered",
		required=False,
		action='store_true',  
		help="Declare that the array order you provided should be \
		optimized slightly before plotting. Optimization involves \
		switching order of adjacent arrays in list as long as that \
		increases the total number of shared spacers among all \
		neighbours."
		)
	run_params.add_argument(
		"--seed",
		metavar=' ',
		required=False, 
		type=int,
		default=2,
		help="The order of outline and fill colours assigned to \
		spacers is semi-random. Change it by providing a number here \
		to change which colours are assigned to each spacer."
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
		"--line-width",
		metavar=' ',
		required=False,
		default=1,
		type=float,  
		help="Control the width of lines connecting shared spacers. Line width \
		will be multiplied by the given number. Default = 1.0"
		)
	plot_params.add_argument(
		"--dpi",
		metavar=' ',
		required=False,
		type=int,
		default=300,
		help="Resolution of output image. Only relevant for bitmap \
		formats such as PNG. Has no effect on SVG outputs."
		)
	plot_params.add_argument(
		"--connection-outline",
		required=False,
		action='store_true',
		help="Add outline colour to lines connecting identical spacers."
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
		"--font-size",
		type=float,
		metavar="",
		default=10,
		help="Set font size. Defualt 10pt."
		)

	parser.add_argument(
		"arrays_to_align",
		nargs="*",  
		help="Specify the arrays for which you want to plot an \
		alignment. **Must come at the end of your command after all \
		other arguments.**"
		)

	return parser


def main(args):

	array_dict = file_handling.read_array_file(args.array_file)

	if len(args.arrays_to_align) == 0:
		arrays_of_interest_dict = array_dict
		arrays_of_interest = [a for a in array_dict.keys()]
	else:
		arrays_of_interest_dict = {}
		arrays_of_interest = args.arrays_to_align
		for array in arrays_of_interest:
			arrays_of_interest_dict[array] = array_dict[array]

	if args.preordered:
		# Reverse order so order of list plots top to bottom
		array_order = arrays_of_interest[::-1]
	elif args.approx_ordered:
		array_order, _ = jiggle_list_to_local_max(
			arrays_of_interest_dict, arrays_of_interest[::-1])
	else:
		if len(arrays_of_interest) < 9:
			array_order = decide_array_order_global_best(
				arrays_of_interest_dict)
		else:
			array_order, local_best_score = decide_array_order_local_best(
				arrays_of_interest_dict, 100, args.iterations)
			print(
				"The score (sum of spacers shared between all neighbouring "
				"arrays) of the best ordering found was: {}".format(
					local_best_score))

	## Find spacers present in more than one array
	all_spacers = list(arrays_of_interest_dict.values())
	occurrences = defaultdict(int)
	for array in all_spacers:
		for spacer in set(array):
			occurrences[spacer]+=1

	imp_spacers = []

	for k,v in occurrences.items():
		if v > 1:
			imp_spacers.append(k)

	sys.stderr.write(
		"Identified {} spacers present in more than one array.\n".format(
			len(imp_spacers)))

	# Process colour file related args
	spacer_colours = colour_schemes.process_colour_args(
		args, imp_spacers)

	plotting.plot_diffplot(arrays_of_interest_dict, array_order, imp_spacers,
		spacer_colours, text_size=args.font_size,
		plot_width=args.plot_width, plot_height=args.plot_height,
		dpi=args.dpi, outfile=args.out_file, line_width=args.line_width,
		 connection_outline=args.connection_outline)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Given an array file and a list of arrays to "
		"align, produces a plot showing shared spacers among arrays "
		"and their location within each array. Array file format is "
		"whitespace separated columns where column 1 is the array ID "
		"and columns 3 onwards are spacer IDs, names, or sequences."
		"\ne.g. CRISPRdiff.py -a Array_IDs.txt -o "
		"test.png -c cols.txt -iter 100 155 9 204 73 97",
		)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)
