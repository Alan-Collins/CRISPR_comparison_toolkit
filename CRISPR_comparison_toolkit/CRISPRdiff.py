#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import collections
from itertools import combinations, permutations
from random import sample, randrange, seed
import argparse
import json

from cctk import colour_schemes, file_handling


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

	Given a dict of arrays and the spacers they contain and a list of 
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


def find_indices(lst, element):
	"""Return all indices of list at which an specified element is found.
	
	Given a list and an element found in that list, return all of the 
	indices at which that element is found.
	e.g. for a list ['apple', 'tomatoe', 'apple', 'banana']
	Returns [0,2] for 'apple'

	lst.index() only returns the first instance by default. 
	The second argument provided to index is the position to start
	searching. This approach starts looking again from the index after
	the last found index.


	Args:
	  lst (list): 
	    a list of anything
	  element (any type):
	    An element you expect to find in the list
	
	returns:
	  result (list)
	    A list of indices at which the element was found in the list.
	    Returns an empty list if no indices were found.
	"""
	result = []
	offset = -1
	while True:
		try:
			offset = lst.index(element, offset+1) 
		except ValueError:
			return result
		result.append(offset)


def cmdline_args():
	parser = argparse.ArgumentParser(
		description="Given an array file and a list of arrays to "
		"align, produces a plot showing shared spacers among arrays "
		"and their location within each array. Array file format is "
		"whitespace separated columns where column 1 is the array ID "
		"and columns 3 onwards are spacer IDs, names, or sequences."
		"\ne.g. CRISPRdiff.py -a Array_IDs.txt -o "
		"test.png -c cols.txt -iter 100 155 9 204 73 97",
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument(
		"-a", 
		dest="array_file",
		required=True,
		help="Specify array representatives file."
		)
	parser.add_argument(
		"-o",
		dest="out_file",
		required=True, 
		help="Specify output file."
		)
	parser.add_argument(
		"-c",
		dest="colour_file",
		required=False, 
		help="Specify file with custom colour list (Optional). "
		"Colours must be hex codes. One colour per line with no header "
		"line in file. e.g. #fd5925."
		)
	parser.add_argument(
		"-m",
		dest="colour_scheme_outfile",
		required=False, 
		help="Specify output file to store json format dictionary of "
		"the colour schemes used for spacers in this run."
		)
	parser.add_argument(
		"-s",
		dest="colour_scheme_infile",
		required=False,
		help="Specify input file containing json format dictionary of "
		"the colour scheme to be used for spacers in this run. Any "
		"spacers not in the input file will be coloured according to "
		"the normal process."
		)
	parser.add_argument(
		"-i", "--iter",
		dest="iterations",
		required=False,
		default=10,
		type=int,
		help="(Default = 10) If you are aligning fewer than 9 arrays, "
		"a good order will be found using a repeated shuffling search "
		"method. Set the number of replicates of this search you want "
		"performed. Higher numbers more likely to find the best "
		"possible ordering, but take longer."
		)
	parser.add_argument(
		"-l",
		dest="legend",
		required=False,
		action='store_true',  
		help="Include a legend in the output plot (Highly recommended "
		"to use spacer IDs rather than sequences with this setting). "
		"N.B. Spacer order in the legend is the same as the order of "
		"first instance of spacers working from bottom to top, right "
		"to left along your plotted arrays."
		)
	parser.add_argument(
		"--leader_align",
		dest="leader_align",
		required=False,
		action='store_true',  
		help="Declare that you want the plot to line up all the leader "
		"ends instead of the default of lining up the trailer ends."
		)
	parser.add_argument(
		"--preordered",
		dest="preordered",
		required=False, 
		action='store_true',  
		help="Declare that the array order you provided is the one you "
		"want plotted."
		)
	parser.add_argument(
		"--approxordered",
		dest="approxordered",
		required=False,
		action='store_true',  
		help="Declare that the array order you provided should be "
		"optimized slightly before plotting. Optimization involves "
		"switching order of adjacent arrays in list as long as that "
		"increases the total number of shared spacers among all "
		"neighbours."
		)
	parser.add_argument(
		"--seed",
		dest="seed",
		required=False, 
		type=int,
		default=2,
		help="The order of outline and fill colours assigned to "
		"spacers is semi-random. Change it by providing a number here "
		"to change which colours are assigned to each spacer."
		)
	parser.add_argument(
		"--dpi",
		dest="dpi",
		required=False,
		type=int,
		default=300,
		help="Resolution of output image. Only relevant for bitmap "
		"formats such as PNG. Has no effect on SVG outputs."
		)
	parser.add_argument(
		"--connection_outline",
		dest="connection_outline",
		required=False,
		action='store_true',
		help="Identical spacers in arrays plotted adjacent to one "
		"another are connected by a line with the fill colour of the "
		"spacer by default. If you would like the line connecting "
		"those spacers to have the same outline as the spacers as well "
		"then use this option."
		)
	parser.add_argument(
		"arrays_to_align",
		nargs="+",  
		help="Specify the arrays for which you want to plot an "
		"alignment. **Must come at the end of your command after all "
		"other arguments.**"
		)

	args = parser.parse_args()

	return args


def main(args):

	array_dict = file_handling.read_array_file(args.array_file)

	if len(args.arrays_to_align) == 0:
		arrays_of_interest_dict = array_dict
		array_network = [a for a in array_dict.keys()]
	else:
		arrays_of_interest_dict = {}
		array_network = args.arrays_to_align
		for array in array_network:
			arrays_of_interest_dict[array] = array_dict[array]

	if args.preordered:
		array_order = array_network
	elif args.approxordered:
		array_order, _ = jiggle_list_to_local_max(
			arrays_of_interest_dict, array_network)
	else:
		if len(array_network) < 9:
			array_order = decide_array_order_global_best(
				arrays_of_interest_dict)
		else:
			array_order, local_best_score = decide_array_order_local_best(
				arrays_of_interest_dict, 100, args.iterations)
			print(
				"The score (sum of spacers shared between all neighbouring"
				"arrays) of the best ordering found was: {}".format(
					local_best_score))

	## Find spacers present in more than one array
	all_spacers = list(arrays_of_interest_dict.values())
	occurrences = dict(collections.Counter(
		[item for sublist in all_spacers for item in sublist]
		))

	imp_spacers = []

	for k,v in occurrences.items():
		if v > 1:
			imp_spacers.append(k)

	print("Identified {} spacers present in more than one array.".format(
		len(imp_spacers)))

	## Figure out colour scheme to use

	if args.colour_file:
		with open(args.colour_file, 'r') as fin:
			cf_list = [i.strip() for i in fin.readlines()]
		colours = colour_schemes.choose_col_scheme(
			len(imp_spacers),
			args.seed,
			cf_list)
	else:
		colours = colour_schemes.choose_col_scheme(
			len(imp_spacers),
			args.seed)

	# if provided, read in colour scheme.

	if args.colour_scheme_infile:
		with open(args.colour_scheme_infile, 'r') as fin:
			spacer_colours = json.load(fin)
		if any([s not in spacer_colours.keys() for s in imp_spacers]):
			colour_idx = 0
			for s in imp_spacers:
				if s not in spacer_colours.keys():
					while colours[colour_idx] in spacer_colours.values():
						colour_idx += 1
					spacer_colours[s] = colours[colour_idx]
					colour_idx += 1
	else:
		# build a dictionary with colours assigned to each spacer.
		spacer_colours  = {}

		for i, spacer in enumerate(sorted(imp_spacers)):
			spacer_colours[spacer] = colours[i]


	if args.colour_scheme_outfile:
		with open(args.colour_scheme_outfile, 'w', encoding='utf-8') as fout:
			json.dump(spacer_colours, fout, ensure_ascii=False, indent=4)


	largest_array_size = max(
		[len(x) for x in [array_dict[y] for y in array_network]])

	dim_y = max(len(array_network)*0.5, 2)
	if args.legend:
		ratio_of_heights = ((len(imp_spacers)+1)*0.2)/dim_y
		ncols = int(1+ ratio_of_heights)
		dim_x = max(largest_array_size*0.25 + 1.2*ncols, 3)
	else:
		dim_x = max(largest_array_size*0.25, 2)
	


	plt.rcParams.update({'font.size': 7})
	fig, ax = plt.subplots()

	fig.set_size_inches(dim_x, dim_y)

	### Plot array alignments

	if not args.leader_align:
		longest_array = max(
			[len(arrays_of_interest_dict[i]) for i in array_order])
		pad_dict = {
		array: longest_array-len(arrays_of_interest_dict[array]) for array in array_order}

	else:
		pad_dict = {array:0 for array in array_order}

	# count of which array we are on (i.e. y axis value to plot at)
	arcount = 1 
	for array in array_order:
		# count of which spacer in array we are on 
		# (i.e. which x-axis value to plot at)
		spcount = 0+pad_dict[array] 
		for spacer in arrays_of_interest_dict[array]:
			if spacer in imp_spacers:
				line_width= 0.1
				spcolour = spacer_colours[spacer]
			else: 
				spcolour = ("#000000", "#000000") #black
				line_width = 0.01
			ax.fill_between(
				[spcount+0.60, spcount+1.40], 
				arcount-line_width, 
				arcount+line_width, 
				color = spcolour[0], 
				edgecolor=spcolour[1],
				linewidth=1.5,
				joinstyle='miter',
				zorder=2
				)

			spcount+=1
		arcount+=1

	### Plot lines connecting identical spacers in neighbouring arrays

	arcount = 0
	for i in range(len(array_order)-1):
		array1 = arrays_of_interest_dict[array_order[i]]
		array2 = arrays_of_interest_dict[array_order[i+1]]
		for spacer in list(set(array1).intersection(array2)):
			array1_indices = find_indices(array1, spacer)
			array2_indices = find_indices(array2, spacer)
			for a in array1_indices:
				for b in array2_indices:
					sp_x1 = a+1+pad_dict[array_order[i]]
					sp_y1 = i+1
					sp_x2 = b+1+pad_dict[array_order[i+1]]
					sp_y2 = i+2
					spcolour = spacer_colours[spacer]
					if args.connection_outline:
						ax.plot(
							[sp_x1, sp_x2],
							[sp_y1+0.1, sp_y2-0.1],
							color=spcolour[1],
							linewidth=2.0,
							label=spacer, 
							zorder=0
							)
					ax.plot(
						[sp_x1, sp_x2],
						[sp_y1+0.1, sp_y2-0.1],
						color=spcolour[0],
						linewidth=1.0,
						solid_capstyle="butt",
						label=spacer,
						zorder=1
						)



	axes = plt.gca()
	plt.yticks(range(1, len(arrays_of_interest_dict.keys())+1, 1))
	if args.leader_align:
		plt.tick_params(
			axis='x',
			which='both',
			bottom=False,
			labelbottom=False
			)
	else:
		plt.xticks(
			range(1, max(
					[len(x) for x in arrays_of_interest_dict.values()])+1, 1))
		plt.xlabel("Spacer in array", fontsize=11)

	ax.set_yticklabels(array_order)
	plt.ylabel("Array ID", fontsize=11)

	if args.legend:
		#Approx height of each legend entry is 0.2 inches. Each array
		# is allocated 0.5 inches. Add 1 to number of spacers to
		# account for legend title.
		ratio_of_heights = ((len(imp_spacers)+1)*0.2)/dim_y
		# Split the legend over multiple columns if it would be taller
		# than the y axis. Plus 1 to account for rounding down and a bit of extra security.
		ncols = int(1+ ratio_of_heights)
		# What is the height of the legend now,
		# having split it over columns?
		h_leg = 0.2 + (len(imp_spacers)*0.2)/ncols 
		# Add blank space half the difference in height between legend
		# and y axis to approximately center the legend relative
		# to y axis. Then scale that on 0-1 relative to y axis size.
		vert_pad = ((dim_y - h_leg)/2)/dim_y 
		# Build the legend manually as the automatic legend included
		# duplicates as spacers are plotted once per array they are in.
		colors = list(spacer_colours.values()) 
		lines = [
			plt.fill_between(
				[1, 1],
				1,
				1,
				color=c[0],
				edgecolor=c[1]) 
			for c in colors
			]
		labels = list(spacer_colours.keys())
		plt.legend(
			lines,
			labels,
			ncol=ncols,
			loc=(1.01,vert_pad),
			title="Spacer IDs"
			) 


	fig.tight_layout()

	plt.savefig(args.out_file, dpi=args.dpi)
	plt.close(fig)


if __name__ == '__main__':
	args = cmdline_args()

	main(args)
