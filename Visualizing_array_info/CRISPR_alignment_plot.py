#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import collections
from itertools import combinations
from itertools import permutations
from random import sample
from random import randrange
import argparse

# hex values from this website http://phrogz.net/css/distinct-colors.html

Cols_hex_27 = ["#fd5925", "#dbc58e", "#008d40", "#304865", "#934270", "#f7b8a2", "#907500", "#45deb2", "#1f4195", "#d67381", "#8e7166", "#afb200", "#005746", "#a598ff", "#8f0f1b", "#b96000", "#667f42", "#00c7ce", "#9650f0", "#614017", "#59c300", "#1a8298", "#b5a6bd", "#ea9b00", "#bbcbb3", "#00b0ff", "#cd6ec6"]

#hex values from https://mokole.com/palette.html

Cols_hex_40 = ["#696969","#556b2f","#a0522d","#800000","#006400","#808000","#483d8b","#3cb371","#008080","#bdb76b","#4682b4","#000080","#9acd32","#32cd32","#daa520","#7f007f","#ff4500","#00ced1","#ff8c00","#c71585","#0000cd","#00ff00","#9400d3","#dc143c","#00bfff","#f4a460","#adff2f","#da70d6","#ff00ff","#1e90ff","#db7093","#fa8072","#ffff54","#dda0dd","#7b68ee","#afeeee","#98fb98","#7fffd4","#ffe4c4","#ffc0cb"]

# hex values from http://medialab.github.io/iwanthue/

Cols_hex_50 = ["#66a939","#e1893e","#2b1d31","#9f7be3","#5787cd","#5ddad5","#80447d","#a28687","#d93e6c","#666030","#71d39e","#de5527","#b0d7cc","#3c1c6a","#5c8e77","#db3c3d","#78e049","#437635","#46281a","#565ab9","#9b3a4d","#d97dd2","#804e38","#3c5e70","#393569","#e2b8d0","#bbd094","#59a3bc","#644ddc","#dcb23d","#da7f76","#ae93ce","#9b6923","#949535","#5fdd83","#b49261","#df42a7","#9833b6","#cf4ce1","#e5c5a0","#5e1d35","#95391e","#d3e037","#a8bee3","#9a2b7a","#725972","#4f2096","#2b412a","#cde175","#d3799e"]

Cols_hex_50_2 = ["#00e9b9","#5b00b8","#00da50","#ed00c5","#39a600","#7900a5","#e5ff67","#000b71","#eac900","#015ce5","#b6a700","#a266ff","#01da87","#f6009b","#9effa6","#8b0081","#ffd362","#6d81ff","#fd4b00","#1ab5ff","#ff8025","#017bd5","#ba4200","#00cdb5","#ff5749","#02a69f","#9d004c","#bdffea","#0a003f","#feffae","#e091ff","#a98100","#d2bbff","#2f4300","#ffa9e6","#766800","#a8ddff","#6a2400","#dcfeff","#230000","#ffdbe7","#202000","#ffd7b5","#00334a","#ff929f","#00604a","#5c001f","#00a3bb","#463a00","#007a9e"]

Cols_hex_100 = ["#82cd00","#2f3ce3","#81fd3c","#8058ff","#77ff6a","#023acd","#efe100","#570099","#59c000","#eb61ff","#00c64a","#ae00b4","#c5ff6d","#28006d","#fffa63","#0141ab","#bbbd00","#c181ff","#88af00","#7c0081","#49ffaf","#fa004e","#01c360","#ff64b6","#009d42","#ff99ff","#106b00","#ff83d7","#9affa1","#3f004f","#d8ff7d","#001c60","#ffa329","#0061bf","#f37b00","#8793ff","#999400","#dea9ff","#547100","#9d0062","#00c085","#ba0056","#00c5a2","#ff5d27","#02d9fd","#9c0013","#00dcdc","#b00036","#01c0b5","#ff5589","#00742e","#ff6a82","#00caff","#b24700","#94f9ff","#990021","#c6ffaa","#340032","#d1ffc4","#0f002c","#ffee9d","#001d49","#ffc96b","#002647","#ffa762","#006bac","#aa7100","#0081b3","#ff7163","#02a3c6","#954400","#bdfff6","#2e0018","#c3ffde","#600028","#fff5c9","#001521","#fff4f0","#1b0d00","#eaccff","#273f00","#bbb9ff","#977b00","#b1daff","#792600","#00857d","#ff977f","#003c27","#ffd2e4","#371100","#ffcea5","#003b4b","#ffb8af","#006158","#4e0900","#006f82","#7a4900","#015466","#572c00","#3b3000"]



def get_list_score(arrays_dict, arrays_order):
	"""
	Given an ordered list of arrays and a dict of array spacer content.
	For each set of neighbouring arrays in the list, finds the number of spacers shared among those two arrays.
	Sums number of shared spacers between all neighbouring arrays in list.
	e.g. 
	Array1 has spacers 1,2,3,4,5
	Array2 has spacers 9,2,3,5,6,7,8,
	Array3 has spacers 1,6,7,8
	Array1 and Array2 share spacers 2,3,5 = 3 spacers
	Array2 and Array3 share spacers 6,7,8 = 3 spacers
	Array1 and Array3 are not neighbours so their shared spacers are not counted.
	The score for this order would therefore be 3 + 3 = 6

	:param arrays_dict (dict): dict object containing all arrays in your dataset as keys and a list of the spacers in those arrays as values. Format: {Array1: [Spacer1, Spacer2, ...], Array2: [Spacer45, Spacer2, ...], ...}
	:param arrays_order (list): ordered list of arrays to be scored. Format: [Array1, Array2, Array3, Array4]
	:returns: Score (int) (The total number of shared spacers among neighbouring arrays in the provided order list of arrays.)
	"""
	overlap_dict = {}
	for i in combinations(list(arrays_dict.keys()), 2):
		overlap_dict[i]=len(list(set(arrays_dict[i[0]]).intersection(arrays_dict[i[1]])))
	score = 0
	for i in range(len(arrays_order)-1):
		if (arrays_order[i], arrays_order[i+1]) in overlap_dict.keys():
			score+= overlap_dict[(arrays_order[i], arrays_order[i+1])]
		elif (arrays_order[i+1], arrays_order[i]) in overlap_dict.keys():
			score+= overlap_dict[(arrays_order[i+1], arrays_order[i])]
		else:
			print("can't find combination %s, %s" %(arrays_order[i+1], arrays_order[i]))
	return score


def decide_array_order_global_best(arrays_dict):
	"""
	Given a dict of arrays and their spacers, where the arrays in that dict are to be ordered such that the neighbours share the most possible spacers, 
	generates a list of all possible orders of spacers and calls get_list_score to find the best possible array order.

	:param arrays_dict (dict): dict object containing all arrays in your dataset as keys and a list of the spacers in those arrays as values. Format: {Array1: [Spacer1, Spacer2, ...], Array2: [Spacer45, Spacer2, ...], ...}
	:returns: best_order (list) (The order of arrays that results in the highest total number of shared spacers among neighbouring arrays in the list.)
	"""
	elements = list(arrays_dict.keys())
	possible_orders = list(permutations(elements, len(elements)))
	score = 0
	for i in possible_orders:
		new_score = get_list_score(arrays_dict, i)
		if new_score > score:
			score = int(new_score)
			best_order = list(i)
	return best_order


def jiggle_list_to_local_max(arrays_dict, order):
	"""
	Given a dict of arrays and the spacers they contain and a list of arrays, first calculates the score for the order of arrays using get_list_score function.
	Then swaps the order to the first and second array in the list and calculates score again. If score improved, keeps the new order. Otherwise discards new order.
	Proceeds through the list, swapping each array with the preceding and then the following array in the list, keeping the new order whenever the score improves.
	Repeats this process until it goes through the list of arrays from first to last array without any of the rearrangements increasing the score.

	:param arrays_dict (dict): dict object containing all arrays in your dataset as keys and a list of the spacers in those arrays as values. Format: {Array1: [Spacer1, Spacer2, ...], Array2: [Spacer45, Spacer2, ...], ...}
	:param order (list): ordered list of arrays to be scored. Format: [Array1, Array2, Array3, Array4]
	:returns: 	best_score (int) (The score of the best order found)
				best_order (list) (the best order found)

	"""
	repeat = True
	best_order = list(order)
	best_score = get_list_score(arrays_dict, order)
	while repeat:
		repeat = False
		for i in range(len(order)):
			try:
				order[i], order[i-1] = order[i-1], order[i] # switch list element with the element before it 
				new_score = get_list_score(arrays_dict, order)
				if new_score > best_score:
					repeat = True
					best_score = int(new_score)
					best_order = list(order)
			except:
				pass
			try:
				order[i], order[i+1] = order[i+1], order[i] # switch list element with the element after it
				new_score = get_list_score(arrays_dict, order)
				if new_score > best_score:
					repeat = True
					best_score = int(new_score)
					best_order = list(order)
			except:
				pass
	return best_order, best_score


def shuffle_random_arrays(order):
	"""
	Given a list of arrays, picks two random indices in the list and swaps those elements of the list.

	:param order (list): ordered list of arrays to be scored. Format: [Array1, Array2, Array3, Array4]
	:returns: order (list) (new order of arrays with two random indices swapped)
	"""
	ran = len(order)
	a = randrange(ran)
	b = randrange(ran)
	while b == a:
		b = randrange(ran)
	order[a], order[b] = order[b], order[a]

	return order


def decide_array_order_local_best(arrays_dict, reps, nits):
	"""
	Given a dict of arrays and the spacers they contain and a list of arrays, creates a random initial order of arrays in a list. Then calls shuffle_random_arrays function to switch the position of two random arrays in the list.
	If this shuffling improves the score calculated by get_list_score, keep the new order. If not discard the new order and proceed with the old order.
	Repeat this shuffling step until it has failed to improve the order how ever many times is defined by the "reps" argument.
	Repeat this whole process with a new initial random order however many times is defined by the "nits" argument (n iterations).
	Returns the order and score that were highest among all of the iterations

	:param arrays_dict (dict): dict object containing all arrays in your dataset as keys and a list of the spacers in those arrays as values. Format: {Array1: [Spacer1, Spacer2, ...], Array2: [Spacer45, Spacer2, ...], ...}
	:param reps (int): Number of consecutive shuffling steps that should result in no increase in score before the function stops optimizing this order of arrays
	:param nits (int): Number of times to repeat the whole process with a different initial random array order
	:returns: 	overall_best_order (list) (Highest scoring order of arrays found among all the iterations)
				overall_best_score (list) (Highest score found among all the iterations)
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



# def decide_array_order_local_best(arrays_dict, reps): Old function that started with random order shuffled neighbours and stored best alignment.
# 	elements = list(arrays_dict.keys())
# 	rep = 0
# 	best_order = ''
# 	best_score = 0
# 	while rep < reps:
# 		rep+=1
# 		order = sample(elements, len(elements))
# 		order, score = jiggle_list_to_local_max(arrays_dict, order) #shuffle(elements))
# 		if score > best_score:
# 			best_score = int(score)
# 			best_order = list(order)
# 	return best_order


def find_indices(lst, element):
	"""
	Given a list and an element found in that list, return all of the indices at which that element is found.
	e.g. for a list ['apple', 'tomatoe', 'apple', 'banana']
	Returns [0,2] for 'apple'

	:param lst (list): a list of anything
	:param element (any type): An element you expect to find in the list
	:returns: result (list) (A list of indices at which the element was found in the list. Returns an empty list if no indices were found.)
	"""
	result = []
	offset = -1
	while True:
		try:
			offset = lst.index(element, offset+1) # lst.index() only returns the first instance by default. The second argument provided to index is the position to start searching. This approach starts looking again from the index after the last found index.
		except ValueError:
			return result
		result.append(offset)


parser = argparse.ArgumentParser(
	description="Given an array file and a list of arrays to align, produces a plot showing shared spacers among arrays and their location within each array. Array file format is whitespace separated columns where column 1 is the array ID and columns 3 onwards are spacer IDs, names, or sequences.\ne.g. python3 CRISPR_alignment_plot.py -a Array_IDs.txt -o test.png -c cols.txt -iter 100 155 9 204 73 97")
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"-o", dest="out_file", required = True, 
	help="Specify output file."
	)
parser.add_argument(
	"-c", dest="colour_file", required = False, 
	help="Specify file with custom colour list (Optional). Colours must be hex codes. One colour per line with no header line in file. e.g. #fd5925 or rbg(10,100,5)."
	)
parser.add_argument(
	"-iter", dest="iterations", nargs="?", default = 10, type=int,
	help="(Default = 10) If you are aligning fewer than 9 arrays, a good order will be found using a repeated shuffling search method. Set the number of replicates of this search you want performed. Higher numbers more likely to find the best possible ordering, but take longer."
	)
parser.add_argument(
	"-l", dest="legend", action='store_true',  
	help="Include a legend in the output plot (Highly recommended to use spacer IDs rather than sequences with this setting). N.B. Spacer order in the legend is the same as the order of first instance of spacers working from bottom to top, right to left along your plotted arrays."
	)
parser.add_argument(
	"-preordered", dest="preordered", action='store_true',  
	help="Declare that the array order you provided is the one you want plotted."
	)
parser.add_argument(
	"-approxordered", dest="approxordered", action='store_true',  
	help="Declare that the array order you provided should be optimized slightly before plotting. Optimization involves switching order of adjacent arrays in list as long as that increases the total number of shared spacers among all neighbours."
	)
parser.add_argument(
	"arrays_to_align", nargs="+",  
	help="Specify the arrays for which you want to plot an alignment. **Must come at the end of your command after all other arguments.**"
	)

args = parser.parse_args(sys.argv[1:])

infile = args.array_file
outfile = args.out_file
array_network = args.arrays_to_align



array_dict = {}
with open(infile, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

arrays_of_interest_dict = {}
for array in array_network:
	arrays_of_interest_dict[array] = array_dict[array]

if args.preordered:
	array_order = array_network
elif args.approxordered:
	array_order, score_dump = jiggle_list_to_local_max(arrays_of_interest_dict, array_network)
else:
	if len(array_network) < 9:
		array_order = decide_array_order_global_best(arrays_of_interest_dict)
	else:
		array_order, local_best_score = decide_array_order_local_best(arrays_of_interest_dict, 100, args.iterations)
		print("The score (sum of spacers shared between all neighbouring arrays) of the best ordering found was: %i" % local_best_score)

## Find spacers present in more than one array

occurrences = dict(collections.Counter([item for sublist in list(arrays_of_interest_dict.values()) for item in sublist]))

imp_spacers = []

for k,v in occurrences.items():
	if v > 1:
		imp_spacers.append(k)

print("Identified %i spacers present in more than one array." % len(imp_spacers))

## Figure out colour scheme to use

if args.colour_file:
	with open(args.colour_file, 'r') as fin:
		colour_scheme = [i.strip() for i in fin.readlines()]
else:
	if len(imp_spacers) > 27:
		if len(imp_spacers) > 40:
			if len(imp_spacers) > 50:
				if len(imp_spacers) <100:
					colour_scheme = Cols_hex_100
					print("Using 100 colour colour scheme")
				else:
					print("Don't have enough colours for the number of spacers you are looking at. Colour scheme will be repeated")
					colour_scheme = []
					for i in range((len(imp_spacers)//100) + 1):
						colour_scheme += Cols_hex_100
					print("Using %i colour colour scheme" % len(colour_scheme))
			else:
				colour_scheme = Cols_hex_50
				print("Using 50 colour colour scheme")
		else:
			colour_scheme = Cols_hex_40
			print("Using 40 colour colour scheme")
	else:
		colour_scheme = Cols_hex_27
		print("Using 27 colour colour scheme")

spacer_colours = {}


largest_array_size = max([len(x) for x in [array_dict[y] for y in array_network]])

if args.legend:
	dim_x = largest_array_size*0.4
else:
	dim_x = largest_array_size*0.2
dim_y = len(array_network)*0.5


plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots()

fig.set_size_inches(dim_x, dim_y)

### Plot array alignments


colcount = 0 # count of which colour we are on in list
arcount = 1 # count of which array we are on (i.e. y axis value to plot at)
for array in array_order:
	spcount = 0 # count of which spacer in array we are on (i.e. which x-axis value to plot at)
	for spacer in arrays_of_interest_dict[array]:
		if spacer in imp_spacers:
			line_width= 3.0
			if spacer in spacer_colours.keys():
				spcolour = spacer_colours[spacer]
			else:
				spcolour = colour_scheme[colcount]
				colcount +=1
				spacer_colours[spacer] = str(spcolour)
		else: 
			spcolour = "#000000" #black
			line_width = 1.0
		ax.plot([spcount+0.5, spcount+1.5],[arcount, arcount], color = spcolour, linewidth = line_width, solid_capstyle="butt")
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
				sp_x1 = a+1
				sp_y1 = i+1
				sp_x2 = b+1
				sp_y2 = i+2
				spcolour = spacer_colours[spacer]
				ax.plot([sp_x1, sp_x2],[sp_y1, sp_y2], color = spcolour, linewidth = 1.0, solid_capstyle="butt", label = spacer)



axes = plt.gca()
plt.yticks(range(1, len(arrays_of_interest_dict.keys())+1, 1))
plt.xticks(range(1, max([len(x) for x in arrays_of_interest_dict.values()])+1, 1))
ax.set_yticklabels(array_order)


plt.xlabel("Spacer in array", fontsize=11)
plt.ylabel("Array ID", fontsize=11)

if args.legend:
	ratio_of_heights = ((len(imp_spacers)+1)*0.18)/dim_y #Approx height of each legend entry is 0.18 inches. Each array is allocated 0.5 inches. Add 1 to number of spacers to account for legend title.
	ncols = int(1+ ratio_of_heights) # Split the legend over multiple columns if it would be taller than the y axis. Plus 1 to account for rounding down.
	h_leg = 0.18 + (len(imp_spacers)*0.18)/ncols # What is the height of the legend now, having split it over columns?
	vert_pad = ((dim_y - h_leg)/2)/dim_y # Add blank space half the difference in height between legend and y axis to approximately center the legend relative to y axis. Then scale that on 0-1 relative to y axis size.
	# Build the legend manually as the automatic legend included duplicates as spacers are plotted once per array they are in.
	colors = list(spacer_colours.values()) 
	lines = [plt.Line2D([0], [0], color=c, linewidth=3, linestyle='--') for c in colors]
	labels = list(spacer_colours.keys())
	plt.legend(lines, labels, ncol=ncols, loc=(1.01,vert_pad), title="Spacer IDs") # Alternative legend positioning methods that don't work: bbox_to_anchor=(1.05, 1) loc="best",


fig.tight_layout()

plt.savefig(outfile, dpi=300)
plt.close(fig)

