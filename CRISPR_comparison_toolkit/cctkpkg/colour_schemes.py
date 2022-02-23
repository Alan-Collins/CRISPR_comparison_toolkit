import sys
from itertools import combinations, permutations
from random import sample, randrange, seed, randint
from math import sqrt, ceil

from . import file_handling


Cols_wong = [
	"#000000",
	"#e69f00",
	"#56B4E9",
	"#009E73",
	"#F0E442",
	"#0072B2",
	"#D55E00",
	"#CC79A7"
	]


def choose_col_scheme(ncolours, s=None, cf_list=None):
	"""Identify colour scheme with enough colours for all spacers.

	Args:
	  ncolours (int):
		Number of colours desired in produced colour scheme.
	  s (int or None):
		Seed to control random ordering of colours in colour scheme.
	  cf_list (list):
	    User-provided colours.
	  
	Returns:
	  colours (list of tuples):
	    list of two-colour combinations for use as fill and outline.
	    e.g.: colours = [
			('#07001c', '#00c8ee'),
			('#810087', '#00c8ee'),
			('#490046', '#92ffa9'),
			...,
			]
	"""

	if cf_list:
		cs_name = "user-provided"
		if ncolours <= len(cf_list):
			colours = [(i, "#000000") for i in cf_list]
			return colours

		else:	
			col_scheme = cf_list	
	
	elif ncolours <= 8:
		colours = [(i, "#000000") for i in Cols_wong]
		return colours


	elif ncolours <= 64:
		col_scheme = Cols_wong
	else:
		cs_name = "built-in"
		col_scheme = Cols_wong
	
	colours = []
	# Repeat the same colour scheme for fill and outline combos.
	seed(s)
	combos = [i for i in permutations(col_scheme, 2)]
	# Permutations doesn't return self to self so add those.
	combos += [(i,i) for i in col_scheme]
	colours = sample(combos, len(combos))

	# If colour scheme is insufficient, repeat until enough.
	if ncolours > len(colours):
		sys.stderr.write(
			"\nWARNING!!! "
			"There are not enough colours in the {} colour scheme.\n"
			"Random colours will be used. "
			"\nIn order to colour each spacer uniquely, "
			"{} colours are needed\n\n".format(
				cs_name,
				ceil(sqrt(ncolours))))
		colours += random_colour_pairs(ncolours-len(colours), s)

	# Reset seed
	seed(None)
		
	return colours


def random_colour_pairs(n, s=None):
	"""Return pairs of random hex code colours with seed control
	
	Args:
	  n (int):
	    Number of pairs of colours to return.
	  s(int):
	    value to control random number generation.

	Returns:
	  colours (list of tuples):
	    A list of tuples of randomly generated pairs of colours.

	"""

	seed(s)
	colours = [
		("#%06x" % randint(0,0xFFFFFF), "#%06x" % randint(0,0xFFFFFF))
		for i in range(n)
		]
	seed(None)

	return colours


def fill_col_scheme_gaps(col_scheme, spacers, seed):
	for s in spacers:
		if s not in col_scheme.keys():
			rand_col = random_colour_pairs(1, seed)[0]
			while rand_col in col_scheme.values():
				seed+=1
				rand_col = random_colour_pairs(1, seed)[0]
			col_scheme[s] = rand_col

	return col_scheme


def process_colour_args(args, non_singleton_spacers):
	if args.colour_scheme_infile:
		spacer_cols_dict = file_handling.read_colour_scheme(
			args.colour_scheme_infile)
		# If necessary add colours to colour scheme for missing spacers
		if any([s not in spacer_cols_dict for s in non_singleton_spacers]):
			spacer_cols_dict = fill_col_scheme_gaps(
				spacer_cols_dict, non_singleton_spacers, args.seed)

	else:
		# Else, figure out colour scheme to use
		if args.colour_file:
			cf_list = file_handling.read_colour_file(args.colour_file)
			colours = choose_col_scheme(
				len(non_singleton_spacers),
				args.seed,
				cf_list)
		else:
			colours = choose_col_scheme(
				len(non_singleton_spacers),
				args.seed)

		# build a dictionary with colours assigned to each spacer.
		spacer_cols_dict = {}
		for i, spacer in enumerate(sorted(non_singleton_spacers)):
			spacer_cols_dict[spacer] = colours[i]


	if args.colour_scheme_outfile:
		file_handling.write_colour_scheme(
			args.colour_scheme_outfile, spacer_cols_dict)

	return spacer_cols_dict
