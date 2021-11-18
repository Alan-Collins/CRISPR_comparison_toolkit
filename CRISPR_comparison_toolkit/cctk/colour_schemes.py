import sys
from itertools import combinations, permutations
from random import sample, randrange, seed
from math import sqrt, ceil

Cols_tol = [
	"#332288",
	"#117733",
	"#44AA99",
	"#88CCEE",
	"#DDCC77",
	"#CC6677",
	"#AA4499",
	"#882255",
	]

Cols_hex_12 = [
	"#07001c",
	"#ff6f8d",
	"#4c62ff",
	"#92ffa9",
	"#810087",
	"#bcffe6",
	"#490046",
	"#00c8ee",
	"#b53900",
	"#ff8cf7",
	"#5b5800",
	"#14d625",
	]

Cols_hex_27 = [
	"#fd5925",
	"#dbc58e",
	"#008d40",
	"#304865",
	"#934270",
	"#f7b8a2",
	"#907500",
	"#45deb2",
	"#1f4195",
	"#d67381",
	"#8e7166",
	"#afb200",
	"#005746",
	"#a598ff",
	"#8f0f1b",
	"#b96000",
	"#667f42",
	"#00c7ce",
	"#9650f0",
	"#614017",
	"#59c300",
	"#1a8298",
	"#b5a6bd",
	"#ea9b00",
	"#bbcbb3",
	"#00b0ff",
	"#cd6ec6",
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
		colours = [(i, "#000000") for i in Cols_tol]
		return colours


	elif ncolours <= 12:
		colours = [(i, "#000000") for i in Cols_hex_12]
		return colours


	elif ncolours <= 27:
		colours = [(i, "#000000") for i in Cols_hex_27]
		return colours


	elif ncolours <= 40:
		colours = [(i, "#000000") for i in Cols_hex_40]
		return colours


	elif ncolours < 65:
		col_scheme = Cols_tol
	elif ncolours < 145:
		col_scheme = Cols_hex_12
	else:
		cs_name = "built-in"
		col_scheme = Cols_hex_27
	
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
			"\n\nWARNING!!! "
			"There are not enough colours in the {} colour scheme.\n"
			"The colour scheme will be repeated so some spacers will "
			"be coloured the same.\nIn order to colour each spacer uniquely, "
			"{} colours are needed\n".format(
				cs_name,
				ceil(sqrt(ncolours))))
	while ncolours > len(colours):
		colours += colours

	# Reset seed
	seed(None)
		
	return colours
