#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-7-6
# DESCRIPTION :  Perform neighbour joining on a given group of CRISPR arrays

import sys
import argparse
import numpy as np


parser = argparse.ArgumentParser(
	description="Perform neighbour joining on a given group of CRISPR arrays"
	)
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"arrays_to_join", nargs="+",  
	help="Specify the IDs of the arrays you want to join. **Must come at the end of your command after all other arguments.**"
	)



args = parser.parse_args(sys.argv[1:])



def needle(seq1, seq2, match = 1, mismatch = -1, gap = -1):
	"""
	Args:
		seq1 (str or list): First sequence of items to align.
		seq2 (str or list): Second sequence of items to align
		match (int): Score for match at a position in alignment.
		mismatch(int): Penalty for mismatch at a position in alignment.
		gap (int): Penalty for a gap at a position in alignment.
	
	Returns:
		(tuple of str lists) Returns a tuple containing the input seq1 and seq2 aligned with '-' added as gaps.
		If strings were given then strings are returned. If lists were given then lists are returned.

	"""

