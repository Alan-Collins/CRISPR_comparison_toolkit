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

	# Make a list of lists of 0s with dimensions x by y: list containing x lists of y 0s each.
	grid = np.zeros((len(seq2)+1, len(seq1)+1))

	# Fill in grid with scores for all possible alignments
	# First score for no alignment (i.e. all gaps)
	for i in range(len(seq1)+1):
		grid[0][i] = gap*i
	for i in range(len(seq2)+1):
		grid[i][0] = gap*i

	# Then score for each cell if you came to it from the nearest best cell/
	for i in range(len(seq1)):
		for j in range(len(seq2)):
			if seq1[i] == seq2[j]:
				score = match
			else:
				score = mismatch
			grid[j+1][i+1] = max([grid[j][i]+score, grid[j+1][i]+gap, grid[j][i+1]+gap])

	print(grid)

print(needle("ATCG", "CCATGG"))
