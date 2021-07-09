#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-7-7
# DESCRIPTION :  Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship.

import sys
import argparse
import numpy as np


parser = argparse.ArgumentParser(
	description="Perform maximum parsimony analysis on CRISPR arrays to infer a tree representing their evolutionary relationship."
	)
parser.add_argument(
	"-a", dest="array_file", required = True,
	help="Specify array representatives file."
	)
parser.add_argument(
	"-p", dest="print_tree", action='store_true',  
	help="Print a graphical representation of the tree using ascii characters (required ete3 to be installed)."
	)
parser.add_argument(
	"arrays_to_join", nargs="+",  
	help="Specify the IDs of the arrays you want to join. **Must come at the end of your command after all other arguments.**"
	)


args = parser.parse_args(sys.argv[1:])

class Array():
	"""
	Class to store information about extant and inferred ancestral CRISPR arrays to aid in their comparisons.
	
	Attributes:
		extant (bool): A boolean indicating if the array is extant in our dataset or if it was inferred as a hypothetical ancestral state.
		chunks (list): A list of the contiguous blocks of spacers with common features (e.g. consecutive spacers that are absent in aligned array).
	"""
	def __init__(self, extant=True):
		self.extant = extant
		self.chunks = []

class Spacer_Chunk(object):
	"""
	Class to store information about spacers in CRISPR arrays.
	
	Attributes:
		singleton (bool): A boolean indicating if this spacer chunk is only found in a single array.
		type (str): A string indicating the nature of this chunk. e.g. indel, aqcuisition, shared...
		spacers (list): A list of the spacer IDs in this chunk.
		linked (bool): A boolean indicating whether these spacers are always found together in arrays when one of them is found.
	"""
	def __init__(self):
		self.singleton = arg
		self.type = arg
		self.spacers = arg
		self.linked = arg
		


def needle(seq1, seq2, match = 1, mismatch = -1, gap = -2):
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

	i = len(seq2)
	j = len(seq1)

	# Read back through the grid along the best path to create the best alignment
	align1, align2 = [], []
	while i > 0 and j > 0: # end when it reaches the top or the left edge
		score_current = grid[i][j]
		score_diagonal = grid[i-1][j-1]
		score_up = grid[i][j-1]
		score_left = grid[i-1][j]
		if seq1[j-1] == seq2[i-1]:
			score = match
		else:
			score = mismatch
		# Check to figure out which cell the current score was calculated from,
		# then update i and j to correspond to that cell.
		if score_current == score_diagonal + score:
			align1.append(seq1[j-1])
			align2.append(seq2[i-1])
			i -= 1
			j -= 1
		elif score_current == score_up + gap:
			align1.append(seq1[j-1])
			align2.append('-')
			j -= 1
		elif score_current == score_left + gap:
			align1.append('-')
			align2.append(seq2[i-1])
			i -= 1

	# Finish tracing up to the top left cell
	while j > 0:
		align1.append(seq1[j-1])
		align2.append('-')
		j -= 1
	while i > 0:
		align1.append('-')
		align2.append(seq2[i-1])
		i -= 1
	
	# Since we traversed the score matrix backwards, need to reverse alignments.
	align1 = align1[::-1]
	align2 = align2[::-1]

	if isinstance(seq1, str) and isinstance(seq2, str):
		align1 = ''.join(align1)
		align2 = ''.join(align2)
	
	return align1, align2


def CRISPR_evol_model_dist(s1, s2):
	"""
	Args:
		s1 (str or list): The first aligned sequence for which you want a distance.
		s2 (str or list): The second aligned sequence for which you want a distance.
	
	Returns:
		(int) Distance.
	"""

	event_cost_dict = {
	'sp_aq' : 1, # Acquisition of a spacer at the leader end
	'ectopic_sp_aq' : 1, # Acquisition of a spacer inside the array
	'deletion' : 1 # Loss of some number of spacers
	}

	leader = True
	gap = False
	mismatch = False

	n_aq = 0
	n_ec_aq = 0
	n_del = 0

	gaps = []
	gaps_n = []
	gap_size = 0
	mismatch_size = 0

	for n, (i,j) in enumerate(zip(s1,s2)):
		if i == '-' or j == '-':
			if leader:
				n_aq += 1
			else:
				if gap:
					gap_size += 1
				else:
					if mismatch: # Mismatched spacer before a gap may be an ectopic acquisition that aligns wierdly due to the gap.
						if mismatch_size == 1:  
							n_ec_aq += 1
						else: # If it's more than 1 spacer long it's unlikely to be multiple ectopic acquisitions in the same place.
							n_del += 2
						mismatch = False
					gap = True
					gaps_n.append(n)
					gap_size += 1
		elif i == j:
			if gap:
				gaps.append(gap_size)
				gap = False
				gap_size = 0
			if leader:
				leader = False
			if mismatch: # If two spacers mismatch before a region of identity, then there must have been an indel in each
						 # or ectopic acquisitions at the same position. Deletions may be more likely.
				n_del += 2
				mismatch = False
		else: # if i != j
			if leader:
				n_aq += 2
			else:
				if mismatch:
					mismatch_size += 1
				else:
					mismatch = True
					mismatch_size = 1
		if i == s1[-1]:
			if mismatch:
				n_del += 2 # If the trailer end has different spacers then perhaps either they recombined or both come from an ancestral longer array

	for n, gap in zip(gaps_n, gaps):
		if gap == 1:
			if n < 10:
				n_ec_aq += 1
			else:
				n_del += 1
		else:
			n_del += 1

	return n_aq, n_ec_aq, n_del


def infer_ancestor(seq1, seq2, all_spacers):
	"""
	Args:
		seq1 (str or list): The first sequence to be compared.
		seq2 (str or list): The second sequence to be compared.
		all_spacers (list): The list of all arrays so that indels can be resolved to favour keeping spacers found in other arrays.
	
	Returns:
		(str or list) A hypothesis of the ancestral state of the provided sequences.
	"""

	ancestor = None

	aln1, aln2 = needle(seq1, seq2)

	print(aln1)
	print(aln2)

	return ancestor


array_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]

arrays = [array_dict[i] for i in args.arrays_to_join]
labels = args.arrays_to_join

print(infer_ancestor(array_dict['1338'], array_dict['1285'], arrays))


if args.print_tree:
	from ete3 import Tree
	#print(Tree(tree))
