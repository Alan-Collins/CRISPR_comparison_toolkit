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


array_dict = {}
with open(args.array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_dict[bits[0]] = bits[2:]



def hamming(seq1, seq2):
	"""
	Args:
		seq1 (str or list): First sequence to compare.
		seq2 (str or list): Second sequence to compare.
	
	Returns:
		(int) The hamming distance between the two sequences.
	"""

	dist = 0

	for i in range(max(len(seq1), len(seq2))):
		if i < len(seq1) and i < len(seq2):
			if seq1[i] != seq2[i]:
				dist += 1
		else:
			dist += 1

	return dist


def make_dist_mat(seqs):
	"""
	Args:
		seqs (list): List of sequences to be compared.
	
	Returns:
		(numpy.array) numpy array of the pairwise distances between arrays.
	"""

	n = len(seqs)

	grid = np.zeros((n,n))
	
	# Iterate across all combinations of arrays
	for i in range(n):
		for j in range(i+1, n):
			seq1 = seqs[i]
			seq2 = seqs[j]

			aln1, aln2 = needle(seq1, seq2)
			d = hamming(aln1, aln2)
			# Assign distance to both parts of symetrical distance matrix
			grid[i][j] = d
			grid[j][i] = d

	return grid

def make_q_mat(mat):
	""" 
	Make a matrix of Q values based on the formula on the wikipedia page https://en.wikipedia.org/wiki/Neighbor_joining
	Q[i,j] = (n-2)d[i,j] - sum(d[i]) - sum(d[j]) 
	Args:
		mat (numpy.array): Matrix of pairwise distances.
	
	Returns:
		(numpy.array) Matrix of pairwise Q values.
	"""

	# make matrix of 0s to fill in
	q = np.zeros(mat.shape)
	n = mat.shape[0]

	# Iterate across matrix
	for i in range(n):
		for j in range(i+1, n):
			qval = (n-2)*mat[i,j] - sum(mat[i]) - sum(mat[j])
			q[i,j] = qval
			q[j,i] = qval


	return q




arrays = [array_dict[i] for i in args.arrays_to_join]
dmat = make_dist_mat(arrays)

print(make_q_mat(dmat))
