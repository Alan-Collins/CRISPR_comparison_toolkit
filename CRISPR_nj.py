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
	"-p", dest="print_tree", action='store_true',  
	help="Print a graphical representation of the tree using ascii characters (required ete3 to be installed)."
	)
parser.add_argument(
	'-d', dest="distance_model", choices=['hamming', 'evol'], default='hamming'
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
			if args.distance_model == 'hamming':
				d = hamming(aln1, aln2)
			elif args.distance_model == 'evol':
				d = CRISPR_evol_model_dist(aln1, aln2)
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


def find_best_pair(qmat):
	"""
	Find the closest pair of items in a Q matrix
	Args:
		qmat (numpy.array): Matrix of Q values.
	
	Returns:
		(tuple of ints) Index of the best scoring pair.
	"""

	best_score = 999999
	best_idx = 9999999

	n = qmat.shape[0]
	for i in range(n):
		for j in range(i+1, n):
			if qmat[i,j] < best_score:
				best_score = qmat[i,j]
				best_idx = (i,j)

	return best_idx


def make_new_dist_mat(old_mat, index):
	"""
	Args:
		oldmat (numpy.array): distance matrix used to find best pair.
		index (tuple of ints): indices of best pair.
	
	Returns:
		(numpy.array) New distance matrix with joined pair collapsed to a single row/column in the matrix.
	"""
	old_n = old_mat.shape[0]
	new_mat = np.zeros((old_n-1, old_n-1))
	new_n = old_n-1

	x = 1
	y = 1

	for i in range(old_n):
		if i != index[0] and i != index[1]:
			for j in range(old_n):
				if j != index[0] and j != index[1]:
					new_mat[x][y] = old_mat[i][j]
					y += 1
			y = 1
			x += 1

	for i in range(new_n):
		new_mat[0,i] = (old_mat[index[0]][i+1] + old_mat[index[1]][i+1] - old_mat[index[0],index[1]])/2
		new_mat[i,0] = (old_mat[index[0]][i+1] + old_mat[index[1]][i+1] - old_mat[index[0],index[1]])/2
	
	return new_mat



def calc_dist_to_node(dmat, index):
	i,j = index
	if dmat.shape[0] != 2:
		dist_i = (dmat[i][j])/2 + 1/(2*(dmat.shape[0]-2))*(sum(dmat[i])-sum(dmat[j]))
	else:
		dist_i = 0
	dist_j = dmat[i][j] - dist_i

	return dist_i, dist_j


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

	dist = n_aq*event_cost_dict['sp_aq'] + n_ec_aq*event_cost_dict['ectopic_sp_aq'] + n_del*event_cost_dict['deletion']

	return dist


arrays = [array_dict[i] for i in args.arrays_to_join]
labels = args.arrays_to_join
dmat = make_dist_mat(arrays)


while dmat.shape[0] > 1: # Keep joining neighbours until the tree is fully resolved
	qmat = make_q_mat(dmat)
	best_idx = find_best_pair(qmat)
	dists = calc_dist_to_node(dmat, best_idx)
	# To form a newick tree structure, combine the labels of the best_idx into a sublist in the list of labels
	a = str(labels[best_idx[0]]).replace('[','(').replace(']',')').replace("'","")
	b = str(labels[best_idx[1]]).replace('[','(').replace(']',')').replace("'","")
	newnode_labels = [["{}:{}".format(a, dists[0]), "{}:{}".format(b, dists[1])]]
	labels = newnode_labels + [i for i in labels if i != labels[best_idx[0]] and i != labels[best_idx[1]]]
	dmat = make_new_dist_mat(dmat, best_idx)

# Convert sublists into newick format by replacing square brackets with parentheses
tree = str(labels)[1:-1].replace('[','(').replace(']',')').replace("'","")+';'

print(tree)

if args.print_tree:
	from ete3 import Tree
	print(Tree(tree))
