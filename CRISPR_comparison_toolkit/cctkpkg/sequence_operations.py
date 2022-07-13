import sys
import os
import numpy as np
from collections import defaultdict, Counter
from itertools import combinations
import subprocess
from copy import deepcopy
import multiprocessing

from . import (
	file_handling
	)

SEQUENCE_DICT = { # IUPAC DNA alphabet as regex
	"A": "A",
	"T": "T",
	"C": "C",
	"G": "G",
	"N": "[ATCG]",
	"W": "[AT]",
	"S": "[CG]",
	"M": "[AC]",
	"K": "[TG]",
	"R": "[AG]",
	"Y": "[TC]",
	"B": "[CTG]",
	"D": "[ATG]",
	"H": "[ACT]",
	"V": "[ACG]"
}


class NetworkEdge():
	"""Class to store edge attributes in Network
	
		Initialized with either array_parsimony.Array or 
		file_handling.FoundArray instance
	
		Attributes:
		  a (str):
		    First array (node)
		  b (str):
		    Second array (node)
		  nshared (int):
		    number of shared spacers
		  jaccard (float):
		    Jaccard similarity index between arrays a and b
		  a_type (str):
		    Array a CRISPR type
		  b_type (str):
		     Array b CRISPR type
		"""
	def __init__(self, a, b):
		self.a = a.id
		self.b = b.id
		self.nshared = len(set(a.spacers).intersection(set(b.spacers)))
		self.jaccard = len(set(a.spacers).intersection(
			set(b.spacers)))/len(set(a.spacers).union(set(b.spacers)))
		self.a_type = a.repeat_id if "repeat_id" in vars(a) else ''
		self.b_type = b.repeat_id if "repeat_id" in vars(b) else ''


def rev_comp(string):
	"""Reverse complement a string of nucleotide sequence
	
	Args:
	  string (str):
		  Nucleotide sequence

	Returns:
	  str:
		reverse complement of given nucleotide sequence
	
	Raises:
	  TypeError: If string is not str.

	"""
	if type(string) is not str:
		raise TypeError(
			"string must be str, not {}.".format(type(string).__name__))

	rev_str = ''
	rev_comp_lookup = {
	"A" : "T", 
	"T" : "A", 
	"C" : "G", 
	"G" : "C", 
	"a" : "t", 
	"t" : "a", 
	"c" : "g", 
	"g" : "c",
	}
	for i in reversed(string):
		if i in "ATCGatcg":
			rev_str += rev_comp_lookup[i]
		else:
			rev_str += i
	return rev_str


def hamming(string1, string2):
	"""Calculate hamming distance between two sequences. 

	Compare them as provided, and with one sequence	shifted one position
	to the left, and with the other sequence shifted one to the left.
	
	Args:
	  string1 (str):
		First sequence to be compared
	  string2 (str):
		Second sequence to be compared

	Returns:
	  tuple:
		dist (int): 
		  The hamming distance of the unadjusted sequences
		distplus1 (int):
		  The hamming distance of the first sequence shifted
		distminus1 (int):
		  The hamming distance of the econd sequence shifted

	Raises:
	  TypeError: If string1 and string2 are not str.

	"""
	if type(string1) is not str:
		raise TypeError(
			"Inputs must be str, not {}.".format(type(string1).__name__))
	if type(string2) is not str:
		raise TypeError(
			"Inputs must be str, not {}.".format(type(string2).__name__))
	dist, distplus1, distminus1 = 0,0,0
	string1plus1 = string1[1:]
	string2plus1 = string2[1:]
	for i in range(max(len(string1), len(string2))):
		if i < len(string1) and i < len(string2):
			if string1[i] != string2[i]:
				dist += 1
		else:
			dist += 1
	
	for i in range(max(len(string1plus1), len(string2))):
		if i < len(string1plus1) and i < len(string2):
			if string1plus1[i] != string2[i]:
				distplus1 += 1
		else:
			distplus1 += 1
	
	for i in range(max(len(string1), len(string2plus1))):
		if i < len(string1) and i < len(string2plus1):
			if string1[i] != string2plus1[i]:
				distminus1 += 1
		else:
			distminus1 += 1

	return dist, distplus1, distminus1


def needle(seq1, seq2, match = 100, mismatch = -1, gap = -2):
	"""
	Perform Needleman-Wunsch pairwise alignment of two sequences.
	Args:
	  seq1 (str, list, or tuple):
		First sequence of items to align.
	  seq2 (str, list, or tuple):
		Second sequence of items to align
	  match (int):
		Score for match at a position in alignment.
	  mismatch(int):
		Penalty for mismatch at a position in alignment.
	  gap (int):
		Penalty for a gap at a position in alignment.
	
	Returns:
	  (tuple of str lists) Returns a tuple containing the input 
	  seq1 and seq2 aligned with '-' added as gaps.
	  If strings were given then strings are returned.
	  If lists were given then lists are returned.

	Raises:
	  TypeError: If seq1 and seq2 are not str, list, or tuple.
	  TypeError: If match, mismatch, and gap are not int or float.
	"""

	if type(seq1) not in [str, list, tuple]:
		raise TypeError(
			"Inputs must be str, list, or tuple, not {}.".format(
				type(seq1).__name__))
	if type(seq2) not in [str, list, tuple]:
		raise TypeError(
			"Inputs must be str, list, or tuple, not {}.".format(
				type(seq2).__name__))

	if type(match) not in [int, float]:
		raise TypeError(
			"match must be int or float, not {}.".format(
				type(match).__name__))
	if type(mismatch) not in [int, float]:
		raise TypeError(
			"mismatch must be int or float, not {}.".format(
				type(mismatch).__name__))
	if type(gap) not in [int, float]:
		raise TypeError(
			"gap must be int or float, not {}.".format(
				type(gap).__name__))


	# Make a list of lists of 0s with dimensions x by y: 
	# list containing x lists of y 0s each.
	grid = np.zeros((len(seq2)+1, len(seq1)+1))

	# Fill in grid with scores for all possible alignments
	# First score for no alignment (i.e. all gaps)
	for i in range(len(seq1)+1):
		grid[0][i] = gap*i
	for i in range(len(seq2)+1):
		grid[i][0] = gap*i

	# Then score for each cell if you came to it from the 
	# nearest best cell
	for i in range(len(seq1)):
		for j in range(len(seq2)):
			if seq1[i] == seq2[j]:
				score = match
			else:
				score = mismatch
			grid[j+1][i+1] = max([grid[j][i]+score,
				grid[j+1][i]+gap,
				grid[j][i+1]+gap])

	i = len(seq2)
	j = len(seq1)

	# Read back through the grid along the best path to create the
	# best alignment
	align1, align2 = [], []
	# end when it reaches the top or the left edge
	while i > 0 and j > 0:
		score_current = grid[i][j]
		score_diagonal = grid[i-1][j-1]
		score_up = grid[i][j-1]
		score_left = grid[i-1][j]
		if seq1[j-1] == seq2[i-1]:
			score = match
		else:
			score = mismatch
		# Check to figure out which cell the current score was 
		# calculated from, then update i and j to correspond to
		# that cell.
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
	
	# Since we traversed the score matrix backwards, need to
	# reverse alignments.
	align1 = align1[::-1]
	align2 = align2[::-1]

	if isinstance(seq1, str) and isinstance(seq2, str):
		align1 = ''.join(align1)
		align2 = ''.join(align2)
	
	return align1, align2


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


def get_repeat_info(CRISPR_types_dict, repeat):
	"""Find best repeat match and assign CRISPR type accordingly.

	Figure out which of the repeats the user provides in the repeats 
	fasta file is the most similar to the repeat found by minced

	Args:
	  CRISPR_types_dict (dict):
		Fasta dict of the repeats fasta file 
		as returned by fasta_to_dict
	  repeat (str):
		The repeat identified by minced for which you want to find 
		a match

	Returns:
	  tuple: 
		best_match (str):
		  The fasta header of the user-provided repeat that is most 
		  similar to that found by minced
		best_score (int):
		  The hamming distance between the best match repeat and the 
		  minced-identified repeat
		reverse (bool):
		  Is the best match repeat the reverse complement of 
		  the identified repeat?
	"""
	best_score = 1000
	best_match = ''
	reverse = False
	for k,v in CRISPR_types_dict.items():
		score = min(hamming(repeat, v))
		if score < best_score:
			best_score = score
			best_match = k
			reverse = False

		score = min(hamming(
			rev_comp(repeat), v))
		if score < best_score:
			best_score = score
			best_match = k
			reverse = True


	return best_match, best_score, reverse


def non_redundant_CR(
	all_assemblies,
	snp_thresh=0,
	prev_spacer_id_dict={},
	prev_array_dict={},
	outdir="./"
	):
	""" Identify non-redundant spacers and arrays

		Given a list of file_handling.AssemblyCRISPRs instances and some
		other parameters, return dereplicated CRISPR objects

		Args:
		  all_assemblies (list of file_handling.AssemblyCRISPRs):
		    The CRISPR information to be dereplicated
		  snp_thresh (int):
		    The threshold of SNPs between CRISPR spacers below which
		    spacers should be considered the same
		  prev_spacer_id_dict (dict):
		    If appending to existing dataset, dict of spacers associated
		    with each CRISPR type. E.g.,:
		    {"CRISPR_type": ["list", "of", "spacer_ids"]}
		  prev_array_dict (dict):
		    If appending to existing dataset, dict of array IDs and the
		    spacers in those arrays. E.g.,:
		    {"sp1 sp2 sp3 sp4": "Array_ID"}

		Returns:
		  tuple:
		    non_red_spacer_dict (dict):
		      All spacers found for each CRISPR type.
		        {"CRISPR_type": ["list", "of", "spacer", "sequences"]}
		    non_red_spacer_id_dict (dict):
		      All spacers and their associated ID.
		        {"spacer_sequence": "spacer_ID"}
		    non_red_array_dict (dict):
		      Arrays found for each CRISPR type.
		        {"CRISPR_type": ["list", "of", "array", "spacers"]}
		    non_red_array_id_dict (dict):
		      Array IDs of identified arrays.
		        {"sp1 sp2 sp3 ...": "Array_ID"}
		    cluster_reps_dict (dict of dicts):
		      If any spacers were identified that differ by fewer than
		      the SNP threshold, this dict describes the representative
		      chosen for each group of similar spacers. E.g.,
		        {"CRISPR_type": {"representative": [list of spacers]}
		    rev_cluster_reps_dict:
		      The reverse of the above dict. For each spacer that was 
		      removed, what is its representative. E.g.,
		        {"CRISPR_type": {"spacer": "Representative"}}



	"""
	non_red_spacer_dict = defaultdict(list)
	all_arrays_dict = defaultdict(list)
	for strain in all_assemblies:
		for k,v in strain.arrays.items():
			CR_type = v.repeat_id
			all_arrays_dict[CR_type].append(" ".join(v.spacers))
			for s in v.spacers:
				if s not in set(non_red_spacer_dict[CR_type]):
					non_red_spacer_dict[CR_type].append(s)

	cluster_reps_dict = {}
	rev_cluster_reps_dict = defaultdict(dict)
	if snp_thresh > 0:
		spacer_snp_network = make_spacer_snp_network(
			non_red_spacer_dict, snp_thresh=snp_thresh, outdir=outdir)
		for CR_type, network in spacer_snp_network.items():
			clusters = identify_network_clusters(network)
			cluster_reps_dict[CR_type] = pick_cluster_rep(
				clusters, non_red_spacer_dict[CR_type], prev_spacer_id_dict)
			# Add to the reverse dict to use for replacing instances of 
			# cluster members with their rep later
			for rep, cluster_members in cluster_reps_dict[CR_type].items():
				for member in cluster_members:
					rev_cluster_reps_dict[CR_type][member] = rep

	# empty if not appending
	non_red_spacer_id_dict = deepcopy(prev_spacer_id_dict)
	
	for CR_type, spacers in non_red_spacer_dict.items():
		if snp_thresh > 0:
			# Remove cluster members if appropriate
			spacers = [
				sp for sp in spacers if sp not in rev_cluster_reps_dict[CR_type]]
		
		# If appending figure out highest spacer number of this CR_type
		sp_num = 0
		if len(prev_spacer_id_dict) > 0:
			for prev_id in prev_spacer_id_dict.values():
				# Skip spacers not of this type
				if CR_type not in prev_id:
					continue
				prev_id = int(prev_id.split("_")[1])
				if prev_id > sp_num:
					sp_num = prev_id
		
		for seq in spacers:
			if seq in non_red_spacer_id_dict:
				continue

			sp_num += 1
			non_red_spacer_id_dict[seq] = CR_type+"_"+str(sp_num)

	if snp_thresh > 0:
		# Replace cluster members with their rep if appropriate
		for CR_type,v in all_arrays_dict.items():
			for n, array in enumerate(v):
				if not any(
					[sp in rev_cluster_reps_dict[CR_type] for sp in array.split()]):
					continue
				new_array = []
				for sp in array.split():
					if sp in rev_cluster_reps_dict[CR_type]:
						new_array.append(rev_cluster_reps_dict[CR_type][sp])
					else:
						new_array.append(sp)
				all_arrays_dict[CR_type][n] = " ".join(new_array)


	non_red_array_dict = {k: set(v) for k,v in all_arrays_dict.items()}
	all_array_list = [a for ar_ls in non_red_array_dict.values() for a in ar_ls]

	# empty if not appending
	non_red_array_id_dict = deepcopy(prev_array_dict)

	# If appending figure out highest array number
	array_num = 0
	if len(prev_array_dict) > 0:
		array_num = max([int(i) for i in prev_array_dict.values()])

	# Sort arrays for reproducibility
	for spacers in sorted(all_array_list):
		if spacers in non_red_array_id_dict:
			continue

		array_num += 1
		non_red_array_id_dict[spacers] = array_num

	return (non_red_spacer_dict,
		non_red_spacer_id_dict,
		non_red_array_dict,
		non_red_array_id_dict,
		cluster_reps_dict,
		rev_cluster_reps_dict)


def add_ids(
	all_assemblies,
	non_red_spacer_id_dict,
	non_red_array_id_dict,
	rev_cluster_reps_dict):
	""" Lookup array and spacer IDs and add them to AssemblyCRISPRs 
	"""
	for assembly in all_assemblies:
		if assembly.array_count == 0:
			continue
		for array in assembly.arrays.values():
			for n, spacer in enumerate(array.spacers):
				if spacer in rev_cluster_reps_dict[array.repeat_id]:
					spacer = rev_cluster_reps_dict[array.repeat_id][spacer]
					array.spacers[n] = spacer
				array.spacer_ids.append(non_red_spacer_id_dict[spacer])
			array.id = non_red_array_id_dict[" ".join(array.spacers)]


def build_network(array_list):
	""" Assess spacer sharing between arrays. Add as edge if any shared
	"""
	network = []
	for a,b in combinations(array_list, 2):
		nshared = len(set(a.spacers).intersection(set(b.spacers)))
		if nshared > 0:
			network.append(NetworkEdge(a,b))

	return network


def percent_id(a, b):
	l = len(a)
	matches = 0
	for x,y in zip(a,b):
		if x==y:
			matches+=1

	pident = 100*matches/l

	return pident


def pool_MP_blastdbcmd(inputs, db, batch_size, threads):
	"""
	Manages the multiprocess worker pool and returns to seqs found by
	run_blastcmd to whatever called it.
	Args:
		inputs (list of tuples): all inputs to be run 
			(n, fstr, batch_locs)
		db (str): path to the blast db you want to query.
		batch_size (int): Max size of batches to be run by blastdbcmd 
		threads (int): Number of threads to use for multiprocessing.
		
	
	Returns:
		(list) List of FoundArray instances representing distinct arrays
	"""
	pool = multiprocessing.Pool(processes=threads)

	batch_size = min(batch_size, len(inputs)/threads)

	chunksize = int((len(inputs)/batch_size)//threads)

	batches = []
	for i in range(0, len(inputs), batch_size):
		ns = []
		fstring = ""
		batch_locations = ""
		for j in range(i, min(i+batch_size, len(inputs))):
			ns.append(inputs[j][0])
			fstring += inputs[j][1]
			batch_locations += inputs[j][2]
		batches.append((ns, db, fstring, batch_locations))

	output = pool.starmap(run_blastcmd, batches, chunksize)
	output = [i for i in output]
	pool.close()
	pool.join()


	# Ensure the data are sorted with respect to the blast results
	output.sort(key=lambda batch: batch[0][0])
	seqs = []
	for batch in output:
		seqs += batch[1]

	return seqs


def run_blastcmd(ns, db, fstring, batch_locations):
	"""
	function to call blastdbcmd in a shell and process the output. Uses
	a batch query of the format provided in the blastdbcmd docs. e.g.
	printf "%s %s %s %s\\n%s %s %s\\n" 13626247 40-80 plus 30 14772189 \
	1-10 minus | blastdbcmd -db GPIPE/9606/current/all_contig \
	-entry_batch -
	
	Args:
		ns (list): indices to reorder the output from this func when
		  running on multiple threads
		db (str): path to the blast db you want to query.
		fstring (str):  The string to give to printf
		batch_locations (str):  the seqid, locations, and strand of all
		  the spacers to be retrieved
	
	Returns:
		(list) List of the sequences of regions returned by blastdbcmd
		  in response to the query locations submitted
	
	Raises:
		ERROR running blastdbcmd: Raises an exception when blastdbcmd
		  returns something to stderr, prints the error as well as the
		  information provided to blastdbcmd and aborts the process.
	"""
	x = subprocess.run(
		'printf "{}" {}| blastdbcmd -db {} -entry_batch -\
			'.format(fstring, batch_locations, db), 
		shell=True,
		universal_newlines=True,
		capture_output=True
		) 
	if x.stderr:
		print("ERROR running blastdbcmd on {} :\n{}".format(
			db, batch_locations, x.stderr))
		sys.exit()
	else:
		return (ns, 
			[i for i in x.stdout.split('\n') if '>' not in i and len(i) > 0])


def determine_regex_length(pattern):
	""" Calculate minimum length of sequence that a regex can describe
	"""
	minlen = 0

	i = 0
	while i < len(pattern):
		c = pattern[i]
		if c.isalpha():
			# If character is a base then 1 position required
			minlen+=1
			i+=1
			continue
		
		if c == "[":
			# Find end of set to determine size
			while c != "]":
				i+=1
				if i > len(pattern)-1:
					sys.stderr.write("Error: Unmatched [ in your regex")
					sys.exit()
				c = pattern[i]

		if c == "\\":
			# Following character should be [A-Za-z] or this is an
			# inappropriate escape character
			i+=1
			c = pattern[i]
			if not c.isalpha():
				sys.stderr.write("Escape character used to include the "
					"following illegal character in your regex: {}".format(c))
				sys.exit()

		if i == len(pattern)-1:
			# If this is the end then no quantifier exists
			minlen+=1
			break
		i+=1
		c = pattern[i]
		if c not in [".", "+", "?", "{"]:
			# If not quantifier then set has length 1. continue
			minlen+=1
			continue

		if c in [".", "?"]:
			i+=1
			continue
		if c == "+":
			minlen+=1
			i+=1
			continue
		# Else quantifier is range or exact count.
		# minimum length of it is first number
		i+=1
		c = pattern[i]
		quant_len = ''
		while c.isdigit():
			quant_len+=c
			i+=1
			c= pattern[i]
		# Ends once , or } reached
		minlen+=int(quant_len)

		# continue until } if not reached
		while c != "}":
			i+=1
			c = pattern[i]

		if i == len(pattern)-1:
			# If we're at the end then stop
			break
		
		i+=1
		c = pattern[i]

	return minlen


def make_spacer_snp_network(spacer_fasta_dict, snp_thresh, outdir):
	"""Compares spacers using BLASTn and makes a network of similar spacers
	
	Args:
	  spacer_fasta_dict (dict of dicts):
		For each CRISPR type (main dict keys), a dict of spacer ID and seqs.
		  {"CRISPR_type": {"spacer_ID": "spacer_sequence"}}
	  snp_thresh (int):
		Number of SNPs between spacer, below which an edge should be
		made in the network
	  outdir (str):
		path to a directoy into which temporary files can be written

	Returns:
	  spacer_network_dict (dict):
	    A dict of each CRISPR type and the edges in the network of arrays
	    of that type.
	      {CR_type: [(node1, node2), (node1, node3) ... ]}
	"""
	
	# Check the outdir exists to write temp files to.
	if not os.path.isdir(outdir):
		sys.stderr.write(
			"Directory {} not found. Making it now.\n\n".format(outdir))
		os.makedirs(outdir)

	spacer_network_dict = defaultdict(list)

	for CR_type in spacer_fasta_dict:
		ID_dict = {}
		fasta_spacers = ""
		
		# If too many spacers, write them to a file and blast that
		if len(spacer_fasta_dict[CR_type]) > 250:
			
			# make file to store spacers to blast.
			# Make sure to create a file that doesn't exist yet...
			temp_file = outdir+"temp_spacers_0.fna"
			file_num = 0
			while os.path.isfile(temp_file):
				file_num +=1
				temp_file = outdir+"temp_spacers_{}.fna".format(file_num)
			
			outcontents = ""
			for n, sp in enumerate(list(spacer_fasta_dict[CR_type])):
				outcontents += ">{}\n{}\n".format(n,sp)
				ID_dict[str(n)] = sp
			with open(temp_file, "w") as fout:
				fout.write(outcontents)
			blastn_command = ("blastn -query {} "
				"-subject {} -task blastn-short "
				"-outfmt '6 std qlen slen' | "
				"awk '$1 != $2'"
				"".format(
					temp_file,
					temp_file))


		else:
			for n, sp in enumerate(list(spacer_fasta_dict[CR_type])):
				fasta_spacers += ">{}\\n{}\\n".format(n,sp)
				ID_dict[str(n)] = sp

			blastn_command = (f"blastn -query <(printf '{fasta_spacers}') "
				f"-subject <(printf '{fasta_spacers}') -task blastn-short "
				"-outfmt '6 std qlen slen' | "
				"awk '$1 != $2'")

		blast_run = subprocess.run(
			blastn_command,
			shell=True,
			universal_newlines=True,
			capture_output=True,
			executable='/bin/bash'
			)

		if blast_run.stderr:
			sys.stderr.write("ERROR running blast to dereplicate spacers:\n\n")
			sys.stderr.write(blast_run.stderr)
			sys.exit()
		blast_lines = [
			file_handling.BlastResult(
				i) for i in blast_run.stdout.split('\n') if len(i) > 0]

		for hit in blast_lines:
			if hit.mismatch > snp_thresh:
				continue
			if hit.length != hit.qlen or hit.length != hit.slen:
				dist = min(hamming(ID_dict[hit.qseqid], ID_dict[hit.sseqid]))
				if dist > snp_thresh:
					continue
			spacer_network_dict[CR_type].append(
				(ID_dict[hit.qseqid], ID_dict[hit.sseqid])
				)
		if len(spacer_fasta_dict[CR_type]) > 250:
			os.remove(temp_file)

	return spacer_network_dict


def identify_network_clusters(network):
	""" Make a list of clusters
	Each cluster is a dict containing the nodes (keys) and the number
	of edges connected to that node (values). Number of edges info is
	used to determine if the cluster is completely connected.
	
	Arg:
	  network (list of tuples):
	    list of edges in network where each element in the list is a
	    tuple of nodes joined by the edge

	Returns:
	  A list of the clusters identified in the network. Each element in
	  the list is a dict of connected nodes (keys) and the number of
	  edges connected to them (values). All nodes in a given dict have
	  a path to all other nodes in the same dict, but have no path to
	  any nodes in another dict in the list. However, the path could
	  have any length (i.e., clusters do not need to be completely
	  connected)

	  e.g.,
	  A network with the edges

	  A B
	  A C
	  D E

	  would be represented as
	  [
	  {"A":2, "B":1, "C":1},
	  {"D":1, "E":1}
	  ]
	"""

	cluster_list = []
	for edge in network:
		a, b = edge
		if len(cluster_list) == 0:
			cluster_list.append({a:1,b:1})
			continue

		a_idx = -1
		b_idx = -1
		for n, subdict in enumerate(cluster_list):
			if a in subdict:
				a_idx = n
			if b in subdict:
				b_idx = n

			if a_idx != -1 and b_idx != -1:
				# Both nodes found in existing clusters. Join clusters
				# if not already in the same cluster
				if a_idx == b_idx:
					cluster_list[a_idx][a] += 1
					cluster_list[a_idx][b] += 1
					break


				for k,v in cluster_list[b_idx].items():
					cluster_list[a_idx][k] = v
				cluster_list[a_idx][a] += 1
				cluster_list[a_idx][b] += 1

				del cluster_list[b_idx]
				break

			if subdict == cluster_list[-1]: 
				# If this is the last subdict and a or b is unfound
				if a_idx == -1 and b_idx == -1:
					cluster_list.append({a:1, b:1})
				
				elif a_idx == -1 and b_idx != -1:
					cluster_list[b_idx][a] = 1
					cluster_list[b_idx][b] += 1

				
				elif a_idx != -1 and b_idx == -1:
					cluster_list[a_idx][a] += 1
					cluster_list[a_idx][b] = 1

				break

	return cluster_list


def pick_cluster_rep(clusters, all_spacers, prev_spacer_id_dict):
	""" Pick a representative from each cluster
	
	The cluster member present in the most assemblies is taken to be the
	most likely original spacer sequence while other cluster members
	have SNPs relative to that.

	N.B. Using this process, some spacers in a cluster may differ by
	more than the specified number of SNPs if they are connected via
	an intermediate spacer with which they differ by fewer than the
	specified SNPs.
	"""

	spacer_count = Counter(all_spacers)

	rep_dict = {}

	for clus in clusters:
		max_count = 0
		for sp in clus:
			# prioiritize spacers that were in previous run if appending.
			if sp in prev_spacer_id_dict:
				rep = sp
				break
			if spacer_count[sp] > max_count:
				rep = sp
				max_count = spacer_count[sp]

		rep_dict[rep] = [sp for sp in clus if sp != rep]

	return rep_dict

