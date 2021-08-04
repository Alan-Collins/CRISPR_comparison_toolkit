#!/usr/bin/env python3

# AUTHOR	  :  ALAN COLLINS
# VERSION	 :  v1
# DATE		:  2021-8-3
# DESCRIPTION :  Process BLAST output of spacers against a blastdb. For results that have cut off due to mismatches, extend the hit to the full length and report mismatches. Report up- and down-stream bases for PAM analysis.

import sys
import argparse
import subprocess

class blast_result():
	"""
	A class to store column contents from a blast result in an easily retrieved form.
	Allows improved code readability when interacting with blast result lines

	Attributes:
		qseqid (str):   query (e.g., unknown gene) sequence id
		sseqid (str):   subject (e.g., reference genome) sequence id
		pident (float): percentage of identical matches
		length (int):   alignment length (sequence overlap)
		mismatch (int): number of mismatches
		gapopen (int):  number of gap openings
		qstart (int):   start of alignment in query
		qend (int): end of alignment in query
		sstart (int):   start of alignment in subject
		send (int): end of alignment in subject
		evalue (str):   expect value
		bitscore (float):   bit score
		qlen (int):  length of query sequence
		slen (int):  length of subject sequence
		strand (str):	Whether the blast hit was on the top (plus) or bottom (minus) strand of the DNA
	"""
	def __init__(self, blast_line):
		bits = blast_line.split('\t')
		self.qseqid = bits[0]
		self.sseqid = bits[1]
		self.pident = float(bits[2])
		self.length = int(bits[3])
		self.mismatch = int(bits[4])
		self.gapopen = bits[5]
		self.qstart = int(bits[6])
		self.qend = int(bits[7])
		self.sstart = int(bits[8])
		self.send = int(bits[9])
		self.evalue = bits[10]
		self.bitscore = bits[11]
		self.qlen = int(bits[12])
		self.slen = int(bits[13])
		self.strand = 'plus' if self.sstart < self.send else 'minus'


class protospacer():
	"""
	Class to store information about predicted protospacers.
	
	Attributes:
		up5 (str): The 5 bases upstream of the protospacer on the targeted strand.
		down5 (str): The 5 bases downstream of the protospacer on the targeted strand.
		protoseq (str): The sequence of the protospacer.
		pid (float): Percent identity between spacer and protospacer over whole length.
		mismatch (int): The number of mismatch positions between the spacer and protospacer.
		strand (str): ("plus"|"minus") orientation of protospacer relative to the strand represented in your fasta sequence.
		start (int): Index of first base in protospacer in the target sequence.
		stop  (int): Index of last base in protospacer in the target sequence.
		target (str): ID of the target sequence.

	"""
	def __init__(self, up5="", down5="", protoseq="", pid=0., mismatch=0, strand = "", start=0, stop=0, target=""):
		self.up5 = up5
		self.down5 = down5
		self.protoseq = protoseq
		self.pid = pid
		self.mismatch = mismatch
		self.strand = strand
		self.start = start
		self.stop = stop
		self.target = target
		

def run_blastcmd(db, fstring, batch_locations):
	"""
	function to call blastdbcmd in a shell and process the output. Uses a batch query of the format provided in
	the blastdbcmd docs. e.g.
	printf "%s %s %s %s\\n%s %s %s\\n" 13626247 40-80 plus 30 14772189 1-10 minus | blastdbcmd -db GPIPE/9606/current/all_contig -entry_batch -
	
	Args:
		db (str): path to the blast db you want to query.
		fstring (str):  The string to give to printf
		batch_locations (str):  the seqid, locations, and strand of all the spacers to be retrieved
	
	Returns:
		(list) List of the sequences of regions returned by blastdbcmd in response to the query locations submitted
	
	Raises:
		ERROR running blastdbcmd: Raises an exception when blastdbcmd returns something to stderr, 
		prints the error as well as the information provided to blastdbcmd and aborts the process.
	"""
	x = subprocess.run('printf "{}" {}| blastdbcmd -db {} -entry_batch -'.format(fstring, batch_locations, db), 
						shell=True, universal_newlines = True, capture_output=True) 
	if x.stderr:
		print("ERROR running blastdbcmd on {} :\n{}".format(db, batch_locations, x.stderr))
		sys.exit()
	else:
		return [i for i in x.stdout.split('\n') if '>' not in i and len(i) > 0]



def run_blastn(args):
	"""
	Runs blastn of provided spacers against provided blastdb. Processes the blast output to find hits and if the hits don't cover the whole length, the hits are extended using the provided blastdb.
	Args:
		args (argparse class): All of the argparse options given by the user.

	
	Returns:
		(list):  List of lines of a blast output file.
	
	Raises:
		ERROR running blast.
	"""	
	blastn_command = "blastn -query {} -db {} -task blastn-short -outfmt '6 std qlen slen' -num_threads {} -max_target_seqs {} -evalue {} {}".format(args.spacer_file, args.blast_db_path, args.num_threads, args.max_target_seqs, args.evalue, args.other_blast_options)
	blast_run = subprocess.run(blastn_command, shell=True, universal_newlines = True, capture_output=True)
	if blast_run.stderr:
		print("ERROR running blast on {}:\n{}".format(args.blast_db_path, blast_run.stderr))
		sys.exit()
	blast_lines = [blast_result(i) for i in blast_run.stdout.split('\n') if len(i) > 0]
	return blast_lines


def rev_comp(string):
	"""Converts As, Ts, Cs, and Gs, to their reverse complement nucleotide. Skips any characters that are not A, T, C, or G. Can work with upper or lower case sequence and returns the same case that was provided.
	Args:
		string (str): The sequence for which you want the reverse complement.
	
	Returns:
		(str) The reverse complement of the input sequence..
	"""
	rev_str = ''
	rev_comp_lookup = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c"}
	for i in reversed(string):
		if i in "ATCGatcg":
			rev_str += rev_comp_lookup[i]
		else:
			rev_str += i
	return rev_str


def hamming(seq1, seq2):
	"""
	Args:
		seq1 (str): First sequence to compare.
		seq2 (str): Second sequence to compare.

	Returns:
		(int) count of positions that differ between the two sequences.

	Raises:
		ValueError: Sequences must be the same length.
	"""

	if len(seq1) != len(seq2):
		raise ValueError("Sequences must be the same length.")

	diff = 0
	for a,b in zip(seq1, seq2):
		if a != b:
			diff+=1

	return diff


def fasta_to_dict(FASTA_file):
	""" Given a file in fasta format, opens the file and reads the headers and sequences into a dict.
	Args:
		FASTA_file (str): Path to the fasta file.
	
	Returns:
		(dict) Fasta sequences as values associated with their headers as keys
	"""
	fasta_dict = {}
	with open(FASTA_file, 'r') as f:
		multifasta = f.read()
	f.close()
	fastas = multifasta.split(">")
	trimmed_fastas = []
	for i in fastas:
		if len(i) != 0:
			trimmed_fastas.append(i)

	fastas = trimmed_fastas

	for i in fastas:
		header = i.split("\n")[0].split()[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def fill_info(db, spacer, result):
	"""
	Args:
		result (blast_result class): The blast_result instance to process.
	
	Returns:
		(protospacer class) This is a description of what is returned.
	"""
	# {'qseqid': 'C1S1', 'sseqid': 'DMS3', 'pident': 100.0, 'length': 13, 'mismatch': 0, 'gapopen': '0', 'qstart': 19, 'qend': 31, 'sstart': 31243, 'send': 31231, 'evalue': '0.010', 'bitscore': '26.3', 'qlen': 32, 'slen': 36415, 'strand': 'minus', 'trunc': False, 'sstart_mod': False, 'send_mod': False}
	
	proto = protospacer(target=result.sseqid, strand=result.strand)

	# Adjust for missing edges due to blast not extending through mismatches. If there are no missing edges this won't add anything.
	if proto.strand == 'plus':
		proto.start = result.sstart - (result.qstart - 1) # Adjust for 1 base counting
		proto.stop = result.send + (result.qlen - result.qend)
		up5_coords = (proto.start - 5, proto.start - 1)
		down5_coords = (proto.stop + 1, proto.stop + 5)
	else:
		proto.start = result.send - (result.qlen - result.qend)
		proto.stop = result.sstart + (result.qstart -1) # Adjust for 1 base counting
		up5_coords = (proto.stop + 1, proto.stop + 5)
		down5_coords = (proto.start - 5, proto.start - 1)	

	# Build blastdbcmd inputs
	fstring = ""
	batch_locations = ""

	for loc in [up5_coords, (proto.start, proto.stop), down5_coords]:
		fstring += '%s %s %s\n'
		batch_locations += '{} {}-{} {} '.format(proto.target, loc[0], loc[1], proto.strand)

	proto.up5, proto.protoseq, proto.down5 = run_blastcmd(db, fstring, batch_locations)
	proto.mismatch = hamming(spacer, proto.protoseq)
	proto.pid = 100-(100*proto.mismatch/result.qlen)

	return proto



def main():

	parser = argparse.ArgumentParser(
		description="Process BLAST output of spacers against a blastdb. For results that have cut off due to mismatches, extend the hit to the full length and report mismatches. Report up- and down-stream bases for PAM analysis."
		)
	parser.add_argument(
		"-d", dest="blast_db_path", required = True,
		help="path to blast db files (not including the file extensions). The blastdb must have been made with the option '-parse_seqids' for this script to function."
		)
	parser.add_argument(
		"-s", dest="spacer_file", required = True,
		help="The file with your spacers in fasta format."
		)
	parser.add_argument(
		"-o", dest="outfile", required = True,
		help="path to output file."
		)
	parser.add_argument(
		"-e", dest="evalue", required = False, default='10',
		help="DEFAULT: 10. set the evalue cutoff below which blastn will keep blast hits when looking for CRISPR repeats in your blast database. Useful for reducing inclusion of low quality blast hits with big databases in combination with the -m option."
		)
	parser.add_argument(
		"-m", dest="max_target_seqs", required = False, default='10000',
		help="DEFAULT: 10000. Set the max_target_seqs option for blastn when looking for CRISPR repeats in your blast database. Blast stops looking for hits after finding and internal limit (N_i) sequences for each query sequence, where N_i=2*N+50. These are just the first N_i sequences with better evalue scores than the cutoff, not the best N_i hits. Because of the nature of the blast used here (small number of queries with many expected hits) it may be necessary to increase the max_target_seqs value to avoid blast ceasing to search for repeats before all have been found. The blast default value is 500. The default used here is 10,000. You may want to reduce it to increase speed or increase it to make sure every repeat is being found. If increasing this value (e.g. doubling it) finds no new spacers then you can be confident that this is not an issue with your dataset."
		)
	parser.add_argument(
		"-t", dest="num_threads", required = False, default=1, type=int,
		help="DEFAULT: 1. Number of threads you want to use for the blastn step of this script."
		)
	parser.add_argument(
		"-x", dest="other_blast_options", required = False, default='',
		help="DEFAULT: none. If you want to include any other options to control the blastn command, you can add them here. Options you should not provide here are: blastn -query -db -task -outfmt -num_threads -max_target_seqs -evalue"
		)


	args = parser.parse_args(sys.argv[1:])

	spacer_dict = fasta_to_dict(args.spacer_file)

	blast_output = run_blastn(args)

	test = blast_output[0]

	first = fill_info(args.blast_db_path, spacer_dict[test.qseqid], test)

	print(vars(first))


if __name__ == '__main__':
	main()

