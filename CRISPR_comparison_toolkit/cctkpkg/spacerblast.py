#!/usr/bin/env python3


import sys
import argparse
import subprocess
import re
from collections import defaultdict

class BlastResult():
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


class Protospacer():
	"""
	Class to store information about predicted protospacers.
	
	Attributes:
		spacer (str): The fasta header of the spacer.
		up (str): The bases upstream of the protospacer on the targeted strand.
		down (str): The bases downstream of the protospacer on the targeted strand.
		protoseq (str): The sequence of the protospacer.
		pid (float): Percent identity between spacer and protospacer over whole length.
		mismatch (int): The number of mismatch positions between the spacer and protospacer.
		strand (str): ("plus"|"minus") orientation of protospacer relative to the strand represented in your fasta sequence.
		start (int): Index of first base in protospacer in the target sequence.
		stop  (int): Index of last base in protospacer in the target sequence.
		target (str): ID of the target sequence.
		length (int): The length of the spacer.
	"""
	def __init__(self, spacer="", up="", down="", protoseq="", pid=0., mismatch=0, strand = "", start=0, stop=0, target="", length=0):
		self.spacer = spacer
		self.up = up
		self.down = down
		self.protoseq = protoseq
		self.pid = pid
		self.mismatch = mismatch
		self.strand = strand
		self.start = start
		self.stop = stop
		self.target = target
		self.length = length
		

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
	blast_lines = [BlastResult(i) for i in blast_run.stdout.split('\n') if len(i) > 0]
	return blast_lines


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


def fill_initial_info(result, flanking_n):
	"""
	Args:
		result (BlastResult class): The BlastResult instance to process.
		flanking_n (tuple): How many bases upstream and downstream do you want?
	
	Returns:
		tuple(Protospacer class, str, str) Protospacer class instance with some info added and the fstring and batch_locations for the protospacer and flanking sequence
	"""

	if int(result.gapopen) > 0: # If there are gaps between spacer and protospacer then ignore.
		return None, None, None
	
	proto = Protospacer(spacer=result.qseqid,target=result.sseqid, strand=result.strand, length=result.qlen)
	up_n, down_n = flanking_n
	# Adjust for missing edges due to blast not extending through mismatches. If there are no missing edges this won't add anything.
	if proto.strand == 'plus':
		proto.start = result.sstart - (result.qstart - 1) # Adjust for 1 base counting
		proto.stop = result.send + (result.qlen - result.qend)
		up_coords = (proto.start - up_n, proto.start - 1)
		down_coords = (proto.stop + 1, proto.stop + down_n)

	else:
		proto.start = result.send - (result.qlen - result.qend)
		proto.stop = result.sstart + (result.qstart -1) # Adjust for 1 base counting
		up_coords = (proto.stop + 1, proto.stop + up_n)
		down_coords = (proto.start - down_n, proto.start - 1)	

	if any([i < 1 for i in [up_coords[0], up_coords[1], down_coords[0], down_coords[1]]]) or any([i > result.slen-1 for i in [up_coords[0], up_coords[1], down_coords[0], down_coords[1]]]): # If the requested flanking region goes beyond the ends of the target sequence then it can't be returned.
		return None, None, None

	# Build blastdbcmd inputs
	fstring = ""
	batch_locations = ""

	if up_n and down_n:
		coords_list = [up_coords, (proto.start, proto.stop), down_coords]
	elif up_n:
		coords_list = [up_coords, (proto.start, proto.stop)]
	elif down_n:
		coords_list = [(proto.start, proto.stop), down_coords]
	else:
		coords_list = [(proto.start, proto.stop)]
	
	for loc in coords_list:
		fstring += '%s %s %s\n'
		batch_locations += '{} {}-{} {} '.format(proto.target, loc[0], loc[1], proto.strand)

	return proto, fstring, batch_locations


def fill_remaining_info(proto, spacer, blastdbcmd_output):
	"""
	Args:
		proto (Protospacer class) Protospacer instance containing initial information that you want to fill.
		spacer (str): sequence of the spacer.
		blastdbcmd_output (list): elements of the blastdbcmd output for this protospacer.
	
	Returns:
		(Protospacer class) Protospacer instance with the rest of the information filled in.
	"""

	proto.up, proto.protoseq, proto.down = blastdbcmd_output 
	proto.mismatch = hamming(spacer, proto.protoseq)
	proto.pid = 100-(100*proto.mismatch/proto.length)

	return proto


def find_pam(pattern, sequence, is_regex):
	"""Search the provided sequence for the provided pattern. Use regex is specified.
	Args:
		pattern (str): The pattern describing your PAM.
		sequence (str): The sequence flanking your protospacer in which you expect the PAM to be located.
		is_regex (bool): Is your pattern a regex pattern?
	
	Returns:
		(Bool) Was a PAM found?.
	
	Raises:
		ValueError: Your regex could not be parsed.
		ValueError: Your pattern and sequence are different lengths.
		ValueError: Your non-regex pattern can only contain the following characters: N, A, T, C, G.
	"""
	if is_regex:
		try:
			pat=re.compile(pattern, re.IGNORECASE)
		except Exception as e:
			raise ValueError("Your regex could not be parsed. Please check the following problem: {}".format(e))
	
	else:
		if len(pattern) != len(sequence):
			raise ValueError("Your pattern and sequence are different lengths. Make sure your pattern is the same number of bases as are retrieved from the flanking region of your protospacer. You provided a pattern of {} bp and a flanking sequence of {} bp".format(len(pattern), len(sequence)))
		
		if any([i.upper() not in ["N", "A", "T", "C", "G"] for i in pattern]):
			raise ValueError("Your non-regex pattern can only contain the following characters: N, A, T, C, G.")

		pattern = pattern.upper().replace("N", ".")
		pat=re.compile(pattern, re.IGNORECASE)

	return bool(pat.match(sequence))


def build_parser(parser):
	req_options = parser.add_argument_group("Required arguments")
	req_options.add_argument(
		"-d", "--blastdb",
		dest="blast_db_path",
		required = True,
		help="path to blast db files (not including the file extensions). The \
		blastdb must have been made with the option '-parse_seqids' for this \
		script to function."
		)
	req_options.add_argument(
		"-s", "--spacers",
		dest="spacer_file",
		required = True,
		help="The file with your spacers in fasta format."
		)
	
	out_options = parser.add_argument_group("Output control arguments")
	out_options.add_argument(
		"-o", "--out",
		metavar=' ',
		dest="outfile",
		required = False,
		help="path to output file. If none provided, outputs to stdout."
		)
	out_options.add_argument(
		"-q", "--no-pam-out",
		metavar=' ',
		dest="no_pam_outfile",
		required = False,
		help="path to output file for protospacers with no PAM if desired."
		)
	out_options.add_argument(
		"-n", "--flanking",
		metavar=' ',
		dest="flanking_n",
		required = False,
		default=0,
		type=int,
		help="DEFAULT: 0. Number of bases you want returned from the 5' and \
		3' sequences flanking your protospacers. N.B. In order to \
		consistently return these sequences, protospacers that are closer to \
		the end of the target sequence than the number specified here will be \
		discarded."
		)
	out_options.add_argument(
		"-u", "--upstream",
		metavar=' ',
		dest="upstream_n",
		required = False,
		default=0,
		type=int,
		help="DEFAULT: 0. Number of bases you want returned from the 5' \
		sequences flanking your protospacers. N.B. In order to consistently \
		return these sequences, protospacers that are closer to the end of \
		the target sequence than the number specified here will be discarded."
		)
	out_options.add_argument(
		"-w", "--downstream",
		metavar=' ',
		dest="downstream_n", 
		required = False,
		default=0,
		type=int,
		help="DEFAULT: 0. Number of bases you want returned from the 3' \
		sequences flanking your protospacers. N.B. In order to consistently \
		return these sequences, protospacers that are closer to the end of \
		the target sequence than the number specified here will be discarded."
		)
	out_options.add_argument(
		"-R", "--regex-pam",
		metavar=' ',
		required = False,
		help='A regex describing the PAM sequence you would like to look for. \
		E.g. if your PAM is C or T then two Gs, you could describe that as \
		"[CT]GG"'
		)
	out_options.add_argument(
		"-P", "--pam",
		metavar=' ',
		required = False,
		help="A pattern describing the PAM sequence you would like to look \
		for. Ns can be used to indicate positions with no sequence \
		requirement. Case insensitive. E.g. if your PAM is 3 of any base \
		followed by CC you could provide NNNCC or nnncc here."
		)
	out_options.add_argument(
		"-l", "--pam-location",
		metavar=' ',
		required = False,
		choices=['up', 'down'],
		help="{'up', 'down'} Where is your PAM relative to your protospacer? \
		'up' indicates 5' and down indicates 3'."
		)
	out_options.add_argument(
		"-p", "--percent-id",
		metavar=' ',
		dest="pid",
		required = False,
		default=0,
		type=float,
		help="DEFAULT: 0. Minimum percent identity between spacer and \
		protospacer in order for the result to be output. Default is output \
		all protospacers that are returned by BLAST."
		)
	out_options.add_argument(
		"-b", "--batch-size",
		metavar=' ',
		dest="blastdbcmd_batch_size",
		required = False,
		default=1000,
		type=int,
		help="DEFAULT: 1000. This runs quicker if it calls blastdbcmd fewer \
		times. To get the sequences of protospacers and flanking sequence, \
		locations are retrieved from blastdbcmd in batches. The larger the \
		batch, the quicker this runs. However, your OS may have a limit on \
		the number that can be used. If you get an error like 'OSError: \
		[Errno 7] Argument list too long: '/bin/sh'' then decrease this value \
		and try again. In my experience, 500 often works."
		)
	out_options.add_argument(
		"-r", "--mask-regions",
		metavar=' ',
		required = False,
		help="If you would like to mask regions in your blastdb to exclude \
		hits in those regions, provide a .bed format file with those regions \
		here. This is useful for example if you know that the sequences in \
		your blastdb contain CRISPR arrays and would like to exclude hits of \
		CRISPR spacers against similar spacers in the array."
		)

	### Options for blastn ####

	blast_options = parser.add_argument_group(
		"BLAST arguments", 
		"Arguments to control the blastn command used by this script.")
	blast_options.add_argument(
		"-e", "--evalue",
		metavar=' ',
		required = False,
		default='10',
		help="DEFAULT: 10. set the evalue cutoff below which blastn will keep \
		blast hits when looking for CRISPR repeats in your blast database. \
		Useful for reducing inclusion of low quality blast hits with big \
		databases in combination with the -m option."
		)
	blast_options.add_argument(
		"-m", "--max-target-seqs",
		metavar=' ',
		required = False,
		default='10000',
		help="DEFAULT: 10000. Set the max_target_seqs option for blastn when \
		looking for CRISPR repeats in your blast database. Blast stops \
		looking for hits after finding and internal limit (N_i) sequences for \
		each query sequence, where N_i=2*N+50. These are just the first N_i \
		sequences with better evalue scores than the cutoff, not the best N_i \
		hits. Because of the nature of the blast used here (small number of \
		queries with many expected hits) it may be necessary to increase the \
		max_target_seqs value to avoid blast ceasing to search for repeats \
		before all have been found. The blast default value is 500. The \
		default used here is 10,000. You may want to reduce it to increase \
		speed or increase it to make sure every repeat is being found. If \
		increasing this value (e.g. doubling it) finds no new spacers then \
		you can be confident that this is not an issue with your dataset."
		)
	blast_options.add_argument(
		"-t", "--threads",
		dest="num_threads", 
		metavar=' ',
		required = False,
		default=1,
		type=int,
		help="DEFAULT: 1. Number of threads you want to use for the blastn \
		step of this script."
		)
	blast_options.add_argument(
		"-x", "--other-options",
		dest="other_blast_options",
		metavar=' ',
		required = False,
		default='',
		help="DEFAULT: none. If you want to include any other options to \
		control the blastn command, you can add them here. Options you should \
		not provide here are: blastn -query -db -task -outfmt -num_threads \
		-max_target_seqs -evalue"
		)

	return parser


def main(args):

	if any([args.pam_location, args.pam, args.regex_pam]):
		if all([args.pam, args.regex_pam]):
			print("Please provide either a PAM pattern using -P/--pam or a \
				regex using -R/--regex_pam but not both.")
			sys.exit()
		if not all(
			[args.pam_location, args.pam]) and not all(
			[args.pam_location, args.regex_pam]):
			print("Please provide both a PAM pattern using -P/--pam or \
				-R/--regex_pam AND a location to search for the PAM using \
				-l/--pam_location.")
			sys.exit()

	spacer_dict = fasta_to_dict(args.spacer_file)
	if args.flanking_n:
		flanking_n = (args.flanking_n, args.flanking_n)
	else:
		flanking_n = (args.upstream_n, args.downstream_n)

	mask_dict = defaultdict(list)
	if args.mask_regions:
		with open(args.mask_regions, 'r') as fin:
			for line in fin.readlines():
				if "#" in line:
					continue
				bits = line.split()
				if len(bits) < 3:
					print("Error: Number of columns in your mask_regions bed \
						file must be at least 3:\nFasta_header\tstart\tstop")
					sys.exit()
				mask_dict[bits[0]].append([int(i) for i in bits[1:3]])


	blast_output = run_blastn(args)

	outcontents = ["Spacer_ID\tTarget_contig\tProtospacer_start\t\
	Protospacer_end\tPercent_identity\tmismatches\tprotospacer_sequence\t\
	upstream_bases\tdownstream_bases\ttarget_strand"]
	no_pam_outcontents = ["Spacer_ID\tTarget_contig\tProtospacer_start\t\
	Protospacer_end\tPercent_identity\tmismatches\tprotospacer_sequence\t\
	upstream_bases\tdownstream_bases\ttarget_strand"]
	fstring = batch_locations = ""
	protos = []

	count = 0
	all_protospacer_infos = []
	for result in blast_output:
		p, f, b = fill_initial_info(result, flanking_n)
		# If the match couldn't be extended because it is at the end of
		# the contig then skip this one.
		if p == None:
			continue
		protos.append(p)
		fstring += f
		batch_locations += b
		count += 1

		if count == args.blastdbcmd_batch_size:
			count = 0
			all_protospacer_infos += run_blastcmd(
				args.blast_db_path,
				fstring,
				batch_locations)
			fstring = batch_locations = ""
	if count != 0:
		all_protospacer_infos += run_blastcmd(
			args.blast_db_path,
			fstring,
			batch_locations)

	p_count = 0

	# Want to work through the 
	increment = 1 + int(bool(flanking_n[0])) + int(bool(flanking_n[1]))
	for i in range(0,len(all_protospacer_infos),increment):
		p = protos[p_count]
		# Depending on which flanking sequence was requested, fill in unrequested sequence with empty string.

		if flanking_n[0] and flanking_n[1]:
			p = fill_remaining_info(
				p,
				spacer_dict[p.spacer],
			all_protospacer_infos[i:i+3])
		elif flanking_n[0]:
			p = fill_remaining_info(
				p,
				spacer_dict[p.spacer],
				all_protospacer_infos[i:i+2]+[''])
		elif flanking_n[1]:
			p = fill_remaining_info(
				p,
				spacer_dict[p.spacer],
				['']+all_protospacer_infos[i:i+2])
		else:
			p = fill_remaining_info(
				p,
				spacer_dict[p.spacer],
				['', all_protospacer_infos[i], ''])

		if p.pid >= args.pid:
			if p.target in mask_dict.keys():
				mask_region = []
				for region in mask_dict[p.target]:
					start, stop = region
					if start > stop:
						start, stop = stop, start
					mask_region += [i for i in range(start, stop)]
				if p.start in mask_region or p.stop in mask_region:
					p_count += 1
					continue
			# Check for PAM if user requested
			if args.pam or args.regex_pam:
				if args.pam_location == 'up':
						seq = p.up
				else:
					seq = p.down
				if args.pam:
					if not find_pam(args.pam, seq, False):
						no_pam_outcontents.append(
							"\t".join([str(_) for _ in [
								p.spacer,
								p.target,
								p.start,
								p.stop,
								p.pid,
								p.mismatch,
								p.protoseq,
								p.up,
								p.down,
								p.strand]]))
						p_count += 1
						continue
				else:
					if not find_pam(args.regex_pam, seq, True):
						no_pam_outcontents.append(
							"\t".join([str(_) for _ in [
								p.spacer,
								p.target,
								p.start,
								p.stop,
								p.pid,
								p.mismatch,
								p.protoseq,
								p.up,
								p.down,
								p.strand]]))
						p_count += 1
						continue
			outcontents.append(
				"\t".join([str(_) for _ in [
					p.spacer,
					p.target,
					p.start,
					p.stop,
					p.pid,
					p.mismatch,
					p.protoseq,
					p.up,
					p.down,
					p.strand]]))

		p_count += 1

	if args.outfile:
		with open(args.outfile, 'w') as fout:
			fout.write("\n".join(outcontents) + "\n")
	else:
		for line in outcontents:
			try: # Prevents Broken pipe error if the user pipes output
				print(line)
			except:
				sys.exit()

	if args.no_pam_outfile:
		with open(args.no_pam_outfile, 'w') as fout:
			fout.write("\n".join(no_pam_outcontents) + "\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Process BLAST output of spacers against a blastdb. For \
		results that have cut off due to mismatches, extend the hit to the \
		full length and report mismatches. Report up- and down-stream bases \
		for PAM analysis.")
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)
