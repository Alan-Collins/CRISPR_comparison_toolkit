#!/usr/bin/env python3


import sys
import argparse
import subprocess
import re
from copy import copy
from collections import defaultdict

from . import sequence_operations, file_handling

description = """
usage: cctk spacerblast [-h] -d <path (no extension)> -s <path> [-o <path>] \
[-q <path>] [-n <int>] [-u <int>] [-w <int>] [-R <string>] [-P <string>] \
[-l {'up', 'down'}] [-p <float>] [-b <int>] [-r <path>] [-e <float>] \
[-m <int>] [-t <int>] [-x <string>]

optional arguments:
  -h, --help        show this help message and exit

required arguments:
  -d, --blastdb     path to blast db excluding extension
  -s, --spacers     spacers in fasta format

output control arguments:
  -o, --out         path to output file. Default: stdout
  -q, --no-pam-out  path to output file for protospacers with no PAM, if desired
  -n, --flanking    number of bases to return from both sides of protospacers
  -u, --upstream    number of bases to return from the 5'side of protospacers
  -w, --downstream  number of bases to return from the 3'side of protospacers
  -R, --regex-pam   regex describing the PAM sequence
  -P, --pam         nucleotide pattern describing the PAM sequence e.g. NNNCC
  -l, --pam-location
                    {'up', 'down'} PAM location relative to protospacer
  -p, --percent-id  minimum percent identity between spacer and protospacer
  -b, --batch-size  size of batch to submit to blastdbcmd. Default: 1000
  -r, --mask-regions
                    file in bed format listing regions to ignore
  -E, --seed-region
                    Specify part of protospacer in which mismatches should not be 
                    tolerated. Format: start:stop, 0-base coordinates, 5':3'. E.g., "0:6" 
                    or ":6" specifies first 6 bases (0,1,2,3,4,5). "-6:-1" or "-6:" 
                    specifies last 6 bases.

BLAST arguments:
  arguments to control the blastn command used by this script

  -e, --evalue      blast e value. Default: 10
  -m, --max-target-seqs
                    blast max target seqs. Default: 10000
  -t, --threads     threads to use. Default: 1
  -x, --other-options
                    additional blastn options. Forbidden options: \
blastn -query -db -task -outfmt -num_threads -max_target_seqs -evalue
"""

class Protospacer():
	"""
	Class to store information about predicted protospacers.
	
	Attributes:
		spacer (str): The fasta header of the spacer.
		up (str): The bases upstream of the protospacer on the targeted
		  strand.
		down (str): The bases downstream of the protospacer on the
		  targeted strand.
		protoseq (str): The sequence of the protospacer.
		pid (float): Percent identity between spacer and protospacer
		  over whole length.
		mismatch (int): The number of mismatch positions between the
		  spacer and protospacer.
		strand (str): ("plus"|"minus") orientation of protospacer
		  relative to the strand represented in your fasta sequence.
		start (int): Index of first base in protospacer in the target
		  sequence.
		stop  (int): Index of last base in protospacer in the target
		  sequence.
		target (str): ID of the target sequence.
		length (int): The length of the spacer.
	"""
	def __init__(self,
		spacer="",
		up="",
		down="",
		protoseq="",
		aligned="",
		pid=0.,
		mismatch=0,
		strand="",
		start=0,
		stop=0,
		target="",
		length=0
		):
		self.spacer = spacer
		self.up = up
		self.down = down
		self.protoseq = protoseq
		self.aligned = aligned
		self.pid = pid
		self.mismatch = mismatch
		self.strand = strand
		self.start = start
		self.stop = stop
		self.target = target
		self.length = length
		

def run_blastn(args):
	"""
	Runs blastn of provided spacers against provided blastdb. Processes
	the blast output to find hits and if the hits don't cover the whole
	length, the hits are extended using the provided blastdb.
	Args:
		args (argparse class): All of the argparse options given by the
		  user.

	
	Returns:
		(list):  List of lines of a blast output file.
	
	Raises:
		ERROR running blast.
	"""	
	blastn_command = "blastn -query {} -db {} -task blastn-short \
	-outfmt '6 std qlen slen' -num_threads {} -max_target_seqs {} -evalue {} \
	{}".format(args.spacer_file,
		args.blast_db_path,
		args.num_threads,
		args.max_target_seqs,
		args.evalue,
		args.other_blast_options)
	blast_run = subprocess.run(blastn_command,
		shell=True,
		universal_newlines=True,
		capture_output=True)
	if blast_run.stderr:
		sys.stderr.write("ERROR running blast on {}:\n{}".format(
			args.blast_db_path, blast_run.stderr))
		sys.exit()
	blast_lines = [
		file_handling.BlastResult(
			i) for i in blast_run.stdout.split('\n') if len(i) > 0]
	return blast_lines


def fill_initial_info(result, flanking_n):
	"""
	Args:
		result (BlastResult class): The BlastResult instance to process.
		flanking_n (tuple): How many bases upstream and downstream do
		  you want?
	
	Returns:
		tuple(Protospacer class, str, str) Protospacer class instance
		  with some info added and the fstring and batch_locations for
		  the protospacer and flanking sequence
	"""

	if int(result.gapopen) > 0: 
	  # If there are gaps between spacer and protospacer then ignore.
		return None, None, None
	
	proto = Protospacer(
		spacer=result.qseqid,
		target=result.sseqid,
		strand=result.strand,
		length=result.qlen)
	up_n, down_n = flanking_n
	# Adjust for missing edges due to blast not extending through
	# mismatches. If there are no missing edges this won't add anything.
	if proto.strand == 'plus':
		proto.start = result.sstart - (result.qstart-1) # Adjust 1 base
		proto.stop = result.send + (result.qlen-result.qend)
		up_coords = (proto.start-up_n, proto.start-1)
		down_coords = (proto.stop+1, proto.stop+down_n)
		
		if (any(
			[i < 1 for i in [
				up_coords[0],
				up_coords[1],
				down_coords[0],
				down_coords[1]]]
			) and up_n != 0) or (any(
				[i > result.slen-1 for i in [
					up_coords[0],
					up_coords[1],
					down_coords[0],
					down_coords[1]]]
			) and down_n != 0):
			# If the requested flanking region goes beyond the ends of the
			# target sequence then it can't be returned.
			return None, None, None

	else:
		proto.start = result.send - (result.qlen-result.qend)
		proto.stop = result.sstart + (result.qstart-1) # Adjust 1 base
		up_coords = (proto.stop+1, proto.stop+up_n)
		down_coords = (proto.start-down_n, proto.start-1)	
		if (any(
			[i < 1 for i in [
				up_coords[0],
				up_coords[1],
				down_coords[0],
				down_coords[1]]]
			) and down_n != 0) or (any(
				[i > result.slen-1 for i in [
					up_coords[0],
					up_coords[1],
					down_coords[0],
					down_coords[1]]]
			) and up_n != 0):
				# If the requested flanking region goes beyond the ends of the
				# target sequence then it can't be returned.
				return None, None, None
	
	if any(
			[i < 1 for i in [
				up_coords[0],
				up_coords[1],
				down_coords[0],
				down_coords[1]]]
			) or any(
				[i > result.slen-1 for i in [
					up_coords[0],
					up_coords[1],
					down_coords[0],
					down_coords[1]]]
			):
		# requested region goes beyond ends of target sequence and
		# Cannot be returned
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
		batch_locations += '{} {}-{} {} '.format(
			proto.target,
			loc[0],
			loc[1],
			proto.strand)
	
	return proto, fstring, batch_locations


def fill_remaining_info(proto, spacer, blastdbcmd_output):
	"""
	Args:
		proto (Protospacer class) Protospacer instance containing
		  initial information that you want to fill.
		spacer (str): sequence of the spacer.
		blastdbcmd_output (list): elements of the blastdbcmd output for
		  this protospacer.
	
	Returns:
		(Protospacer class) Protospacer instance with the rest of the
		  information filled in.
	"""

	proto.up, proto.protoseq, proto.down = blastdbcmd_output 
	proto.mismatch, _, _ = sequence_operations.hamming(spacer, proto.protoseq)
	proto.aligned = "".join(
		[p if s==p else '.' for s,p in zip(spacer, proto.protoseq)])
	proto.pid = 100-(100*proto.mismatch/proto.length)

	return proto


def compile_pam(pattern, is_regex):
	"""Reformat non-regex pam into regex if needed and compile.
	Args:
		pattern (str): The pattern describing your PAM.
		is_regex (bool): Is your pattern a regex pattern?
	
	Returns:
		(re.compile) Compiled pattern.
	
	Raises:
		ValueError: Your regex could not be parsed.
		ValueError: Your non-regex pattern can only contain the
		  following characters: [list of characters].
	"""
	if is_regex:
		try:
			pat=re.compile(pattern, re.IGNORECASE)
		except Exception as e:
			raise ValueError("Your regex could not be parsed. Please check "
				"the following problem: {}".format(e))
	
	else:
		if any([
			i.upper() not in sequence_operations.SEQUENCE_DICT for i in pattern]):
			raise ValueError("Your non-regex pattern can only contain the "
				"following characters: {}".format(", ".join([
					i for i in sequence_operations.SEQUENCE_DICT.keys()
					])))

		pattern = "".join([
			sequence_operations.SEQUENCE_DICT[i] for i in pattern.upper()])

	pat = re.compile(pattern, re.IGNORECASE)

	return pat


def check_pam_length(args):
	if args.pam:
		minlen = len(args.pam)
	else:
		# First check if the user has given a match group:
		if "(" in args.regex_pam:
			print("Error: Cannot process match groups. Please do not use "
				"parentheses () in your regex.")
			sys.exit()
		minlen = sequence_operations.determine_regex_length(args.regex_pam)

	return minlen


def set_flanking_n(args, min_pam_len):
	if args.flanking_n:
		flanking_n = (args.flanking_n, args.flanking_n)
	else:
		flanking_n = (args.upstream_n, args.downstream_n)

	# Make sure flanking seq long enough to catch pam
	if args.pam_location == "up":
		if flanking_n[0] < min_pam_len:
			sys.stderr.write(
				"Your specified PAM is at least {} bases, but you only "
				"requested {} upstream bases. {} bases will now be "
				"retrieved on the upstream side.\n".format(
					min_pam_len, flanking_n[0], min_pam_len))
			flanking_n = (min_pam_len, flanking_n[1])
	if args.pam_location == "down":
		if flanking_n[1] < min_pam_len:
			sys.stderr.write(
				"Your specified PAM is at least {} bases, but you only "
				"requested {} downstream bases. {} bases will now be "
				"retrieved on the downstream side.\n".format(
					min_pam_len, flanking_n[1], min_pam_len))
			flanking_n = (flanking_n[0], min_pam_len)

	return flanking_n


def adjust_pam(pam, min_pam_len, pam_location, flanking_n):
	if pam_location == 'up':
		if min_pam_len == flanking_n[0]:
			return pam

		new_pattern = "[ATCG]+" + pam.pattern
		pam = re.compile(new_pattern, re.IGNORECASE)
		return pam

	else:
		if min_pam_len == flanking_n[1]:
			return pam

		new_pattern = pam.pattern + "[ATCG]+"
		pam = re.compile(new_pattern, re.IGNORECASE)
		return pam


def check_seed_arg(arg):
	if arg == None:
		return None
	
	if (
		(":" not in arg) 
		or (len(arg.split(":")) > 2)
		or any([i not in {"0","1","2","3","4","5","6","7","8","9",":", "-", " "} for i in arg])
	):
		sys.stderr.write('ERROR: Seed region must be specified as a range of start:stop with no spaces. '
		+ 'i.e., 0-base integer start coordinate, a ":" character, and a 0-base stop '
		+ 'coordinate. See the documentation for details.\n')
		sys.exit()
	
	start, stop = arg.split(":")
	# convert to ints or default if blank
	if start.strip() == '':
		start = 0
	else:
		start = int(start)
	
	if stop.strip() == '':
		stop = -1
	else:
		stop = int(stop)
	
	return (start, stop)


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
	out_options.add_argument(
		"-E", "--seed-region",
		metavar=' ',
		required = False,
		nargs="?",
		help='Specify part of protospacer in which mismatches should not be \
		tolerated. Format: start:stop 0-base coordinates 5\':3\'. E.g., "0:6" \
		or ":6" specifies first 6 bases (0,1,2,3,4,5). "-6:-1" or "-6:" \
		specifies last 6 bases.'
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

	# Check seed region arg format if specified
	seed_reg = check_seed_arg(args.seed_region)

	spacer_dict = file_handling.fasta_to_dict(args.spacer_file)

	min_pam_len = 0
	pam = None
	if any([args.pam_location, args.pam, args.regex_pam]):
		if all([args.pam, args.regex_pam]):
			print("Please provide either a PAM pattern using -P/--pam or a "
				"regex using -R/--regex_pam but not both.")
			sys.exit()
		if not all(
			[args.pam_location, args.pam]) and not all(
			[args.pam_location, args.regex_pam]):
			print("Please provide both a PAM pattern using -P/--pam or "
				"-R/--regex_pam AND a location to search for the PAM using "
				"-l/--pam_location.")
			sys.exit()

		min_pam_len = check_pam_length(args)

		if args.pam:
			pattern = args.pam
			is_regex = False

		else:
			pattern = args.regex_pam
			is_regex = True

		pam = compile_pam(pattern, is_regex)

	# Adjust flanking_n in case it is too few bases for PAM
	flanking_n = set_flanking_n(args, min_pam_len)

	if pam:
		# Adjust PAM with extra "[ATCG]+" if it is too short
		pam = adjust_pam(pam, min_pam_len, args.pam_location, flanking_n)

	mask_dict = defaultdict(list)
	if args.mask_regions:
		with open(args.mask_regions, 'r') as fin:
			for line in fin.readlines():
				if "#" in line:
					continue
				bits = line.split()
				if len(bits) < 3:
					print("Error: Number of columns in your mask_regions bed "
						"file must be at least 3:\nFasta_header\tstart\tstop")
					sys.exit()
				mask_dict[bits[0]].append([int(i) for i in bits[1:3]])


	blast_output = run_blastn(args)

	outcontents = ["Spacer_ID\tTarget_contig\tProtospacer_start\t\
Protospacer_end\tPercent_identity\tmismatches\tprotospacer_sequence\t\
mismatch_locations"]
	if flanking_n[0] != 0:
		outcontents[0] += "\tupstream_bases"
	if flanking_n[1] != 0:
		outcontents[0] += "\tdownstream_bases"
	outcontents[0] += "\ttarget_strand"

	no_pam_outcontents = copy(outcontents)
	fstring = batch_locations = ""
	protos = []

	count = 0
	pool_input = [] # grouped information to give to blastdbcmd
	# Want to work through the flanking sequences in increments that
	# depend upon which flanking sequences were retrieved
	increment = 1 + int(bool(flanking_n[0])) + int(bool(flanking_n[1]))
	args.blastdbcmd_batch_size //= increment
	for result in blast_output:
		p, f, b = fill_initial_info(result, flanking_n)
		
		# If the match couldn't be extended because it is at the end of
		# the contig then skip this one.
		if p == None:
			continue
		protos.append(p)
		pool_input.append((count, f, b))
		count += 1

	all_protospacer_infos = sequence_operations.pool_MP_blastdbcmd(
		pool_input,
		args.blast_db_path,
		args.blastdbcmd_batch_size,
		args.num_threads
		)


	p_count = 0
	for i in range(0,len(all_protospacer_infos),increment):
		p = protos[p_count]
		# Depending on which flanking sequence was requested, fill in
		# unrequested sequence with empty string.

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
			# Check for seed region mismatches if user specified.
			if seed_reg:
				if "." in p.aligned[seed_reg[0]:seed_reg[1]]:
					continue
			# Check for PAM if user requested
			if args.pam or args.regex_pam:
				if args.pam_location == 'up':
						seq = p.up
				else:
					seq = p.down
				if args.pam or args.regex_pam:
					if not bool(pam.fullmatch(seq)):
						no_pam_outcontents.append(
							"\t".join([str(i) for i in [
								p.spacer,
								p.target,
								p.start,
								p.stop,
								p.pid,
								p.mismatch,
								p.protoseq,
								p.aligned,
								p.up,
								p.down,
								p.strand]]).replace( #remove double tab if missing col
									"\t\t", "\t"))
						p_count += 1
						continue
			outcontents.append(
				"\t".join([str(i) for i in [
					p.spacer,
					p.target,
					p.start,
					p.stop,
					p.pid,
					p.mismatch,
					p.protoseq,
					p.aligned,
					p.up,
					p.down,
					p.strand]]).replace( #remove double tab if missing col
						"\t\t", "\t"))

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
