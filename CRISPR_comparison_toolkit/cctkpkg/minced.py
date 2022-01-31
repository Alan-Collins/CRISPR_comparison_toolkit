#!/usr/bin/env python3

import argparse
import re
import os
import sys
import subprocess
from shutil import which
from collections import Counter, defaultdict

from . import file_handling, sequence_operations

def mince_it(minced_path, genomes_dir, outdir):
	"""Run minced on the files in your input directory
	
	Args:
	  minced_path (str):
		Path to the minced executable
	  genomes_dir (str):
		Path to the genomes you want to run minced on
	  outdir (str):
		Path to the directory you want to save minced output files. 
		A dir called MINCED/ will be created and files stored in there.
	"""
	minced_out = outdir + "MINCED_OUT/"
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	if not os.path.isdir(minced_out):
		os.mkdir(minced_out)

	for file in os.listdir(genomes_dir):
		genomefile = str(genomes_dir+file)
		accession = file.split(".f")[0]

		command = " ".join([
			minced_path,
			genomefile,
			minced_out
			]) + accession + "_minced_out.txt"
		p = subprocess.Popen(command, shell=True)
		p_status = p.wait()


def process_minced_out(CRISPR_types_dict, outdir):
	""" Process minced output into human-friendly/useful formats

	Read in minced output files and aggregate the data into a more 
	easily analysed form.
	
	Args:
	  CRISPR_types_dict (dict):
		fasta dict of the repeats fasta file as returned by fasta_to_dict
	  outdir (str):
		Path to the output directory into which minced files were saved 
		in the MINCED/ dir.
	
	A dir called PROCESSED will be created and the files generated from
	processed minced output will be stored there.
	"""
	Minced_outdir = outdir + "MINCED_OUT/" 
	processed_out = outdir + "PROCESSED/"

	if not os.path.isdir(processed_out):
		os.mkdir(processed_out)

	all_assemblies = []

	for file in os.listdir(Minced_outdir):
		infile = Minced_outdir+file
		all_assemblies.append(file_handling.AssemblyCRISPRs(
			CRISPR_types_dict, infile))

	(non_red_spacer_dict,
		non_red_spacer_id_dict,
		non_red_array_dict,
		non_red_array_id_dict
	) = sequence_operations.non_redundant_CR(all_assemblies)

	# Fill in spacer and array ID info
	sequence_operations.add_ids(
		all_assemblies,
		non_red_spacer_id_dict,
		non_red_array_id_dict)

	# Output summary info
	sys.stderr.write("Total unique spacers: %i\n" %len(non_red_spacer_id_dict))
	sys.stderr.write("Total unique arrays: %i\n" %len(non_red_array_id_dict))

	# Save array files
	file_handling.write_CRISPR_files(all_assemblies,
	non_red_spacer_id_dict,
	non_red_array_id_dict,
	processed_out)

	# Make list of FoundArray instances representing all unique arrays
	array_dict = {}
	for assembly in all_assemblies:
		for array in assembly.arrays.values():
			if array.id not in array_dict:
				array_dict[array.id] = array
	array_list = [a for a in array_dict.values()]

	# Build network of array spacer sharing
	network = sequence_operations.build_network(array_list)

	# Write network file
	file_handling.write_network_file(network, processed_out+"Array_network.txt")


def build_parser(parser):
	
	req_options = parser.add_argument_group("Required arguments")
	req_options.add_argument(
		"-o", "--outdir",
		metavar=" ",
		required = True,
		help="directory for minced output and processed minced files"
		)

	other_options = parser.add_argument_group("Other inputs")
	other_options.add_argument(
		"-i", "--indir",
		metavar=" ",
		required=False,
		help="Specify the input directory containing genome fastas."
		)
	other_options.add_argument(
		"-l", "--minced-path",
		metavar=" ",
		required=False,
		default="minced",
		help="If running minced, specify the path to the minced \
		executable. Default is to run minced as if it were in your PATH"
		)
	other_options.add_argument(
		"-r", "--repeats",
		metavar=" ",
		dest="repeats_file",
		required = False,
		help="FASTA format file containing the CRISPR repeats you want \
		to look for. If you don't provide one, this script will look \
		for type 1F, 1E, and 1C repeats from Pseudomonas aeruginosa. \
		These are used to classify arrays and orient them consistently \
		with one another."
		)

	run_options = parser.add_argument_group("Specify run type")
	run_options.add_argument(
		"-m", "--run-minced",
		action='store_true',
		help="Indicate that you want minced to be run"
		)
	run_options.add_argument(
		"-p", "--process-minced",
		action='store_true',
		help="Indicate that you want to run processing steps on minced \
		output"
		)

	return parser


def main(args):
	indir = args.indir
	outdir = args.outdir

	if indir and indir[-1] != "/":
		indir += "/"
	if outdir[-1] != "/":
		outdir+= "/"

	if args.repeats_file:
		crispr_types_dict = file_handling.fasta_to_dict(args.repeats_file)
	else:
		# Default repeats to make sure identical repeats are
		# oriented consistently when processing minced output.
		crispr_types_dict = {
		'1E':'GTGTTCCCCACGGGTGTGGGGATGAACCG',
		'1F':'GTTCACTGCCGTGTAGGCAGCTAAGAAA',
		'1C':'GTCGCGCCCCGCACGGGCGCGTGGATTGAAAC'
		}

	if args.run_minced:
		if which(args.minced_path) is None:
			sys.exit("ERROR: minced executable not found. \
				To run minced you must provide the path to the \
				minced executable.")
		mince_it(args.minced_path, indir, outdir)
	if args.process_minced:
		process_minced_out(crispr_types_dict, outdir)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description="Use minced to identify CRISPR arrays.")
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)
	
