#! /usr/bin/env python3

import argparse
import subprocess
import sys
import multiprocessing
import os
import textwrap as _textwrap



def run_CRISPR_cas_finder(infile):
	newdir = infile.split('/')[-1].split('.')[0]
	
	print("Running {}".format(newdir))

	if not os.path.isdir(newdir):
		os.mkdir(newdir)
	else:
		print("It appears you have already run {}. Skipping...".format(newdir))
		return

	with open(newdir + '/' + newdir + "_stdout.log", 'w') as logfile1:
		with open(newdir + '/' + newdir + "_stderr.log", 'w') as logfile2:

			command = "cd {}; perl /home/alanc/Downloads/CRISPRCasFinder/CRISPRCasFinder-release-4.2.20/CRISPRCasFinder.pl -in {} -so /home/alanc/Downloads/CRISPRCasFinder/CRISPRCasFinder-release-4.2.20/sel392v2.so -cas".format(newdir, infile)

			p = subprocess.run(command, shell=True, stdout=logfile1, stderr=logfile2)#check=True
		logfile2.close()
	logfile1.close()



def pool_MP_CRISPR_cas_finder(filelist, cores, chunksize):

	pool = multiprocessing.Pool(processes=int(cores))

	pool.imap_unordered(run_CRISPR_cas_finder, filelist, chunksize)
	pool.close()
	pool.join()

def build_parser(parser):
	parser.add_argument(
		"-i",  dest="indir", required = False,
		help="Specify input directory containing fasta format genomes. Cannot be combined with input file(s)."
		)
	parser.add_argument(
		"-f",  dest="infile", required = False, nargs="+",
		help="Specify input file(s) containing fasta format genomes. Cannot be combined with input directory."
		)
	parser.add_argument(
		"-c",  dest="cores", type=int, nargs="?", default = 1,
			help="Specify number of cores to use. Default: 1"
		)
	return parser


def main(args):

	if args.indir:
		indir = args.indir + '/' if args.indir[-1] != '/' else args.indir
		filelist = [indir + i for i in os.listdir(indir)]
	elif args.infile:
		filelist = args.infile
	else:
		print("You must provide either an input directory or a list of one or more files on which you want this program to run.")
		quit()
		

	chunksize = len(filelist)//args.cores

	print("Running with {} cores on {} files. Estimated time to completion: {} minutes (assuming approx 1 minute per genome)".format(args.cores, len(filelist), len(filelist)/args.cores))

	pool_MP_CRISPR_cas_finder(filelist, args.cores, chunksize)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(
	description="Given a directory containing fasta genomes, runs \
	CRISPR_cas_finder on all the genomes and stores the outputs in your \
	current working directory.",
	)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)
