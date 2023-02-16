#!/usr/bin/env python3

import sys
import os
import argparse

from . import minced, CRISPRdiff, CRISPRtree

description = """
usage: cctk quickrun -i <path> -o <path> [-h]

optional arguments:
  -h, --help  show this help message and exit

Required arguments:
  -i, --indir         	input directory containing genome fastas.
  -o, --outdir        	directory to write output files
  -m, --max-cluster     Largest cluster size to plot. Default: 15
"""

def build_parser(parser):
	
	req_options = parser.add_argument_group("Required arguments")
	req_options.add_argument(
		"-o", "--outdir",
		metavar=" ",
		required = True,
		help="directory to write output files"
		)

	req_options.add_argument(
		"-i", "--indir",
		metavar=" ",
		required=False,
		help="Specify the input directory containing genome fastas."
		)

	other_options = parser.add_argument_group("Other inputs")
	other_options.add_argument(
		"-m", "--max_cluster",
		metavar=" ",
		type=int,
		default=15,
		required = False,
		help="Largest cluster size to plot. Default: 15"
		)

	return parser


def run_minced(indir, outdir):
	parser = argparse.ArgumentParser(
		description="temp parser")
	minced_parser = minced.build_parser(parser)
	minced_args = minced_parser.parse_args([
		"-i", indir,
		"-o", outdir,
		"--min-shared", "2",
		"-mp"
		])
	sys.stderr.write("\nrunning cctk minced on the provided assemblies...\n")
	minced.main(minced_args)

def run_crisprdiff(clusters, array_ids, plotdir):
	sys.stderr.write("\nrunning cctk crisprdiff on the identified array clusters...\n")
	for n, cluster_members in enumerate(clusters):
		sys.stderr.write(f"\nprocessing cluster {n+1}...\n")
		parser = argparse.ArgumentParser(
			description="temp parser")
		crisprdiff_parser = CRISPRdiff.build_parser(parser)
		crisprdiff_args = crisprdiff_parser.parse_args([
			"-a", array_ids,
			"-o", f"{plotdir}cluster_{n+1}_diffplot.png",
			] + cluster_members)
		
		CRISPRdiff.main(crisprdiff_args)

def run_crisprtree(clusters, array_ids, plotdir):
	sys.stderr.write("\nrunning cctk crisprtree on the identified array clusters...\n")
	for n, cluster_members in enumerate(clusters):
		sys.stderr.write(f"\nprocessing cluster {n+1}...\n")
		parser = argparse.ArgumentParser(
			description="temp parser")
		crisprtree_parser = CRISPRtree.build_parser(parser)
		crisprtree_args = crisprtree_parser.parse_args([
			"-a", array_ids,
			"-o", f"{plotdir}cluster_{n+1}_tree.png",
			] + cluster_members)
		
		CRISPRtree.main(crisprtree_args)

def main(args):

	indir = args.indir
	outdir = args.outdir
	
	if indir and indir[-1] != "/":
		indir += "/"
	if outdir[-1] != "/":
		outdir+= "/"

	plotdir = outdir + "PLOTS/"

	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	if not os.path.isdir(plotdir):
		os.makedirs(plotdir)

	# identify CRISPRs with minced
	run_minced(indir, outdir)

	# Load clusters
	clusters = []
	with open(outdir + "PROCESSED/Array_clusters.txt") as cluster_file:
		for line in cluster_file:
			cluster = line.split()
			clusters.append(line.split())
	
	# filter clusters so they run in a sensible amount of time
	filt_clusters = [
		c for c in clusters if len(c) >=3 and len(c) <= args.max_cluster]

	sys.stderr.write(f"Found {len(filt_clusters)} clusters with between 3 and {args.max_cluster} arrays.\n")
	if len(filt_clusters) == 0:
		sys.stderr.write("Exiting.\n")
		sys.exit()

	# plot diffplots
	run_crisprdiff(filt_clusters, outdir + "PROCESSED/Array_IDs.txt", plotdir)

	# plot crisprtrees
	run_crisprtree(filt_clusters, outdir + "PROCESSED/Array_IDs.txt", plotdir)
