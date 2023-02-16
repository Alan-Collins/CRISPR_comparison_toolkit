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
	return parser


def run_minced(args):
	parser = argparse.ArgumentParser(
		description="temp parser")
	minced_parser = minced.build_parser(parser)
	minced_args = minced_parser.parse_args([
		"-i", args.indir,
		"-o", args.outdir,
		"--min-shared", "2",
		"-mp"
		])
	sys.stderr.write("\nrunning cctk minced on the provided assemblies...\n")
	minced.main(minced_args)

def load_network(network_file):
	# dict where the key is an array and the value is its parent array.
	# Whenever an array gets absorbed into an existing network, the
	# array it shared spacers with becomes the parent array. Following
	# the chain of parents should terminate at the array that is the
	# namesake of the cluster in network_dict.
	parent_dict = {}
	# dict where the key is a namesake array for the cluster and value
	# is a list of arrays in that cluster.
	network_dict = {}
	for line in network_file.readlines()[1:]:
		array1, array2 = line.split()[:2]	
		if array1 in parent_dict.keys():
			# If array is already in a cluster, follow the chain of
			# parents to the namesake at the top.
			parent1 = parent_dict[array1]
			while parent1 in parent_dict.keys():
				# Follow the parents up the chain to find the top of the pile.
				parent1 = parent_dict[parent1]
		else:
			parent1 = array1
		if array2 in parent_dict.keys():
			parent2 = parent_dict[array2]
			while parent2 in parent_dict.keys():
				parent2 = parent_dict[parent2]
		else:
			parent2 = array2
		if parent1 == parent2:
			# These are already clustered together so we can move on
			continue
		if parent1 in network_dict.keys() and parent2 in network_dict.keys():
			 # Add the list of cluster 2 members to cluster 1
			network_dict[parent1] = list(
				set(network_dict[parent1] + network_dict[parent2])
				)
			# Remove the now absorbed cluster 2 from the dict
			del network_dict[parent2]
			parent_dict[array2] = parent1
		elif parent1 in network_dict.keys():
			# If parent2 not in network_dict.keys() then parent2 must
			# be the same as array2 so no need to worry about adding
			# a whole chain to the existing array
			network_dict[parent1].append(array2)
			parent_dict[array2] = parent1
		elif parent2 in network_dict.keys():
			network_dict[parent2].append(array1)
			parent_dict[array1] = parent2
		else:
			# If neither are in the network_dict, then create a new
			# cluster with just the 2 arrays in it
			network_dict[parent1] = [parent1, parent2]
			parent_dict[array2] = parent1
	for k,v in network_dict.items():
		network_dict[k] = list(set(v))

	return network_dict


def run_crisprdiff(network, array_ids, plotdir):
	sys.stderr.write("\nrunning cctk crisprdiff on the identified array clusters...\n")
	for cluster_id, cluster_members in network:
		sys.stderr.write(f"\nprocessing cluster {cluster_id}...\n")
		parser = argparse.ArgumentParser(
			description="temp parser")
		crisprdiff_parser = CRISPRdiff.build_parser(parser)
		crisprdiff_args = crisprdiff_parser.parse_args([
			"-a", array_ids,
			"-o", f"{plotdir}{cluster_id}_diffplot.png",
			] + cluster_members)
		
		CRISPRdiff.main(crisprdiff_args)

def run_crisprtree(network, array_ids, plotdir):
	sys.stderr.write("\nrunning cctk crisprtree on the identified array clusters...\n")
	for cluster_id, cluster_members in network:
		sys.stderr.write(f"\nprocessing cluster {cluster_id}...\n")
		parser = argparse.ArgumentParser(
			description="temp parser")
		crisprtree_parser = CRISPRtree.build_parser(parser)
		crisprtree_args = crisprtree_parser.parse_args([
			"-a", array_ids,
			"-o", f"{plotdir}{cluster_id}_tree.png",
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
	run_minced(args)

	# Identify clusters
	with open(outdir + "PROCESSED/Array_network.txt") as network_file:
		network = load_network(network_file)

	# filter clusters so they run in a sensible amount of time
	filt_network = [(k,v) for k,v in network.items() if len(v) >3 and len(v) < 15]

	# plot diffplots
	run_crisprdiff(filt_network, outdir + "PROCESSED/Array_IDs.txt", plotdir)

	# plot crisprtrees
	run_crisprtree(filt_network, outdir + "PROCESSED/Array_IDs.txt", plotdir)
