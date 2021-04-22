#!/usr/bin/env python3

import subprocess
import argparse
import sys
import os
from collections import Counter
from collections import OrderedDict

def fasta_to_dict(FASTA_file): # takes opened and read fasta file and returns a dict of format {fasta_header:sequence}
	fasta_dict = {}
	fastas = FASTA_file.split(">")
	trimmed_fastas = []
	for i in fastas:
		if len(i) != 0:
			trimmed_fastas.append(i)
	fastas = trimmed_fastas
	for i in fastas:
		header =  i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq
	return fasta_dict

def find_network_ends(edge_list): # takes list of edges from network file from narble output and finds the ends of each array and returns them as a list of tuples
	all_edges =[]
	for i in edge_list:
		all_edges += i[1:3]
	dodecamer_counts_dict = Counter(all_edges)
	ends =[]
	for k,v in dodecamer_counts_dict.items():
		if v ==1:
			ends.append(k)
		if v >2:
			print("Branch point identified with coverage above provided threshold for %s when running repeat %s. Manual curation of this output will be necessary." %(k, repname))
	paired = []
	distinct_array_ends = []
	for i in ends:
		for j in edge_list:
			if j[1] == i:
				if j[-1] == "LR":
					distinct_array_ends.append(i)
					break

	return distinct_array_ends


def find_array_order(twelvebp, edge_list): #ends[i], links)
	all_edges =[]
	for i in edge_list:
		all_edges += i[1:3]
	dodecamer_counts_dict = Counter(all_edges)

	ordered_dodecs = [twelvebp]
	chunk_of_interest = twelvebp
	edge_found = False
	while not edge_found:
		for edge in edge_list:
			if chunk_of_interest == edge[1]:
				chunk_of_interest = edge[2]
				edge_list.remove(edge)
				ordered_dodecs.append(chunk_of_interest)
				if dodecamer_counts_dict[chunk_of_interest] == 1:
					edge_found = True
			else:
				if edge == edge_list[-1]:
					print("Could not figure out the array starting with 12bp: %s" % twelvebp)
					edge_found = True
	return ordered_dodecs


def find_ordered_spacers(ordered_dodecs, spacer_file):
	with open(spacer_file, 'r') as sf:
		spacers_coverage_dict = Counter(list(fasta_to_dict(sf.read()).values()))
	ordered_spacer_dict = OrderedDict()
	count = 0
	for i in range(1, len(ordered_dodecs[:-2]), 2):
		count += 1
		nspacers = 0
		spacer_to_add = []
		for spacer in spacers_coverage_dict.keys():
			if spacer[:12] == ordered_dodecs[i] and spacer[-12:] == ordered_dodecs[i+1]:
				nspacers +=1
				spacer_to_add.append(spacer)
		if len(spacer_to_add) == 0:
			print("\n\nNo spacers found for the following 12bp chunks:")
			print(ordered_dodecs[i] + "\t" + ordered_dodecs[i+1])
		if nspacers > 1:
			print("\n\nMultiple spacers match the two 12bp chunks used by crhunt and narbl. Manual curation required to figure out what to do. Only spacer with the highest coverage is kept here. The 2 12bp chunks are: %s \t %s\nThe spacers identified are %s\nThe coverage of those spacers is (in order): %s\n\n" %(ordered_dodecs[i], ordered_dodecs[i+1], "\t".join(spacer_to_add), "\t".join([str(spacers_coverage_dict[x]) for x in spacer_to_add])))
			score = 0
			best = ''
			for conflict in spacer_to_add:
				if spacers_coverage_dict[conflict] > score:
					best = conflict
					score = spacers_coverage_dict[conflict]
			header = "_".join([args.id, repname, str(count), "cov", str(spacers_coverage_dict[best])])
			ordered_spacer_dict[header] = best
		else:
			header = "_".join([args.id, repname, str(count), "cov", str(spacers_coverage_dict[spacer_to_add[0]])])
			ordered_spacer_dict[header] = spacer_to_add[0]
	
	return ordered_spacer_dict



parser = argparse.ArgumentParser(
	description="runs narbl pipeline scripts written by Whitney England.\n\
	This pipeline was written to work with the following versions of scripts (if given):\
	crhunt4d.sh\
	narblWE.sh",
	formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
	"-s", dest="scripts_dir", required = True,
	help="Specify path to crhunt, narbl, and crisprmax scripts."
	)
parser.add_argument(
	"-o", dest="out_dir", required = False, default='', 
	help="Specify output directory path. Default is to write to current working directory"
	)
parser.add_argument(
	"-c", dest="cores", required = True, 
	help="Specify how many cpu cores to use."
	)
parser.add_argument(
	"-r", dest="repfile", required = True, 
	help="Specify file containing CRISPR repeats in fasta format."
	)
parser.add_argument(
	"-g", dest="readsfile", required = True, 
	help="Specify reads file in fasta format."
	)
parser.add_argument(
	"-n", dest="repsize", required = True, 
	help="Specify approximate repeat size."
	)
parser.add_argument(
	"-m", dest="mincov", required = True, type=int, 
	help="Specify minimum coverage of supporting 2 spacers being neighbours."
	)
parser.add_argument(
	"-i", dest="id", required = True, 
	help="Specify strain name or id for use in headers of output spacer fasta."
	)

args = parser.parse_args(sys.argv[1:])

if len(args.out_dir) !=0:
	main_outdir = os.path.abspath(args.out_dir)# + '/' if args.out_dir[-1] != '/' else args.out_dir
	if not os.path.isdir(main_outdir):
		os.mkdir(main_outdir)
else:
	main_outdir = os.path.abspath(os.getcwd())

scripts_dir = os.path.abspath(args.scripts_dir) # + '/' if args.scripts_dir[-1] != '/' else args.scripts_dir



repeats_dir = main_outdir + '/repeats/'
readsfile = os.path.abspath(args.readsfile)
repfile = os.path.abspath(args.repfile)


with open(repfile, 'r') as rin:
	repeats_dict = fasta_to_dict(rin.read())

repeats = list(repeats_dict.keys())

if len(repeats) > 1:
	
	if not os.path.isdir(repeats_dir):
		os.mkdir(repeats_dir)

	for k,v in repeats_dict.items():
		with open(repeats_dir + k + ".fna", 'w+') as repout:
			repout.write("\n".join([">"+k,v]))



for repname in repeats:
	outdir = main_outdir + '/' + repname + '/'

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	repfile = repeats_dir + repname + ".fna" if len(repeats) > 1 else repfile
	crhunt_out = outdir + 'crhunt_out/'
	narbl_out = outdir + 'narbl_out/'
	final_out = outdir + 'final_out/'

	if not os.path.isdir(crhunt_out):
		os.mkdir(crhunt_out)

	if not os.path.isdir(narbl_out):
		os.mkdir(narbl_out)

	if not os.path.isdir(final_out):
		os.mkdir(final_out)


	crhunt_command = "bash %s/crhunt4d.sh %s %s %s %s" %(scripts_dir, readsfile, repfile, crhunt_out[:-1], args.cores)

	print("\nrunning crhunt with the following command\n")
	print(crhunt_command + "\n")

	crhunt_p = subprocess.Popen(crhunt_command.split())
	crhunt_p.wait()

	if os.path.isfile(crhunt_out + repname + ".reads") and os.stat(crhunt_out + repname + ".reads").st_size == 0:
		print("\n\nNo repeats found for %s in %s\n\n" %(repname, args.id))
		continue

	# narbl outputs into the directory in which it is run. 
	# In order to organise the output we will run it in narbl_out and define all paths relative to that.

	narbl_command = "bash %s/narblWE.sh %s%s.reads %s %s" %(scripts_dir, crhunt_out, repname, repfile, args.repsize)

	print("\nrunning narbl with the following command\n")
	print(narbl_command + "\n")

	narbl_p = subprocess.Popen(narbl_command.split(), cwd=narbl_out)
	narbl_p.wait()

	links = []

	with open(narbl_out + repname + ".links.count.txt", 'r') as linksf:
		for line in linksf.readlines():
			if 'X' not in line and int(line.split()[0]) >= args.mincov:
				links.append(line.split())

	ends = find_network_ends(links)

	n_arrays_found = len(ends)

	for i in range(n_arrays_found):
		ordered_dodecs = find_array_order(ends[i], links)
		ordered_spacers_dict = find_ordered_spacers(ordered_dodecs, "%s%s.spacers.final.fasta" %(narbl_out, repname))
		with open(final_out + "_".join([args.id, repname, "array", str(i+1), "spacers.fna"]), 'w+') as spacer_outf:
			spacer_outf.write("\n".join([str(">" + k + "\n" + v) for k,v in ordered_spacers_dict.items()]) + '\n')
