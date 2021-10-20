#!/usr/bin/env python3

# AUTHOR		:	ALAN COLLINS
# VERSION		:	v1.3
# DATE			:	2021/2/16
# DESCRIPTION	:	Process minced output and produce csv and txt summaries of the information
# CHANGELOG
# v1.1
#	Added contig information to array locations.
# v1.2
#	Added option to provide repeats file
# v1.3
#	Added functionality to define path to minced executable


import argparse
import re
import os
import sys
import subprocess


class MincedObj():
	"""Class for reading and storing contents of minced output.

		Args:
			CRISPR_types_dict (dict): dictionary of fasta sequence of CRISPR repeats to compare to minced-identified repeats
			minced_file (str): Path to the minced output file to parse into this class instance

		Attributes:
			accession (str): name of file (up to extension) that minced was run on
			array_count (int): number of distinct arrays identified
			array_ids (dict): ID numbers of each array {Array_number : ID number}
			array_locations (dict): Locations of each array {Array_number : [contig, start, stop]}
			array_orientations (dict): dict of array numbers and bool of whether they are in the reverse direction {Array_number : reverse?}
			array_sizes (dict): Number of spacers in each array {Array_number : Spacer_count}
			arrays (dict): Sequences of spacers in each array {Array_number : [Spacer_sequences]}
			arrays_encoded (dict): ID numbers of spacers in each array {Array_number : [Spacer_IDs]}
			has_CRISPR (bool): Whether any CRISPR arrays were identified
			repeat (dict): For each array, the seqeunce of the mose common repeat {Array_number : repeat_sequence}
			repeat_scores (dict): For each array, the hamming distance when compared to the most similar repeat in the provided repeats file {Array_number : hammind_dist}
			repeat_types (dict): For each array, the most similar repeat in the provided repeats file {Array_number : repeat_type}
	"""
	def __init__(self, CRISPR_types_dict, minced_file):
		self.accession = minced_file.split("/")[-1].split("_minced")[0]
		if os.stat(minced_file).st_size == 0:
			self.has_CRISPR = False		
		else: 
			self.has_CRISPR = True
		self.array_count = 0
		self.arrays = {}
		self.arrays_encoded = {}
		self.array_ids = {}
		self.array_locations = {}
		self.array_sizes = {}
		self.repeat = {}
		self.repeat_types = {}
		self.repeat_scores = {}
		self.array_orientations = {}
		if self.has_CRISPR:
			with open(minced_file, 'r') as f:
				for line in f.readlines():
					if "Sequence" in line and "bp)" in line:
						contig = re.match(r"Sequence '(.*)' \([0-9]+ bp\)", line)[1]
					if "CRISPR" in line:
						self.array_count += 1
						array_num = int(re.match(r'CRISPR ([0-9]+) ', line)[1])
						self.array_locations[array_num] = [contig] + [i for i in re.findall(r'Range:\s+(\d+)\W+(\d+)\W+', line)[0]]
						self.arrays[array_num] = []
						self.arrays_encoded[array_num] = []
						self.array_sizes[array_num] = 0
						repeats = []
						spacers = []
					if "[" in line:
						self.array_sizes[array_num] += 1
						repeats.append(re.match(r'\d+\s+(\w+)\s+', line)[1])
						spacers.append(re.match(r'\d+\s+\w+\s+(\w+)\s+', line)[1])
					if "Repeats" in line:
						repeat = max(set(repeats), key=repeats.count)
						self.repeat_types[array_num], self.repeat_scores[array_num], self.array_orientations[array_num] = get_repeat_info(CRISPR_types_dict, repeat)
						if self.array_orientations[array_num]:
							self.repeat[array_num] = rev_comp(repeat)
							self.arrays[array_num] = [rev_comp(spacer) for spacer in reversed(spacers)]
						else:
							self.repeat[array_num] = repeat
							self.arrays[array_num] = spacers
					

			f.close()


def fasta_to_dict(FASTA_file):
	"""Read a fasta file into a dict with the headers (minus the > symbol) as keys and the associated sequence as values.
	
		Args:
			FASTA_file (str): path to fasta format file
	
		Returns:
			dict: dict of format {fasta_header : sequence}
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
		header = i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def rev_comp(string):
	"""Reverse complement a string of nucleotide sequence
	
		Args:
			string (str): Nucleotide sequence
	
		Returns:
			str: reverse complement of given nucleotide sequence
	"""
	rev_str = ''
	rev_comp_lookup = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c"}
	for i in reversed(string):
		if i in "ATCGatcg":
			rev_str += rev_comp_lookup[i]
		else:
			rev_str += i
	return rev_str


def hamming(string1, string2):
	"""Calculate hamming distance between two sequences. Compare them as provided, and with one sequence
	shifted one position to the left, and with the other sequence shifted one to the left.
	
		Args:
			string1 (str): First sequence to be compared
			string2 (str): Second sequence to be compared
	
		Returns:
			tuple:
				dist (int): the hamming distance of the unadjusted sequences
				distplus1 (int): the hamming distance of the first sequence shifted
				distminus1 (int): the hamming distance of the econd sequence shifted
	"""
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


def mince_it(minced_path, genomes_dir, out_dir):
	"""Run minced on the files in your input directory
	
		Args:
			minced_path (str): Path to the minced executable
			genomes_dir (str): path to the genomes you want to run minced on
			out_dir (str): path to the directory you want to save minced output files. A dir called MINCED/ will be created and files stored in there.
	"""
	minced_out = out_dir + "MINCED_OUT/"
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	if not os.path.isdir(minced_out):
		os.mkdir(minced_out)

	for file in os.listdir(genomes_dir):
		genomefile = str(genomes_dir+file)
		accession = file.split(".f")[0]

		command = minced_path + " " + genomefile + " " + minced_out + accession + "_minced_out.txt"
		p = subprocess.Popen(command, shell=True)
		p_status = p.wait()



def get_repeat_info(CRISPR_types_dict, repeat):
	"""Figure out which of the repeats the user provides in the repeats fasta file is the most similar to the repeat found by minced
	
		Args:
			CRISPR_types_dict (dict): fasta dict of the repeats fasta file as returned by fasta_to_dict
			repeat (str): The repeat identified by minced for which you want to find a match
	
		Returns:
			tuple: 
				best_match (str): the fasta header of the user-provided repeat that is most similar to that found by minced
				best_score (int): the hamming distance between the best match repeat and the minced-identified repeat
				reverse (bool): Is the best match repeat the reverse complement of the identified repeat?
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

		score = min(hamming(rev_comp(repeat), v))
		if score < best_score:
			best_score = score
			best_match = k
			reverse = True


	return best_match, best_score, reverse


def process_minced_out(CRISPR_types_dict, out_dir):
	"""Read in minced output files and aggregate the data into a more easily analysed form.
	
		Args:
			CRISPR_types_dict (dict): fasta dict of the repeats fasta file as returned by fasta_to_dict
			out_dir (str): Path to the output directory into which minced files were saved in the MINCED/ dir.
						   A dir called PROCESSED will be created and the files generated from processed minced output
						   will be stored there.
	"""
	Minced_out_dir = out_dir + "MINCED_OUT/" 
	ProcessedOut = out_dir + "PROCESSED/"

	if not os.path.isdir(ProcessedOut):
		os.mkdir(ProcessedOut)

	allfiles = list()

	for file in os.listdir(Minced_out_dir):
		infile = str(Minced_out_dir+file)
		allfiles.append(MincedObj(CRISPR_types_dict, infile))
	all_spacers = []
	all_arrays = []
	for strain in allfiles:
		for k,v in strain.arrays.items():
			all_arrays.append(" ".join(v))
			for s in v:
				all_spacers.append(s)

	non_red_spacers = set(all_spacers)
	non_red_spacer_dict = {k:v for k,v in zip(non_red_spacers, range(1,len(non_red_spacers)+1))}

	non_red_arrays = set(all_arrays)
	non_red_array_dict = {k:v for k,v in zip(non_red_arrays, range(1,len(non_red_arrays)+1))}

	for i in range(len(allfiles)):
		for k,v in allfiles[i].arrays.items():
			if " ".join(v) in non_red_array_dict.keys():
				allfiles[i].array_ids[k] = non_red_array_dict[" ".join(v)]
			else:
				print("Process_minced_out says: Couldn't find array %s\nfrom %s" %(v, allfiles[i].accession))
			for s in v:
				if s in non_red_spacer_dict.keys():
					allfiles[i].arrays_encoded[k].append(str(non_red_spacer_dict[s]))
				else:
					print("Process_minced_out says: Couldn't find spacer %s\nfrom %s" %(s, allfiles[i].accession))


	print("Total unique spacers: %i" %len(non_red_spacers))
	print("Total unique arrays: %i" %len(non_red_arrays))

	with open(ProcessedOut + "CRISPR_FASTAs.txt", 'w+') as f:
		for k,v in non_red_spacer_dict.items():
			f.write(">%s\n%s\n" %(v,k))
	f.close()

	with open(ProcessedOut + "Array_IDs.txt", 'w+') as f:
		for k,v in non_red_array_dict.items():
			f.write("%s\tspacer_IDs:\t%s\n" %(v," ".join([str(non_red_spacer_dict[i]) for i in k.split()])))
	f.close()

	with open(ProcessedOut + "Array_seqs.txt", 'w+') as f:
		for k,v in non_red_array_dict.items():
			f.write("%s\tspacer_seq:\t%s\n" %(v,k))
	f.close


	with open(ProcessedOut + "CRISPR_summary_table.csv", 'w') as f:
		headers = "Strain,Has_CRISPR,Array_count,Spacers,Spacers_encoded,Array_IDs,Array_locations,Array_sizes,Repeat_sequences,Array_CRISPR_types,Repeat_scores"
		f.write(headers)
		for i in allfiles:
			Strain = i.accession
			Has_CRISPR = str(i.has_CRISPR)
			Array_count = i.array_count
			Spacers = []
			Spacers_encoded = []
			Array_ids = []
			Array_locations = []
			Array_sizes = []
			Repeat_sequences = []
			Array_types = []
			Repeat_scores = []
			for k,v in i.arrays.items():
				Spacers.append(str(k) + ": " + " ".join(v))
			for k,v in i.arrays_encoded.items():
				Spacers_encoded.append(str(k) + ": " + " ".join(v))
			for k,v in i.array_ids.items():
				Array_ids.append("'" + str(k) + ": " + str(v) + "'")
			for k,v in i.array_locations.items():
				Array_locations.append(str(k) + ": " + " ".join(v))
			for k,v in i.array_sizes.items():
				Array_sizes.append("'" + str(k) + ": " + str(v) +"'")
			for k,v in i.repeat.items():
				Repeat_sequences.append(str(k) + ": " + v)
			for k,v in i.repeat_types.items():
				Array_types.append(str(k) + ": " + v)
			for k,v in i.repeat_scores.items():
				Repeat_scores.append("'" + str(k) + ": " + str(v) + "'")
			Spacers = "\t".join(Spacers)
			Spacers_encoded = "\t".join(Spacers_encoded)
			Array_ids = "\t".join(Array_ids)
			Array_locations = "\t".join(Array_locations)
			Array_sizes = "\t".join(Array_sizes)
			Repeat_sequences = "\t".join(Repeat_sequences)
			Array_types = "\t".join(Array_types)
			Repeat_scores = "\t".join(Repeat_scores)
			
			f.write("\n%s,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s" %(Strain,Has_CRISPR,Array_count,Spacers,Spacers_encoded,Array_ids,\
				Array_locations,Array_sizes,Repeat_sequences,Array_types,Repeat_scores))

	with open(ProcessedOut + "CRISPR_summary_table.txt", 'w') as f:
			headers = "Strain\tHas_CRISPR\tArray_count\tSpacers\tSpacers_encoded\tArray_IDs\tArray_locations\tArray_sizes\tRepeat_sequences\tArray_CRISPR_types\tRepeat_scores"
			f.write(headers)
			for i in allfiles:
				Strain = i.accession
				Has_CRISPR = str(i.has_CRISPR)
				Array_count = i.array_count
				Spacers = []
				Spacers_encoded = []
				Array_ids = []
				Array_locations = []
				Array_sizes = []
				Repeat_sequences = []
				Array_types = []
				Repeat_scores = []
				for k,v in i.arrays.items():
					Spacers.append(" ".join(v))
				for k,v in i.arrays_encoded.items():
					Spacers_encoded.append(" ".join(v))
				for k,v in i.array_ids.items():
					Array_ids.append(str(v))
				for k,v in i.array_locations.items():
					Array_locations.append(" ".join(v))
				for k,v in i.array_sizes.items():
					Array_sizes.append(str(v))
				for k,v in i.repeat.items():
					Repeat_sequences.append(v)
				for k,v in i.repeat_types.items():
					Array_types.append(v)
				for k,v in i.repeat_scores.items():
					Repeat_scores.append(str(v))
				Spacers = "|".join(Spacers)
				Spacers_encoded = "|".join(Spacers_encoded)
				Array_ids = "|".join(Array_ids)
				Array_locations = "|".join(Array_locations)
				Array_sizes = "|".join(Array_sizes)
				Repeat_sequences = "|".join(Repeat_sequences)
				Array_types = "|".join(Array_types)
				Repeat_scores = "|".join(Repeat_scores)
				
				f.write("\n%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(Strain,Has_CRISPR,Array_count,Spacers,Spacers_encoded,Array_ids,\
					Array_locations,Array_sizes,Repeat_sequences,Array_types,Repeat_scores))


		




if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description="Runs minced and processes the output into a spacer multifasta and summary tables describing which genomes have which CRISPRs")
	parser.add_argument(
		"-i",  dest="indir", required = False,
		help="Specify the input directory containing genome fastas. "
		)
	parser.add_argument(
		"-l",  dest="minced_path", required = False,
		help="If running minced, specify the path to the minced executable. "
		)
	parser.add_argument(
		"-o", dest="out_dir", required = True,
		help="Specify directory for minced output and processed minced files. If just running minced processing steps then this is the directory above the directory where your minced output is stored. e.g. for dir1/MINCED_OUT/minced_output_files, use -o dir1. Files generated by processing the minced outputs will be stored in dir1/PROCESSED/"
		)
	parser.add_argument(
		"-r", dest="repeats_file", required = False,
		help="FASTA format file containing the CRISPR repeats you want to look for. If you don't provide one, this script will look for type 1F, 1E, and 1C repeats from Pseudomonas aeruginosa. These are used to classify arrays and orient them consistently with one another."
		)
	parser.add_argument(
		"-m", action='store_true',
			help="Indicate that you want minced to be run"
		)
	parser.add_argument(
		"-p", action='store_true',
			help="Indicate that you want to run processing steps on minced output"
		)

	args = parser.parse_args(sys.argv[1:])

	indir = args.indir
	out_dir = args.out_dir

	if indir and indir[-1] != "/":
		indir += "/"
	if out_dir[-1] != "/":
		out_dir+= "/"

	if args.repeats_file:
		CRISPR_types_dict = fasta_to_dict(args.repeats_file)
	else:
		CRISPR_types_dict = { #Specify the desired orientation and type of all repeats you are interested in here or provide them as a fasta file.
		'1E':'GTGTTCCCCACGGGTGTGGGGATGAACCG',
		'1F':'GTTCACTGCCGTGTAGGCAGCTAAGAAA',
		'1C':'GTCGCGCCCCGCACGGGCGCGTGGATTGAAAC'
		}

	if args.m:
		if not args.minced_path:
			print("ERROR: To run minced you must provide the path to the minced executable.")
			sys.exit()
		mince_it(args.minced_path, indir, out_dir)
	if args.p:
		process_minced_out(CRISPR_types_dict, out_dir)
