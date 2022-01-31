import os
import sys
import json
import re
from collections import Counter, defaultdict

# cctkpkg imports
from . import sequence_operations 

class AssemblyCRISPRs():
	"""Class for reading and storing contents of minced and blast output.

	Args:
	  CRISPR_types_dict (dict):
		dictionary of fasta sequence of CRISPR repeats to compare to 
		minced-identified repeats
	  minced_file (str):
	  Path to the minced output file to parse into this class instance

	Attributes:
	  accession (str):
		name of file (up to extension) that minced was run on
	  array_count (int):
		number of distinct arrays identified
	  array_ids (dict):
		ID numbers of each array {Array_number : ID number}
	  array_locations (dict):
		Locations of each array {Array_number : [contig, start, stop]}
	  array_orientations (dict):
		dict of array numbers and bool of whether they are in the 
		reverse direction {Array_number : reverse?}
	  array_sizes (dict):
		Number of spacers in each array {Array_number : Spacer_count}
	  arrays (dict):
		Sequences of spacers in each array 
		{Array_number : [Spacer_sequences]}
	  arrays_encoded (dict):
		ID numbers of spacers in each array 
		{Array_number : [Spacer_IDs]}
	  has_CRISPR (bool):
		Whether any CRISPR arrays were identified
	  repeat (dict):
		For each array, the seqeunce of the mose common repeat 
		{Array_number : repeat_sequence}
	  repeat_scores (dict):
		For each array, the hamming distance when compared to the most 
		similar repeat in the provided repeats file 
		{Array_number : hammind_dist}
	  repeat_types (dict):
		For each array, the most similar repeat in the provided 
		repeats file {Array_number : repeat_type}
	"""
	def __init__(self, crispr_types_dict=None, minced_file=None):
		self.accession = ''
		self.array_count = 0
		self.arrays = {}
		self.has_crispr = False
		if minced_file:
			self.read_minced(crispr_types_dict, minced_file)
			
					

	def read_minced(self, crispr_types_dict, minced_file):
		self.accession = minced_file.split("/")[-1].split("_minced")[0]
		if os.stat(minced_file).st_size == 0:
			self.has_crispr = False
			return
		else: 
			self.has_crispr = True
		with open(minced_file, 'r') as f:
			for line in f.readlines():
				if "Sequence" in line and "bp)" in line:
					contig = re.match(
						r"Sequence '(.*)' \([0-9]+ bp\)",
						line)[1]
					continue
				if "CRISPR" in line and "Range:" in line:
					self.array_count+=1
					array_num = int(re.match(
							r'CRISPR ([0-9]+) ',
							line)[1])
					array = FoundArray()
					array.genome = self.accession
					array.contig = contig
					array.start, array.stop = [
						i for i in re.findall(
							r'Range:\s+(\d+)\W+(\d+)\W+', 
							line)[0]]
					continue
				if "[" in line:
					array.repeats.append(re.match(
						r'\d+\s+(\w+)\s+', 
						line)[1])
					array.spacers.append(re.match(
						r'\d+\s+\w+\s+(\w+)\s+', 
						line)[1])
					continue
				if "Repeats" in line:
					repeat = Counter(array.repeats).most_common(1)[0][0]
					(array.repeat_id, 
						array.repeat_score, 
						array.reverse
					) = sequence_operations.get_repeat_info(
						crispr_types_dict, 
						repeat)

					if array.reverse:
						array.repeats = [
						sequence_operations.rev_comp(
							repeat) for repeat in reversed(array.repeats)]
						array.spacers = [
						sequence_operations.rev_comp(
							spacer) for spacer in reversed(array.spacers)]	

					# If there are lots of mismatches between the
					# repeat and its best match. It's probably wrong.
					if array.repeat_score > 5:
						array.repeat_id = "unknown_CRISPR_type({})".format(
							array.repeat_id)
					self.arrays[array_num] = array
						

class FoundArray():
	"""
	Class to store information about a CRISPR array

	Where the array was found and which spacers it contains.
	
	Attributes:
		genome (str) ID of the genome in which this array was identified
		contig (str) fasta header of the sequence in which this array was identified (blast result sseqid)
		start (int) start position of array in contig
		stop (int) end position of array in contig
		repeats (list of strs) sequence of repeats
		repeat_id (str) with which repeat was this array identified?
		repeat_score (int) number of mismatches between repeat and best match
		reverse (bool) Was the array found in the reverse direction?
		spacers (list of strs) list of the seqeuences of the identified spacers
		spacer_id (list of ints) list of the IDs of the identified spacers as stored in spacer_dict
		id (int) ID of this array as stored in array_dict
	"""
	def __init__(self):
		self.genome = ''
		self.contig = ''
		self.start = 0
		self.stop = 0
		self.repeats = []
		self.repeat_id = ''
		self.repeat_score = 0
		self.reverse = False 
		self.spacers = []
		self.spacer_ids = []
		self.id = 0


def read_array_file(file):
	"""Read array spacers file into dict.
	
	Read Array IDs or sequences from file produced by one of the CRISPR
	finding scripts. Expected format of file is:
	Array_ID\tCRISPR_type\tspacer1 spacer2 etc

	Args:
	  file (str): 
		path to Array file

	Returns:
	  dict: 
		dict of format {Array_ID : [list, of, spacers]}

	Raises:
	  TypeError: If file is not a str
	  OSError: If file is not the path to an existing file

	"""

	if type(file) is not str:
		raise TypeError(
			"file must be str, not {}.".format(type(file).__name__))

	if not os.path.exists(file):
		raise OSError(
			"file must be the path to an existing file.")

	array_spacers_dict = {}
	with open(file, 'r') as fin:
		for line in fin.readlines():
			bits = line.split()
			array_spacers_dict[bits[0]] = bits[1:]
	return array_spacers_dict


def fasta_to_dict(FASTA_file):
	"""Read a fasta file into a dict 

	Dict has headers (minus the > symbol) as keys and the associated 
	sequence as values.
	
	Args:
	  FASTA_file (str): 
		path to fasta format file

	Returns:
	  dict: 
		dict of format {fasta_header : sequence}

	Raises:
	  TypeError: If FASTA_file is not a str
	  OSError: If FASTA_file is not the path to an existing file
	"""
	
	if type(FASTA_file) is not str:
		raise TypeError(
			"FASTA_file must be str, not {}.".format(type(FASTA_file).__name__))

	if not os.path.exists(FASTA_file):
		raise OSError(
			"FASTA_file must be the path to an existing file.")


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


def read_colour_scheme(infile):
	with open(infile, 'r') as fin:
		spacer_colours = json.load(fin)
	return spacer_colours


def write_colour_scheme(outfile, colours):
	with open(outfile, 'w', encoding='utf-8') as fout:
		json.dump(colours,
			fout,
			ensure_ascii=False,
			indent=4)


def read_colour_file(colour_file):
	with open(colour_file, 'r') as fin:
		cf_list = [
		i.strip().replace('"', '').replace("'", "") for i in fin.readlines()]
	return cf_list


def write_array_file(array_dict, outfile):
	with open(outfile, 'w') as fout:
		fout.write('\n'.join(
			["{}\t{}".format(
				k, " ".join(v.spacers)
				) for k,v in array_dict.items()])+"\n")


def write_spacer_fasta(spacer_dict, outfile):
	with open(outfile, 'w') as fout:
		for seq, sp_id in spacer_dict.items():
			fout.write(f">{sp_id}\n{seq}\n")


def write_array_id_seq_file(spacer_id_dict, array_id_dict, outdir):

	# seq file
	with open(outdir + "Array_seqs.txt", 'w') as fout:
		for sps, array_id in array_id_dict.items():
			fout.write(f"{array_id}\t{sps}\n")

	# ID file
	with open(outdir + "Array_IDs.txt", 'w') as fout:
		for sps, array_id in array_id_dict.items():
			fout.write("{}\t{}\n".format(
				array_id,
				" ".join([str(spacer_id_dict[sp]) for sp in sps.split()])))


def write_cr_sum_tabs(all_assemblies, outfile):
	# Determine delimiters based on file extension
	if outfile[-3:] == "csv":
		maj_delim = ","
		min_delim = "\t"
		inc_count = True
		# additional nonsense to stop excel assuming dates
		prepend_stuff = '"=""'
		append_stuff = '"""'
	else:
		maj_delim = "\t"
		min_delim = "|"
		inc_count = False

	outcontents = maj_delim.join([
			"Sequence_ID",
			"Has_CRISPR",
			"Array_count",
			"Spacers",
			"Spacer_IDs",
			"Array_IDs",
			"Array_locations",
			"Repeat_sequences",
			"Array_CRISPR_types",
			]) + "\n"
	for entry in all_assemblies:
		# Construct line components
		if inc_count: #If this is for excel, add the array number to the info
			sequence_id = prepend_stuff + entry.accession + append_stuff
			has_crispr = str(entry.has_crispr)
			array_count = str(entry.array_count)
			spacers = min_delim.join(
				"{}: {}".format(
					array_num,
					" ".join(array.spacers)
					) for array_num, array in entry.arrays.items()
				)
			spacer_ids = prepend_stuff + min_delim.join(
				"{}: {}".format(
					array_num,
					" ".join(array.spacer_ids)
					) for array_num, array in entry.arrays.items()
				) + append_stuff
			array_ids = prepend_stuff + min_delim.join(
				"{}: {}".format(array_num, array.id
					) for array_num, array in entry.arrays.items()
				) + append_stuff
			array_locations = prepend_stuff + min_delim.join(
				"{}: {} {}-{}".format(
					array_num,
					array.contig,
					array.start,
					array.stop
					) for array_num, array in entry.arrays.items()
				) + append_stuff
			repeat_sequences = min_delim.join(
				"{}: {}".format(
					array_num,
					Counter(array.repeats).most_common(1)[0][0]
					) for array_num, array in entry.arrays.items()
				)
			array_crispr_types = prepend_stuff + min_delim.join(
				"{}: {}".format(
					array_num, array.repeat_id
					) for array_num, array in entry.arrays.items()
				) + append_stuff

		else: # If txt, array number just complicates downstream work
			sequence_id = entry.accession
			has_crispr = str(entry.has_crispr)
			array_count = str(entry.array_count)
			spacers = min_delim.join(
				"{}".format(
					" ".join(array.spacers)
					) for array in entry.arrays.values()
				)
			spacer_ids = min_delim.join(
				"{}".format(
					" ".join(array.spacer_ids)
					) for array in entry.arrays.values()
				)
			array_ids = min_delim.join(
				"{}".format(
					array.id
					) for array in entry.arrays.values()
				)
			array_locations = min_delim.join(
				"{} {}-{}".format(
					array.contig,
					array.start,
					array.stop
					) for array in entry.arrays.values()
				)
			repeat_sequences = min_delim.join(
				"{}".format(
					Counter(array.repeats).most_common(1)[0][0]
					) for array in entry.arrays.values()
				)
			array_crispr_types = min_delim.join(
				"{}".format(
					array.repeat_id
					) for array in entry.arrays.values()
				)

		line = [
		sequence_id,
		has_crispr,
		array_count,
		spacers,
		spacer_ids,
		array_ids,
		array_locations,
		repeat_sequences,
		array_crispr_types,
		]
		
		outcontents += maj_delim.join(line)+"\n"

	with open(outfile, 'w') as fout:
		fout.write(outcontents)
		

def write_array_reps_file(all_assemblies, outdir):
	array_reps_dict = defaultdict(list)
	for assembly in all_assemblies:
		acc = assembly.accession
		for array in assembly.arrays.values():
			array_reps_dict[array.id].append(acc)

	outcontents = []
	for arr, accs in array_reps_dict.items():
		outcontents.append("{}\t{}".format(arr, " ".join(accs)))

	with open(outdir+"Array_representatives.txt", 'w') as fout:
		fout.write("\n".join(outcontents))


def write_CRISPR_files(
	all_assemblies,
	spacer_id_dict,
	array_id_dict,
	outdir):
	""" Write output files from CRISPR identification scripts
	"""
	# Spacers in fasta format
	write_spacer_fasta(spacer_id_dict, outdir+"CRISPR_spacers.fna")
	
	# Array IDs and spacers
	write_array_id_seq_file(spacer_id_dict, array_id_dict, outdir)

	# Array representatives
	write_array_reps_file(all_assemblies, outdir)

	# CSV summary table
	write_cr_sum_tabs(all_assemblies, outdir+"CRISPR_summary_table.csv")
	# TXT summary table
	write_cr_sum_tabs(all_assemblies, outdir+"CRISPR_summary_table.txt")


def read_assembly_list_file(filename):
	assemblies = []
	with open(filename, 'r') as fin:
		for line in fin:
			assemblies.append(line.strip())

	return assemblies
