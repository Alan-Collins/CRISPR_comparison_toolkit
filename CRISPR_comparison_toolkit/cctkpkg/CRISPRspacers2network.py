#!/usr/bin/env python3

import sys
import argparse

class Assembly():	
	""" Class to process lines from CRISPR_summary_table.csv

	Takes file output by CRISPR identification scripts and reads array 
	information from them.
	
	Attributes:
	  ass_id (str):
	    Name of assembly or sequence record.
	  present (bool):
	    Were CRISPR arrays found in this assembly?
	  count (int):
	    How many CRISPR arrays were found?
	  array_ids (dict):
	    spacers found in each array of format
	    {'Array number' : ['IDs', 'of', 'spacers', ...]}
	  types (dict):
	    The CRISPR subtype that each array was characterized as. Format:
	    {'Array number' : 'subtype'}
	"""
	def __init__(self, ass):
		self.ass_id = ass[0]
		self.present = True if ass[1] == "True" else False
		self.count = int(ass[2])
		self.array_ids = {}
		self.types = {}
		if self.present:
			for i in range(self.count):
				self.arrays[i+1] = ass[4].split('\t')[i].split()[1:]
				self.array_ids[i+1] = ass[5].split('\t')[i].split()[1][:-1]
				self.types[i+1] = ass[9].split('\t')[i].split()[1]



class CRISPR_info():
	"""Key information about CRISPR arrays.
	
	Attributes:
	  array_id (str):
	    Identifier of array.
	  spacers (list):
	    An integer count of the eggs we have laid.
	"""
	def __init__(self, assembly, array_num):
		self.array_id = assembly.array_ids[array_num]
		self.spacers = assembly.arrays[array_num]
		self.type = assembly.types[array_num]


def build_parser(parser):
	parser.add_argument(
		"-i", "--input", required = True,
		help="Input file. CRISPR_summary_table.csv produced by one of the "
		"CRISPR identification scripts."
		)
	parser.add_argument(
		"-n", "--network-file", default="Array_network",
		help="Output network file. Path to location to write network " 
		"representation of array relationships."
		)
	parser.add_argument(
		"-r", "--reps-file", default="Array_representatives",
		help="Output array representatives file. "
		"Path to location to write a list of which assemblies have each array."
		)

	return parser


def main(args):

	arrays = []
	included_arrays = []
	repped_ass_ars = {}

	with open(args.input, 'r') as CRISPR_infile:
		for line in CRISPR_infile.readlines()[1:]:
			ass_info = Assembly(line.split(','))
			if not ass_info.present:
				continue
			for i in range(ass_info.count):
				array_num = i + 1
				if ass_info.array_ids[array_num] not in included_arrays:
					arrays.append(CRISPR_info(ass_info, array_num))
					included_arrays.append(ass_info.array_ids[array_num])
					repped_ass_ars[
						ass_info.array_ids[array_num]] = [ass_info.ass_id]
				else:
					repped_ass_ars[
						ass_info.array_ids[array_num]].append(ass_info.ass_id)

	edges = []

	for i in range(len(arrays)):
		for j in range(i+1, len(arrays)):
			shared_spacers = list(
				set(arrays[i].spacers) & set(arrays[j].spacers))
			if len(shared_spacers) != 0:
				edges.append("\t".join(
					[arrays[i].array_id,
					str(len(shared_spacers)),
					arrays[j].array_id,
					arrays[i].type,
					arrays[j].type,
					str(len(arrays[i].spacers)),
					str(len(arrays[j].spacers)),
					str(len(set(arrays[i].spacers + arrays[j].spacers))),
					str(len(shared_spacers)/len(
						set(arrays[i].spacers + arrays[j].spacers)))
					]))

	with open(args.rep_file, 'w+') as fout:
		fout.write("Array_ID\tAssemblies_with_array\n")
		for k,v in repped_ass_ars.items():
			fout.write("{}\t{}\n".format(k, " ".join(v)))

	with open(args.network_file, 'w+') as fout:
		fout.write("Source_ID\tShared_Spacers\tTarget_ID\tSource_CRISPR_Type\t"
			"Target_CRISPR_Type\tLen_Source_Array\tLen_Target_Array\t"
			"Total_Unique_Spacers\tProportion_Unique_Spacers_Shared\n")
		fout.write("\n".join(edges))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Create a network representation of the relationships \
		between arrays."
		)
	parser = build_parser(parser)

	if len(sys.argv) == 1:
		parser.parse_args(['--help'])
	else:
		args = parser.parse_args()

	main(args)	
