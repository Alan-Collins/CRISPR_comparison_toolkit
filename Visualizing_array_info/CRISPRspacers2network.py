#!/usr/bin/env python3

# Takes ouput from minced2arrays.py and produces a file describing which assemblies have which arrays and a file describing the network of arrays (nodes) and the spacers they share (edges).

import sys

class assembly():
	""" Class to process lines from CRISPR_summary_table.csv output be minced2arrays.py. Takes line of file split on commas"""
	def __init__(self, ass):
		self.ass_id = ass[0]
		self.present = True if ass[1] == "True" else False
		self.count = int(ass[2])
		self.arrays = {}
		self.array_ids = {}
		self.types = {}
		if self.present:
			for i in range(self.count):
				self.arrays[i+1] = ass[4].split('\t')[i].split()[1:]
				self.array_ids[i+1] = ass[5].split('\t')[i].split()[1][:-1]
				self.types[i+1] = ass[9].split('\t')[i].split()[1]



class CRISPR_info():
	def __init__(self, assembly, array_num):
		self.array_id = assembly.array_ids[array_num]
		self.spacers = assembly.arrays[array_num]
		self.type = assembly.types[array_num]

infile = sys.argv[1]
outfile = sys.argv[2]

arrays = []
included_arrays = []
represented_assembly_arrays = {}

with open(infile, 'r') as CRISPR_infile:
	for line in CRISPR_infile.readlines()[1:]:
		ass_info = assembly(line.split(','))
		if ass_info.present:
			for i in range(ass_info.count):
				array_num = i + 1
				if ass_info.array_ids[array_num] not in included_arrays:
					arrays.append(CRISPR_info(ass_info, array_num))
					included_arrays.append(ass_info.array_ids[array_num])
					represented_assembly_arrays[ass_info.array_ids[array_num]] = [ass_info.ass_id]
				else:
					represented_assembly_arrays[ass_info.array_ids[array_num]].append(ass_info.ass_id)

edges = []

for i in range(len(arrays)):
	for j in range(i+1, len(arrays)):
		shared_spacers = list(set(arrays[i].spacers) & set(arrays[j].spacers))
		if len(shared_spacers) != 0:
			edges.append("\t".join([arrays[i].array_id, str(len(shared_spacers)), arrays[j].array_id, arrays[i].type, arrays[j].type, str(len(arrays[i].spacers)), str(len(arrays[j].spacers)), str(len(set(arrays[i].spacers + arrays[j].spacers))), str(len(shared_spacers)/len(set(arrays[i].spacers + arrays[j].spacers)))]))

with open("Array_ID_representatives.txt", 'w+') as id_file:
	id_file.write("Array_ID\tAssemblies_with_array\n")
	for k,v in represented_assembly_arrays.items():
		id_file.write("%s\t%s\n" %(k, " ".join(v)))

with open(outfile, 'w+') as network_file:
	network_file.write("Source_ID\tShared_Spacers\tTarget_ID\tSource_CRISPR_Type\tTarget_CRISPR_Type\tLen_Source_Array\tLen_Target_Array\tTotal_Unique_Spacers\tProportion_Unique_Spacers_Shared\n")
	network_file.write("\n".join(edges))


	

