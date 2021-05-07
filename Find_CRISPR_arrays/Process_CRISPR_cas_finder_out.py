#!/usr/bin/env python3

# AUTHOR		:	ALAN COLLINS
# VERSION		:	v1
# DATE			:	2021/1/11
# DESCRIPTION	:	Script to process json files output by CRISPR-cas finder and summarize them in various human-readable files.

import sys
import argparse
import json
import os
import xml.etree.ElementTree as Et





parser = argparse.ArgumentParser(
	description="Script to process json files output by CRISPR-cas finder and summarize them in various human-readable files.",
	formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument(
	"-i",  dest="indirs", required = False, nargs="+",
	help="Specify input directories containing output files that were produced by CRISPR-cas finder."
	)
parser.add_argument(
	"-o",  dest="outdir", required = False,
	help="Specify output directory into which you want the script to place its output files. Default is to output in the current working directory"
	)
parser.add_argument(
    "-r", dest="repeats_file", required = False,
    help="FASTA format file containing the CRISPR repeats you want to look for. If you don't provide one, this script will look for type 1F, 1E, and 1C repeats from Pseudomonas aeruginosa. These are used to classify arrays and orient them consistently with one another."
    )
# parser.add_argument(
# 	"-p",  dest="outprefix", required = False,
# 	help="Specify a prefix string you want at the beginning of output file names."
# 	)

def get_repeat_info(CRISPR_types_dict, repeat):

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

def fasta_to_dict(FASTA_file):
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


def hamming(string1, string2):
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


def rev_comp(string):
	rev_str = ''
	rev_comp_lookup = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c"}
	for i in reversed(string):
		if i in "ATCGatcg":
			rev_str += rev_comp_lookup[i]
		else:
			rev_str += i
	return rev_str
	

args = parser.parse_args(sys.argv[1:])


if args.repeats_file:
	CRISPR_types_dict = fasta_to_dict(args.repeats_file)
else:
	CRISPR_types_dict = { #Specify the desired orientation and type of all repeats you are interested in here
	'1E':'GTGTTCCCCACGGGTGTGGGGATGAACCG',
	'1F':'GTTCACTGCCGTGTAGGCAGCTAAGAAA',
	'1C':'GTCGCGCCCCGCACGGGCGCGTGGATTGAAAC'
	}

### Read in CRISPR system Cas protein specifications (stolen from Macsyfinder packaged with CRISPRCasFinder) Files must be stored in ./Data_files/CRISPR_system_defs/ relative to this script's location.

Cas_dict = {}

path_prefix = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
definitions_path = path_prefix + "Data_files/CRISPR_system_defs/"

for file in os.listdir(definitions_path):
	CR_type = file.strip('.xml')
	Cas_dict[CR_type] = {'n_mandatory' : 0}

	tree = Et.parse(definitions_path + file)
	root = tree.getroot()

	for child in root:
		Cas_dict[CR_type][child.attrib["name"]] =  child.attrib["presence"]
		if child.attrib["presence"] == 'mandatory':
			Cas_dict[CR_type]['n_mandatory'] += 1


results = {}
Cas_results = {}
spacers = {}
nspacers_seen = 1
arrays = {}
arrays_encoded = {}
narrays_seen = 1


outdir = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir

for indir in args.indirs:
	indir = indir + '/' if indir[-1] != '/' else indir
	injson = indir + 'result.json'
	incas = indir + 'rawCas.fna'
	with open(injson, 'r') as fin:
		result = json.loads(fin.read())

	isolate_id = indir.split('/')[-2].split('_')[1]

	if isolate_id not in results.keys():
		results[isolate_id] = {}
		array_num = 0
	for i in result['Sequences']:
		if len(i['Crisprs']) > 0:
			for j in i['Crisprs']:
				if int(j['Evidence_Level']) > 1:
					array_num += 1
					results[isolate_id][array_num] = {
						'Contig' : "_".join(j["Name"].split("_")[:-1]), 
						'Location' : str(str(j["Start"]) + " : " + str(j["End"])), 
						'Number_Spacers' : j["Spacers"], 
						'DR_Consensus' : j["DR_Consensus"], 
						'DR_Type' : get_repeat_info(CRISPR_types_dict, j["DR_Consensus"])[0],
						'DR_Hamming_Dist' : get_repeat_info(CRISPR_types_dict, j["DR_Consensus"])[1],
						"Array_encoded" : '',
						'Spacer' : [], 
						'Spacer_encoded' : [], 
						'DR' : [],
						'Leader' : ''
						
						}
					if j['Potential_Orientation'] == '+':
						for k in j['Regions']:
							if 'FLANK' in k['Type']:
								if k['Leader'] == 1:
									results[isolate_id][array_num]['Leader'] = k['Sequence']
							else:	
								results[isolate_id][array_num][k['Type']].append(k['Sequence']) # Adds the sequence to either the DR or spacer list for that array depending on the 'Type' in the json file.
								if k['Type'] == 'Spacer':
									if k['Sequence'] in spacers.keys():
										results[isolate_id][array_num]['Spacer_encoded'].append(str(spacers[k['Sequence']]))
									else:
										spacers[k['Sequence']] = nspacers_seen
										results[isolate_id][array_num]['Spacer_encoded'].append(str(nspacers_seen))
										nspacers_seen += 1

					else:
						results[isolate_id][array_num]['DR_Consensus'] = rev_comp(j["DR_Consensus"])
						for k in reversed(j['Regions']):
							if 'FLANK' in k['Type']:
								if k['Leader'] == 1:
									results[isolate_id][array_num]['Leader'] = rev_comp(k['Sequence'])
							else:	
								results[isolate_id][array_num][k['Type']].append(rev_comp(k['Sequence']))
								if k['Type'] == 'Spacer':
									if rev_comp(k['Sequence']) in spacers.keys():
										results[isolate_id][array_num]['Spacer_encoded'].append(str(spacers[rev_comp(k['Sequence'])]))
									else:
										spacers[rev_comp(k['Sequence'])] = nspacers_seen
										results[isolate_id][array_num]['Spacer_encoded'].append(str(nspacers_seen))
										nspacers_seen += 1
					# array stuff
					if " ".join(results[isolate_id][array_num]['Spacer']) in arrays.keys():
						results[isolate_id][array_num]["Array_encoded"] = arrays[" ".join(results[isolate_id][array_num]['Spacer'])]
					else:
						arrays[" ".join(results[isolate_id][array_num]['Spacer'])] = str(narrays_seen)
						arrays_encoded[" ".join(results[isolate_id][array_num]['Spacer_encoded'])] = str(narrays_seen)
						results[isolate_id][array_num]["Array_encoded"] = str(narrays_seen)
						narrays_seen += 1
	Cas_results[isolate_id] = {
		'Cas_found' : "FALSE",
		'Cas_types' : [],
		'Cas_completeness' : {},
		'Cas_completeness_binary' : '1', # Will be set to 0 if any identified CRISPR subtypes don't have all mandatory genes.
		'Cas_list' : {}
		}
	Cas_headers = fasta_to_dict(incas).keys()
	if len(Cas_headers) > 0:
		Cas_results[isolate_id]['Cas_found'] = "TRUE"
		for header in Cas_headers:
			elements = header.split('|')
			CR_type = elements[2]
			gene = elements[3].split()[0]
			if CR_type not in Cas_results[isolate_id]['Cas_types']:
				Cas_results[isolate_id]['Cas_types'].append(CR_type)
				if CR_type != "CAS":
					Cas_results[isolate_id]['Cas_completeness'][CR_type] = [0, Cas_dict[CR_type]['n_mandatory']]
					if Cas_dict[CR_type][gene] == 'mandatory':
						Cas_results[isolate_id]['Cas_completeness'][CR_type][0] += 1

				Cas_results[isolate_id]['Cas_list'][CR_type] = [gene]
			else:
				if CR_type != "CAS":
					if gene not in Cas_results[isolate_id]['Cas_list'][CR_type] and Cas_dict[CR_type][gene] == 'mandatory':
						Cas_results[isolate_id]['Cas_completeness'][CR_type][0] += 1
				Cas_results[isolate_id]['Cas_list'][CR_type].append(gene)
		for v in Cas_results[isolate_id]['Cas_completeness'].values():
			if v[0] < v[1]:
				Cas_results[isolate_id]['Cas_completeness_binary'] = '0'

for k,v in Cas_results.items():
	print(k,v)

		

		




CRISPR_sum_to_write = [["Strain,Has_CRISPR,Array_count,Spacers,Spacers_encoded,Array_IDs,Array_locations,Array_sizes,Repeat_sequences,Array_CRISPR_types,Repeat_scores, Leader, Cas_found, Cas_types, Cas_completeness, Cas_list"]]

for k,v in results.items():
	if len(v) == 0:
		CRISPR_sum_to_write.append([k, 'FALSE'])
	else:
		Strain = k
		Has_CRISPR = "TRUE"
		Array_count = str(len(v))
		Spacers = []
		Spacers_encoded = []
		Array_ids = []
		Array_locations = []
		Array_sizes = []
		Repeat_sequences = []
		Array_types = []
		Repeat_scores = []
		Leaders = []
		for array_num, feature_dict in v.items():
			Spacers.append(str(array_num) + ": " + " ".join(feature_dict['Spacer']))
			Spacers_encoded.append(str(array_num) + ": " + " ".join(feature_dict['Spacer_encoded']))
			Array_ids.append("'" + str(array_num) + ": " + str(feature_dict['Array_encoded']) + "'")
			Array_locations.append(str(array_num) + ": " + " ".join(feature_dict['Location']))
			Array_sizes.append("'" + str(array_num) + ": " + str(feature_dict['Number_Spacers']) +"'")
			Repeat_sequences.append(str(array_num) + ": " + feature_dict['DR_Consensus'])
			Array_types.append(str(array_num) + ": " + feature_dict['DR_Type'])
			Repeat_scores.append("'" + str(array_num) + ": " + str(feature_dict['DR_Hamming_Dist']) + "'")
			Leaders.append(str(array_num) + ": " + feature_dict['Leader'])
		Spacers = "\t".join(Spacers)
		Spacers_encoded = "\t".join(Spacers_encoded)
		Array_ids = "\t".join(Array_ids)
		Array_locations = "\t".join(Array_locations)
		Array_sizes = "\t".join(Array_sizes)
		Repeat_sequences = "\t".join(Repeat_sequences)
		Array_types = "\t".join(Array_types)
		Repeat_scores = "\t".join(Repeat_scores)
		Leaders = "\t".join(Leaders)
		Cas_found = Cas_results[isolate_id]['Cas_found']
		Cas_types = " ".join(Cas_results[isolate_id]['Cas_types'])
		Cas_completeness = []
		Cas_list = []
		for crtype,feature in Cas_results[isolate_id]['Cas_completeness'].items():
			Cas_completeness.append(crtype + ': ' + str(feature[0]) + '/' + str(feature[1])) 
			Cas_list.append(crtype + ': ' + " ".join(Cas_results[isolate_id]['Cas_list'][crtype]))
		Cas_completeness = '\t'.join(Cas_completeness)
		Cas_completeness_binary = Cas_results[isolate_id]['Cas_completeness_binary']
		Cas_list = '\t'.join(Cas_list)
		CRISPR_sum_to_write.append([Strain,Has_CRISPR,Array_count,Spacers,Spacers_encoded,Array_ids, Array_locations,Array_sizes,Repeat_sequences,Array_types,Repeat_scores, Leaders, Cas_found, Cas_types, Cas_completeness, Cas_completeness_binary, Cas_list])


# with open(outdir + "CRISPR_summary_table.csv", 'w+') as fout:
# 	fout.write("\n".join([','.join(line) for line in CRISPR_sum_to_write]))

