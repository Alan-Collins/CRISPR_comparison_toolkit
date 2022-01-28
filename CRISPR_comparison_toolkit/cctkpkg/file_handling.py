import os
import json

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

