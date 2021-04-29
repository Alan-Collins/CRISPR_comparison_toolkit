#!/usr/bin/env python3

import sys
import argparse
import textwrap as _textwrap
import matplotlib.pyplot as plt
from random import shuffle

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
	"""
	Short function for argparse that wraps text properly when printing to terminal
	"""
	def _split_lines(self, text, width):
		text = self._whitespace_matcher.sub(' ', text).strip()
		return _textwrap.wrap(text, width)


parser = argparse.ArgumentParser(
	description="Given CRISPR_summary_table.txt that is output by minced2arrays.py, plot a rarefaction curve of novel spacers for the assemblies represented in that file.",
	formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument(
	"-i", dest="input_table", required = True,
	help="Specify input CRISPR_summary_table.txt that is output by minced2arrays.py."
	)
parser.add_argument(
	"-o", dest="outplot", required = True,
	help="Specify filename for output plot."
	)

args = parser.parse_args(sys.argv[1:])

CR_dict = {}
presence_dict = {}

with open(args.input_table, 'r') as CRin:
	for line in CRin.readlines()[1:]:
		cols = line.split('\t')
		genome = cols[0].replace('_scaffolds', '').replace('_modified', '')
		spacers = [i.split() for i in cols[4].split('|')]
		presence = cols[1]
		CR_dict[genome] = spacers
		presence_dict[genome] = presence


n = []

ngenomes = 0

spacers = []

x = []
y = []

genome_order = list(CR_dict.keys())
shuffle(genome_order)

for k in genome_order:
	v = CR_dict[k]
	if presence_dict[k] != "False":
		ngenomes += 1
		spacers += [item for sublist in v for item in sublist]
		x.append(ngenomes)
		y.append(len(set(spacers)))
plt.scatter(x,y, s=1)
plt.xlabel("Cumulative number of genomes in dataset")
plt.ylabel("Number of unique spacers in dataset")
plt.savefig(args.outplot, dpi=300)
