#!/usr/bin/env python

import sys
import argparse
import random


class Array():
	"""Store information about arrays.

	Attributes:
		parent (Array):
		  The array from which this array was derived.
		age_weight (int):
		  Age of the array.
	"""
	def __init__(self, parent=None, age_weight=10):
		self.parent = parent
		self.age_weight = age_weight
		self.spacers = []
		

def cmdline_args():

	p = argparse.ArgumentParser(
		description="evolve CRISPR arrays in silico. Starts with a single \
		array and performs a simple in silico evolution process."
		)

	required = p.add_argument_group('Required arguments')
	required.add_argument(
		"-n", "--num-events", required = True, type=int, metavar="",
		help="How many events should be allowed to occur before the \
		simulation ends?"
		)
	required.add_argument(
		"-o", "--outdir", required = True, type=str, metavar="",
		help="Directory in which output files should be written"
		)

	parameters = p.add_argument_group('Evolution parameters', 
		"Specify the relative frequencies with which different events \
		should occur.")
	parameters.add_argument(
		"-i", "--initial-length", type=int, default=5, metavar="",
		help="Length of the starting array")
	parameters.add_argument(
		"-a", "--acquisition", type=int, default=80, metavar="",
		help="Relative frequency of spacer acquisitions")
	parameters.add_argument(
		"-t", "--trailer-loss", type=int, default=15, metavar="",
		help="Relative frequency of trailer spacer decay")
	parameters.add_argument(
		"-d", "--deletion", type=int, default=5, metavar="",
		help="Relative frequency of deletions of one or more spacers")

	return p.parse_args()


def main(args):

	spacer_n = 1
	init_array = Array()
	for _ in range(args.initial_length):
		init_array.spacers.append(spacer_n)
		spacer_n+=1
	



	


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
