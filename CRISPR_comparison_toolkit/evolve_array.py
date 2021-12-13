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

	p.add_argument(
		"-s", "--seed", type=int, metavar="", help="Set seed for consistency")

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


def tick(age_dict, spacer_n, events_dict):
	"""Choose an event and return the resulting spacer.
	
	Picks an existing array to modify and an event at random. Creates a
	modified copy of the chosen array according to the event chosen and
	returns the modified copy of the array

	Args:
	  age_dict (dict of Array):
		dict with all existing Array instances.
	  spacer_n (int):
		Next spacer ID to use.
	  events_dict (dict):
	    Dict to look up to which event a selected random number
	    corresponds

	Returns:
	  A modified form of the chosen array as an Array class instance.
	"""
	if args.seed:
		random.seed(args.seed)
	
	source_array = age_dict[random.randint(1,len(age_dict)+1)]
	event = events_dict[random.randint(1,len(events_dict)+1)]

	if event == "Acquisition":
		array = do_acquisition(source_array, spacer_n)
	elif event == "Trailer_loss":
		array = do_trailer_loss(source_array)
	else:
		array = do_deletion(source_array)

	return array


def do_acquisition(source_array, spacer_n):
	pass

def do_trailer_loss(source_array):
	pass
	
def do_deletion(source_array):
	pass


def main(args):

	spacer_n = 1 # Track ID of spacers
	total_age = 1 # Track 
	age_dict = {}
	events_dict = {}

	init_array = Array()
	for _ in range(args.initial_length):
		init_array.spacers.append(spacer_n)
		spacer_n+=1
	
	for _ in range(init_array.age_weight):
		age_dict[total_age] = init_array
		total_age+=1

	i = 1
	for _ in range(args.acquisition):
		events_dict[i] = 'Acquisition'
		i+=1
	for _ in range(args.trailer_loss):
		events_dict[i] = 'Trailer_loss'
		i+=1
	for _ in range(args.deletion):
		events_dict[i] = 'Deletion'
		i+=1

	x = tick(age_dict, spacer_n, events_dict)

	# print(vars(x))





	


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
