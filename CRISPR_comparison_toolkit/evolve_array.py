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
	def __init__(self, parent=None):
		self.parent = parent
		if self.parent:
			self.age_weight = int(parent.age_weight*1.2)
			self.spacers = [int(i) for i in parent.spacers]
		else:
			self.age_weight = 10
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

	array = Array(parent=source_array)

	if event == "Acquisition":
		array, spacer_n = do_acquisition(array, spacer_n)
	elif event == "Trailer_loss":
		array = do_trailer_loss(array)
	else:
		array = do_deletion(array)

	print(array.spacers)
	sys.exit()

	return array, spacer_n


def do_acquisition(array, spacer_n):
	array.spacers.append(spacer_n)
	spacer_n+=1

	return array, spacer_n


def do_trailer_loss(array):
	del array.spacers[-1]

	return array
	
def do_deletion(array):
	# Define parameters for a normal distribution
	# This distribution will be used to select random spacers that
	# tend to be close to the middle of the array
	array_len = len(array.spacers)
	# Define the central point of the array
	mean = array_len/2

	# Set the level of variation from the middle of the array
	stddev = array_len/6

	a = pick_normal_index(array_len, mean, stddev)
	b = pick_normal_index(array_len, mean, stddev)

	# To ensure a deletion takes place
	while b == a:
		b = pick_normal_index(array_len, mean, stddev)

	if a > b:
		a,b = b,a

	spacers_to_del = [i for i in range(a,b+1)]

	print(a,b)

	array.spacers = [
		i for n,i in enumerate(array.spacers) if n not in spacers_to_del]

	return array


def pick_normal_index(max, mean, stddev):

	while True: # Keep going until a valid index is chosen
		idx = int(random.normalvariate(mean, stddev))
		if 0 <= idx <= max:
			return idx


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

	x, spacer_n = tick(age_dict, spacer_n, events_dict)

	# print(vars(x))





	


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
