#!/usr/bin/env python

import sys
import argparse
import random

import dendropy


class Array():
	"""Store information about arrays.

	Attributes:
		parent (Array):
		  The array from which this array was derived.
		age_weight (int):
		  Age of the array.
	"""
	def __init__(self, name=0, parent=None):
		self.name = str(name)
		self.parent = parent
		if self.parent:
			self.age_weight = int(parent.age_weight*1.5)
			self.spacers = [int(i) for i in parent.spacers]
		else:
			self.age_weight = 2
			self.spacers = []
		self.age_dict_idxs = []
		

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
		"-a", "--acquisition", type=int, default=75, metavar="",
		help="Relative frequency of spacer acquisitions")
	parameters.add_argument(
		"-t", "--trailer-loss", type=int, default=15, metavar="",
		help="Relative frequency of trailer spacer decay")
	parameters.add_argument(
		"-d", "--deletion", type=int, default=10, metavar="",
		help="Relative frequency of deletions of one or more spacers")
	parameters.add_argument(
		"-l", "--loss-rate", type=int, default=50, metavar="",
		help="Rate (percent) at which arrays are lost after spawning a \
		descendant")

	return p.parse_args()


def tick(age_dict, tree, tree_namespace, spacer_n, array_name, events_dict,
	all_arrays, loss_rate):
	"""Choose an event and return the resulting spacer.
	
	Picks an existing array to modify and an event at random. Creates a
	modified copy of the chosen array according to the event chosen and
	returns the modified copy of the array

	Args:
	  age_dict (dict of Array):
		dict with all existing Array instances.
	  tree (Dendropy.Tree):
	    Tree instance containing array histories.
	  tree_namespace (Dendropy.TaxonNamespace):
	    namespace for tree nodes.
	  spacer_n (int):
		Next spacer ID to use.
	  array_name (int):
	    The ID to assign to newly create Array instances.
	  events_dict (dict):
	    Dict to look up to which event a selected random number
	    corresponds
	  all_arrays (list):
	    list of all Array instances
	  loss_rate (int):
	    Percent rate at which to remove old arrays from the age_dict 
	    after they are the source of a new array

	Returns:
	  A modified form of the chosen array as an Array class instance and
	  other updates objects.
	"""

	source_array = age_dict[random.choice([i for i in age_dict.keys()])]
	event = events_dict[random.randint(1,len(events_dict))]

	array = Array(name=array_name, parent=source_array)
	array_name+=1

	if event == "Acquisition":
		array, spacer_n = do_acquisition(array, spacer_n)
	elif event == "Trailer_loss":
		array = do_trailer_loss(array)
	else:
		array = do_deletion(array)

	# Add new array to age_dict and all_arrays
	max_age = max([i for i in age_dict.keys()])
	for i in range(max_age+1, max_age+array.age_weight+1):
		age_dict[i] = array
		array.age_dict_idxs.append(i)

	all_arrays.append(array)

	# Check whether to remove source array from age_dict
	x = random.randint(1,100)
	if x <= loss_rate:
		for i in source_array.age_dict_idxs:
			del age_dict[i]

	# Add array relationships to tree

	new_node = dendropy.Node()
	new_node.taxon = tree_namespace.get_taxon(array.name)

	tree.find_node_with_taxon_label(source_array.name).add_child(new_node)


	return age_dict, tree, spacer_n, array_name, all_arrays


def do_acquisition(array, spacer_n):
	array.spacers.insert(0,spacer_n)
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
	stddev = array_len/4

	a = pick_normal_index(array_len, mean, stddev)
	b = pick_normal_index(array_len, mean, stddev)

	# To ensure a deletion takes place
	while b == a:
		b = pick_normal_index(array_len, mean, stddev)

	if a > b:
		a,b = b,a

	spacers_to_del = [i for i in range(a,b+1)]

	array.spacers = [
		i for n,i in enumerate(array.spacers) if n not in spacers_to_del]

	return array


def pick_normal_index(max, mean, stddev):

	while True: # Keep going until a valid index is chosen
		idx = int(random.normalvariate(mean, stddev))
		if 0 <= idx <= max:
			return idx


def main(args):

	if args.seed:
		random.seed(args.seed)	

	spacer_n = 1 # Track ID of spacers
	total_age = 1 # Track 
	array_name = 1
	age_dict = {}
	events_dict = {}

	init_array = Array(name=array_name)
	for _ in range(args.initial_length):
		init_array.spacers.insert(0,spacer_n)
		spacer_n+=1
	all_arrays = [init_array]
	
	for _ in range(init_array.age_weight):
		age_dict[total_age] = init_array
		init_array.age_dict_idxs.append(total_age)
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

	tree_namespace = dendropy.TaxonNamespace(
		[str(i) for i in range(1, args.num_events+1)])
	tree = dendropy.Tree(taxon_namespace=tree_namespace)

	first_node = dendropy.Node()
	first_node.taxon = tree_namespace.get_taxon(1)
	tree.seed_node.add_child(first_node)

	for _ in range(args.num_events):
		age_dict, tree, spacer_n, array_name, all_arrays = tick(
			age_dict, tree, tree_namespace, spacer_n, array_name, events_dict,
			all_arrays, args.loss_rate)

	for a in all_arrays:
		print(str(a.name) + '\t' + ' '.join([str(i) for i in a.spacers]))


	print(tree.as_ascii_plot(show_internal_node_labels=True))



	


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
