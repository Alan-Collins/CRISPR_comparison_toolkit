#!/usr/bin/env python

import sys
import argparse
import random

import dendropy


class Array():
	"""Store information about arrays.

	Attributes:
		name (str):
		  Array ID
		parent (Array):
		  The array from which this array was derived.
		spacers (list):
		  The spacers in this array.
	"""
	def __init__(self, name, parent=None):
		self.name = str(name)
		self.parent = parent
		if self.parent:
			self.spacers = [int(i) for i in parent.spacers]
		else:
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


def tick(active_arrays, tree, tree_namespace, spacer_n, array_name, events_dict,
	all_arrays, loss_rate):
	"""Choose an event and return the resulting spacer.
	
	Picks an existing array to modify and an event at random. Creates a
	modified copy of the chosen array according to the event chosen and
	returns the modified copy of the array

	Args:
	  active_arrays (list of Array):
		List of Arrays that are still available to use as source arrays.
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
	    Percent rate at which to remove old arrays from the 
	    active_arrays after they are the source of a new array

	Returns:
	  A modified form of the chosen array as an Array class instance and
	  other updates objects.
	"""

	source_array = random.choice(active_arrays)
	event = events_dict[random.randint(1,len(events_dict))]

	# Don't allow deletion if array is only 1 spacer long

	if len(source_array.spacers) == 1:
		while event == "Deletion":
			event = events_dict[random.randint(1,len(events_dict))]

	array = Array(name=array_name, parent=source_array)
	array_name+=1

	if event == "Acquisition":
		array, spacer_n = do_acquisition(array, spacer_n)
	elif event == "Trailer_loss":
		array = do_trailer_loss(array)
	else:
		array = do_deletion(array)

	# Add new array to active_arrays and all_arrays lists
	active_arrays.append(array)
	all_arrays.append(array)

	# Check whether to remove source array from active_arrays
	x = random.randint(1,100)
	if x <= loss_rate:
		del_idx = active_arrays.index(source_array)
		del active_arrays[del_idx]

	# Add array relationships to tree

	new_node = dendropy.Node()
	new_node.taxon = tree_namespace.get_taxon(array.name)

	tree.find_node_with_taxon_label(source_array.name).add_child(new_node)


	return active_arrays, tree, spacer_n, array_name, all_arrays


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

	# To ensure a deletion takes place, but not the whole array
	while b == a or (a in [0, array_len] and b in [0, array_len]):
		b = pick_normal_index(array_len, mean, stddev)

	if a > b:
		a,b = b,a

	spacers_to_del = [i for i in range(a,b)]

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
	array_name = 0
	active_arrays = []
	events_dict = {}

	init_array = Array(name=array_name)
	array_name+=1

	for _ in range(args.initial_length):
		init_array.spacers.insert(0,spacer_n)
		spacer_n+=1
	all_arrays = [init_array]
	active_arrays = [init_array]
	

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
		[str(i) for i in range(args.num_events+1)])
	tree = dendropy.Tree(taxon_namespace=tree_namespace)

	first_node = dendropy.Node()
	first_node.taxon = tree_namespace.get_taxon('0')
	tree.seed_node.add_child(first_node)

	for _ in range(args.num_events):
		active_arrays, tree, spacer_n, array_name, all_arrays = tick(
			active_arrays, tree, tree_namespace, spacer_n, array_name, events_dict,
			all_arrays, args.loss_rate)

	for a in active_arrays:
		print(str(a.name) + '\t' + ' '.join([str(i) for i in a.spacers]))


	print(tree.as_ascii_plot(show_internal_node_labels=True))



	


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
