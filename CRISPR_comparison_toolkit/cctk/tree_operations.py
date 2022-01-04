from string import ascii_lowercase
from math import ceil, log
from itertools import product
import matplotlib.pyplot as plt


def create_internal_node_ids(n_leaves, prefix="", chars="letters"):
	"""Make list of internal node IDs for tree.

	given the number of leaf nodes in a tree, generates enough IDs to
	uniquely identify all internal nodes of a bifurcating tree.

	Args:
	  n_leaves (int): 
	  	The number of leaf nodes in the tree.
	  prefix (str, optional): 
	  	Prefix with which all node IDs should begin.
	  chars (str, optional): 
	  	Options "letters", "numbers". The kind of character that should 
	  	be used to identify each node.
	Returns:
	  list of str: 
	  	List of unique node IDs.

	Raises:
	  ValueError: If n_leaves is less than 2.
	  TypeError: If n_leaves is not int.
	  TypeError: If prefix is not str.
	  ValueError: If chars is not either "letters" or "numbers".
	"""
	if n_leaves < 2:
		raise ValueError("n_leaves must be 2 or more.")
	if type(n_leaves) is not int:
		raise TypeError(
			"n_leaves must be int, not {}.".format(type(n_leaves).__name__))
	if type(prefix) is not str:
		raise TypeError(
			"prefix must be str, not {}.".format(type(prefix).__name__))
	if chars not in ["letters", "numbers"]:
		raise ValueError('chars must be either "letters" or "numbers".')

	if chars == "letters":
		node_ids = [prefix + i for i in ascii_lowercase]
		# Number of internal nodes in tree is n-1 
		# so only need more than 26 if n >= 28
		# Numer of characters needed to generate N combinations 
		# is ceil(log26(N))
		if n_leaves > 27: 
			node_ids += [prefix + "".join(i) for i in product(
				ascii_lowercase, repeat=(ceil(log(n_leaves, 26))))]
		# End up with more IDs than needed so slice appropriate number.
		node_ids = node_ids[:n_leaves]
	else:
		# Start counting at 1 for general audience.
		node_ids = [prefix+str(i) for i in range(1, n_leaves+1)]
	return node_ids


def scale_branches(tree, max_len=10):
	"""Scales tree branch lengths linearly up to a user-defined maximum.
	
	Args:
	  tree (dendropy.Tree instance):
		Tree to scale
	  max_len (int or float):
		Value to set the longest branch to

	Returns:
	  tree:
	  	dendropy.Tree instance with scaled branch lengths.

	Raises:
	  TypeError: If tree is not a dendropy Tree class instance
	  ValueError: If max_len is not >0
	  TypeError: If max_len is not an int of float.
	"""
	edge_lens = []
	for node in tree:
		if node.edge_length != None:
			edge_lens.append(node.edge_length)

	max_branch = max(edge_lens)
	scale = max_len/max_branch

	for node in tree:
		if node.edge_length != None:
			new_branch = node.edge_length * scale
			# Handle ints vs floats
			if new_branch%1 == 0:
				node.edge_length = int(new_branch)
			else:
				node.edge_length = new_branch
		else:
			node.edge_length = 0
	
	 
	return tree


def draw_branches(tree, node_locs, ax, branch_lengths=True, brlen_scale=0.5):

	for name, location in node_locs.items():
		
		# Add label first
		x, y = location
		label_color = "#000000"
		ax.text(x-0.4, y, name, ha='right', va='center_baseline', fontsize=15,
			color=label_color)
		
		# then add branches
		node = tree.find_node_with_taxon_label(name)

		# First add branch lengths if user desires

		if branch_lengths:
			if node.edge_length != 0:
				ax.text(x + (node.edge_length/2)*brlen_scale, y-0.6,
					node.edge_length, ha='center', va='top', fontsize=12)
		
		# Draw branch from this node to its parent
		
		# identify x location for parent depth if the parent is not seed
		if node.parent_node != tree.seed_node:
			x2 = node_locs[node.parent_node.taxon.label][0]
		else:
			x2 = x + brlen_scale

		ax.plot([x, x2], [y, y], color='black', linewidth = 1,
			solid_capstyle="butt")

		# If internal node, draw line connecting lines to children
		if node.is_internal():
			child_locations = [
				node_locs[i.taxon.label] for i in node.child_nodes()]
			# Draw line between the highest and smallest child y value
			y1 = max([i[1] for i in child_locations])
			y2 = min([i[1] for i in child_locations])

			ax.plot([x, x], [y1, y2], color='black', linewidth = 1,
				solid_capstyle="butt")

	return ax


def yield_nodes(node):
	""" Work through a tree or subtree, yield the taxon labels of nodes
	"""
	if node.is_internal():
		# decide when to insert this node among its children
		# Should be about half way through its list of child nodes
		label_loc = int(len(node.child_nodes())/2)-1
		for i, n in enumerate(node.child_nodes()):
			for x in yield_nodes(n):
				yield x
			if i == label_loc:
				yield node.taxon.label

	else:
		yield node.taxon.label


def find_node_locs(tree, brlen_scale=0.5, branch_spacing=2.3):

	node_list = [n for n in yield_nodes(tree.seed_node.child_nodes()[0])]

	node_locs = {}

	# Prepare tree for root distance queries
	tree.calc_node_root_distances(
		return_leaf_distances_only=False)

	for i, n in enumerate(node_list):
		node = tree.find_node_with_taxon_label(n)
		x_position = -node.root_distance*brlen_scale
		y_position = i*branch_spacing
		node_locs[n] = (x_position, y_position)

	return node_locs


def plot_tree_temp(tree, array_dict, filename, spacer_cols_dict, 
	branch_lengths=False, emphasize_diffs=False, dpi=600, 
	no_align_cartoons=False, no_align_labels=False, fade_ancestral=False):

	node_locs = find_node_locs(tree)
	print(node_locs)
	
	fig, ax = plt.subplots()

	ax = draw_branches(tree, node_locs, ax)

	plt.axis('off')
	plt.tight_layout()

	plt.savefig(filename)


def plot_tree(tree, array_dict, filename, spacer_cols_dict, 
	branch_lengths=False, emphasize_diffs=False, dpi=600, 
	no_align_cartoons=False, no_align_labels=False, fade_ancestral=False):
	""" Plot dendropy Tree with array cartoons
	Args:
	  tree (dendopy Tree class instance):
	    The tree you want to plot.
	  array_dict (dict):
	    Dict of Array class instances with information about the nodes
	    of your tree.
	  filename (str):
	    Path to the file you want created with this plot.
	  spacer_cols_dict (dict):
	    Dict describing which colours have been assigned as the fill and
	    outline colours for spacers.
	  branch_lengths (bool):
	    Should branch lengths be labeled?
	  emphasize_diffs (bool):
	    Should annotations be added to highlight indels and acquisitions
	    in spacer cartoons?
	  dpi (int):
	    The resolution of the output plot
	  no_align_cartoons (bool):
	    Should cartoons of arrays be aligned such that their
	    trailer-most spacers have the same x position? ### FIX
	  no_align_labels (bool):
	    Should labels of arrays be aligned? If not they are placed next
	    to corresponding node. algin_cartoons must also be True for this
	    setting to work.
	  fade_ancestral (bool):
	    Should ancestral array cartoons be made slightly transparent to
	    de-emphasize them?
	"""

	# Find tree dimensions
	tree_width = max(tree.calc_node_root_distances(return_leaf_distances_only=True))
	tree_height = len(tree.nodes())

	#Find deepest node so plotting can start at the most distant leaf
	for node in tree.postorder_node_iter():
		if node.root_distance == tree_width:
			start_node = node

	node_locs = {} # Store where each node is located to draw lines to it.

	start_position = [0,0]
	node_locs[start_node.taxon.label] = start_position

	
	# dim_x, dim_y = 10, 10

	# hscale = (dim_x+1)/(tree_width + max([len(array.spacers) for array in array_dict.values()])) # Factor to scale all branch lengths to fit them in the plot
	# vscale = (dim_y+1)/tree_height
	branch_spacing = 2.3 # Space between branches in tree.
	spacer_width = 0.3 # Thickness off spacer cartoons.
	outline = 3 # Thickness of spacer_outline
	spacing = outline*0.15
	spacer_size = 2.5
	brlen_scale = 0.5
	label_pad = max([len(i) for i in array_dict.keys()] + [5]) # adds 1 unit of space per character in longest array ID. Min space = 5 units

	if not no_align_cartoons or not no_align_labels:
		dash_shift = label_pad
		if not no_align_labels:
			dash_shift = 0

	max_depth = tree_width

	fig, ax = plt.subplots()

	# fig.set_size_inches(dim_x, dim_y)

	node = start_node
	highest_y = 0

	nodes_to_revisit = {} # Store nodes with subtrees and the y value to start at for them

	while True: # Keep going until reaching the root triggers a break.

		first_node = node
		if not first_node.taxon.label in node_locs.keys():
			node_locs[first_node.taxon.label] = position
		if len(node.sibling_nodes())==1:
			if node.taxon.label not in nodes_to_revisit.keys():
				second_node = node.sibling_nodes()[0]
			else:
				del nodes_to_revisit[node.taxon.label] # Remove that from the list as we've finished its tree
				if len(nodes_to_revisit) == 0:
					break
				else:
					node = tree.find_node_with_taxon_label(list(nodes_to_revisit.keys())[0]) # start the next node to revisit
					first_node = node.leaf_nodes()[0]# start drawing from a random leaf in the subtree
					second_node = first_node.sibling_nodes()[0]
					highest_y = y2 = nodes_to_revisit[node.taxon.label] # Set the y position reserved for this subtree
					second_node = first_node.sibling_nodes()[0]
					node_locs[first_node.taxon.label] = [(max_depth-first_node.root_distance)*brlen_scale,y2]
				
				
				


		else: # We've reached the root. Check if any subtrees need to be figured out.
			if len(nodes_to_revisit) == 0: # No subtrees. Exit the loop
				break
			else:
				node = tree.find_node_with_taxon_label(list(nodes_to_revisit.keys())[0]) # start the first node to revisit
				first_node = node.child_nodes()[0].leaf_nodes()[0] # start drawing from a random leaf in the subtree
				second_node = first_node.sibling_nodes()[0]
				highest_y = y2 = nodes_to_revisit[node.taxon.label] # Set the y position reserved for this subtree
				node_locs[first_node.taxon.label] = [(max_depth-first_node.root_distance)*brlen_scale,y2]
		

		# figure out first branch location
		
		x1 = node_locs[first_node.taxon.label][0]
		x2 = node_locs[first_node.taxon.label][0] + first_node.edge_length*brlen_scale 
		y1 = node_locs[first_node.taxon.label][1]
		y2 = node_locs[first_node.taxon.label][1]

		highest_y = max([highest_y,y2])
		
		num_leaves = len(second_node.leaf_nodes()) # Figure out how much space is needed based on the number of leaves below this node
		num_internal = len([i for i in second_node.levelorder_iter(lambda x: x.is_internal())])

		if num_internal > 0:  
			num_internal += num_leaves # - 1 # Counts self so need to subtract 1.

			node_locs[second_node.parent_node.taxon.label] = [(
			max_depth-second_node.parent_node.root_distance)*brlen_scale,
			highest_y+branch_spacing
			]

			y2 = highest_y+num_internal*branch_spacing

			# Leave space for subtree
			highest_y = highest_y+(num_internal+1)*branch_spacing


			position = [(max_depth-second_node.root_distance)*brlen_scale ,y2]
			if second_node.taxon.label not in node_locs.keys():
				node_locs[second_node.taxon.label] = position

			# figure out second branch location

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length*brlen_scale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			nodes_to_revisit[second_node.taxon.label] = y2-((num_internal-1)*branch_spacing)+branch_spacing # store name of subtree parent and position to start drawing subtree


		else:
			y2 = highest_y+2*branch_spacing

			highest_y = y2

			position = [(max_depth-second_node.root_distance)*brlen_scale ,y2]
			if second_node.taxon.label not in node_locs.keys():
				node_locs[second_node.taxon.label] = position

			# figure out second branch location

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length*brlen_scale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]


			y1 = node_locs[first_node.taxon.label][1]
			position = [(max_depth-second_node.parent_node.root_distance)*brlen_scale ,y2-branch_spacing]
		node = second_node.parent_node

	# Add cartoon arrays to show hypothetical ancestral states
	# plot each array using the coordinates of the array label on the plotted tree.
	rep_indel_report_count = 1 # If repeat_indel found with multiple arrays as possible partners, annotate one on the tree and print the rest to stdout if user wants emphasis.
	rep_indel_message_printed = False # Print a help message the first time a list of possible rep indel partners are found.

	for array, location in node_locs.items():
		# Add label first
		x, y = location
		if "Anc" in array:
			label_color = "#607d8b" #"#cc5500"
		else:
			label_color = "#000000"
		if not no_align_labels:
			ax.text(-0.4, y-0.2, array, ha='right', fontsize=15, color=label_color)
		else:
			ax.text(x-0.4, y-0.2, array, ha='right', fontsize=15, color=label_color)
		# then add branches
		first_node = tree.find_node_with_taxon_label(array)

		# First add branch lengths if user desires

		if branch_lengths:
			if first_node.edge_length != 0:
				ax.text(x+(first_node.edge_length/2)*brlen_scale, y-0.6, first_node.edge_length, ha='center', fontsize=12)
		
		# Draw first branch
		
		x1 = node_locs[first_node.taxon.label][0]
		x2 = node_locs[first_node.taxon.label][0] + first_node.edge_length*brlen_scale
		y1 = node_locs[first_node.taxon.label][1]
		y2 = node_locs[first_node.taxon.label][1]

		ax.plot([x1, x2], [y1, y2], color='black', linewidth = 1, solid_capstyle="butt")

		if len(first_node.sibling_nodes()) == 1:
			second_node = first_node.sibling_nodes()[0]
			# draw second branch

			x1 = node_locs[second_node.taxon.label][0]
			x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length*brlen_scale
			y1 = node_locs[second_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			ax.plot([x1, x2], [y1, y2], color='black', linewidth = 1, solid_capstyle="butt")

			# draw line between branches

			x1 = x2 = node_locs[second_node.taxon.label][0] + second_node.edge_length*brlen_scale
			y1 = node_locs[first_node.taxon.label][1]
			y2 = node_locs[second_node.taxon.label][1]

			ax.plot([x2, x2], [y1, y2], color='black', linewidth = 1, solid_capstyle="butt")

			# Then add spacers and highligh differences

			ancestor = 	array_dict[first_node.parent_node.taxon.label]
			child = array_dict[array]

			if emphasize_diffs:
				if not no_align_cartoons:
					if location[0] != 0:
						ax.plot([-dash_shift, location[0]-dash_shift], [location[1], location[1]], linestyle='--', color='black', linewidth = 1, dashes=(10, 2), alpha=0.5)
					start_pos_x = -label_pad
				else:
					start_pos_x = location[0]-label_pad # Start a bit to the left to leave room for the label
				start_pos_y = location[1]
				spacer_count = 0 # How many spacers have been plotted?
				reshift_loc = 1000
				for n, diff_type in reversed(child.module_lookup.items()):
					spacer = child.aligned[n]

					# Add change info
					if n == reshift_loc:
						start_pos_x-=0.5 # Shift future spacers to make space for line

					if n == diff_type.indices[-1]:
						if diff_type.type != 'shared':
							if diff_type.type == "no_ident":
								# Plot a blue box around the offending arrays
								nspacers = len([child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])
								# First bar
								ax.fill_between([start_pos_x-spacer_size*spacer_count-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2],start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#02a8b1", edgecolor='none')
								# Second bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-spacing/2], start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#02a8b1", edgecolor='none')
								# Top bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y+spacer_width+0.2, start_pos_y+spacer_width+0.4, color="#02a8b1", edgecolor='none')
								# Bottom bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y-spacer_width-0.2, start_pos_y-spacer_width-0.4, color="#02a8b1", edgecolor='none')

								start_pos_x-=0.5 # Shift future spacers a bit to make spacer for this line.
								# Shift again after the indel region
							
							if diff_type.type == "repeated_indel":
								# Plot a red box around repeated indels
								nspacers = len([child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])
								# First bar
								ax.fill_between([start_pos_x-spacer_size*spacer_count-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2],start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#cc3300", edgecolor='none')
								# Second bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-spacing/2], start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#cc3300", edgecolor='none')
								# Top bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y+spacer_width+0.2, start_pos_y+spacer_width+0.4, color="#cc3300", edgecolor='none')
								# Bottom bar
								ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y-spacer_width-0.2, start_pos_y-spacer_width-0.4, color="#cc3300", edgecolor='none')

								# Add Array ID of the array in which the spacers of this predicted repeated_indel can be found
								if len(diff_type.partner) > 2:
									ax.text(start_pos_x-spacer_size*(spacer_count+nspacers/2)-0.5, start_pos_y-spacer_width/2-1.3, "\n".join([diff_type.partner[0], "event {}".format(rep_indel_report_count)]), color="#cc3300", ha='center', fontsize=10)
									if not rep_indel_message_printed:
										print("Repeated indels were identified with multiple possible partners. Those cases will be annotated in the tree png file with the one of the arrays identified as a partner followed by an event number corresponding to one of the lists of partner arrays below:\n\n")
										rep_indel_message_printed = True
									print("Event {}: {}\n\n".format(rep_indel_report_count, " ".join(diff_type.partner)))
									rep_indel_report_count += 1

								else:
									if len(diff_type.partner) == 2:
										ax.text(start_pos_x-spacer_size*(spacer_count+nspacers/2)-0.5, start_pos_y-spacer_width/2-1.3, "\n".join(diff_type.partner), color="#cc3300", ha='center', fontsize=10)
									else:
										ax.text(start_pos_x-spacer_size*(spacer_count+nspacers/2)-0.5, start_pos_y-spacer_width/2-1, diff_type.partner[0], color="#cc3300", ha='center', fontsize=10)

								start_pos_x-=0.5 # Shift future spacers a bit to make spacer for this line.
								# Shift again after the indel region
								reshift_loc = diff_type.indices[0]-1
								
							if diff_type.type == 'indel_gap' or diff_type.type == 'indel_mm' or diff_type.type == 'indel': 
								if spacer == '-':
									ax.plot([start_pos_x-spacer_size*spacer_count-spacing, start_pos_x-spacer_size*spacer_count-spacer_size],[start_pos_y+spacer_width+0.2, start_pos_y-spacer_width-0.2], color="#666666", linewidth=3, solid_capstyle="butt")
									ax.plot([start_pos_x-spacer_size*spacer_count-spacing, start_pos_x-spacer_size*spacer_count-spacer_size],[start_pos_y-spacer_width-0.2, start_pos_y+spacer_width+0.2], color="#666666", linewidth=3, solid_capstyle="butt")
									spacer_count+=1 # Shift future spacers a bit to make spacer for this line.
								else:
									nspacers = len([child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])
									# First bar
									ax.fill_between([start_pos_x-spacer_size*spacer_count-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2],start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#666666", edgecolor='none')
									# Second bar
									ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-spacing/2], start_pos_y+spacer_width+0.4, start_pos_y-spacer_width-0.4, color="#666666", edgecolor='none')
									# Top bar
									ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y+spacer_width+0.2, start_pos_y+spacer_width+0.4, color="#666666", edgecolor='none')
									# Bottom bar
									ax.fill_between([start_pos_x-spacer_size*(spacer_count+nspacers)-0.5-0.5-spacing/2, start_pos_x-spacer_size*spacer_count-spacing/2], start_pos_y-spacer_width-0.2, start_pos_y-spacer_width-0.4, color="#666666", edgecolor='none')

									start_pos_x-=0.5 # Shift future spacers a bit to make spacer for this line.
									# Shift again after the indel region
									reshift_loc = diff_type.indices[0]-1

							elif diff_type.type == "acquisition":
								nspacers = len(diff_type.indices)

								rcParams['path.sketch'] = (25, 60, 1)
								ax.plot(np.linspace(start_pos_x-spacer_size*(spacer_count+nspacers),start_pos_x-spacer_size*spacer_count-spacing,3),[start_pos_y-spacer_width/2-0.4]*3,color="#666666", linewidth=2, solid_capstyle="butt")
								rcParams['path.sketch'] = (0, 0, 0)

							elif diff_type.type == "trailer_loss":
								if spacer == '-': # Draw a single sloped line
									ax.plot([start_pos_x-spacer_size*spacer_count-spacing, start_pos_x-spacer_size*spacer_count-spacer_size],[start_pos_y+spacer_width+0.2, start_pos_y-spacer_width-0.2], color="#666666", linewidth=3, solid_capstyle="butt")
									spacer_count+=1 # Shift future spacers a bit to make spacer for this line.
							elif diff_type.type == "duplication":
								nspacers = len(diff_type.indices)

								ax.plot(np.linspace(start_pos_x-spacer_size*(spacer_count+nspacers),start_pos_x-spacer_size*spacer_count-spacing,3),[start_pos_y-spacer_width/2-0.5]*3,color="#666666", linewidth=3, solid_capstyle="butt")
								


					# Plot spacer cartoon
					if spacer != '-':
						if spacer in spacer_cols_dict.keys():
							spcolour = spacer_cols_dict[spacer]
							sp_width = spacer_width
						else:
							spcolour = ("#000000", "#000000") #black
							sp_width = 0.25*spacer_width
						ax.fill_between([start_pos_x-spacer_size*spacer_count-spacing, start_pos_x-spacer_size*spacer_count-spacing-spacer_size+spacing], start_pos_y-sp_width, start_pos_y+sp_width, color=spcolour[0], edgecolor=spcolour[1], linewidth=outline, joinstyle='miter')
						spacer_count+=1
			if "Anc" in array and fade_ancestral:
				ax.fill_between([start_pos_x-spacer_size*spacer_count-spacing, -label_pad], start_pos_y-sp_width-1, start_pos_y+sp_width+0.6, color="#ffffff", alpha=0.4, joinstyle='miter', zorder=10)


		else: # Draw root ancestral array
			if emphasize_diffs:
			# Then add spacers
				spacers = array_dict[array].spacers
				if not no_align_cartoons:
					if location[0] != 0:
						ax.plot([-dash_shift, location[0]-dash_shift], [location[1], location[1]], linestyle='--', color='black', linewidth = 1, dashes=(10, 2), alpha=0.5)
					start_pos_x = -label_pad
				else:
					start_pos_x = location[0]-label_pad # Start a bit to the left to leave room for the label
				start_pos_y = location[1] 
				for n, spacer in enumerate(reversed(spacers)): # work backwards through the array plotting from right to left
					if spacer in spacer_cols_dict.keys():
						spcolour = spacer_cols_dict[spacer]
						sp_width = spacer_width
					else:
						spcolour = ("#000000", "#000000") #black
						sp_width = 0.25*spacer_width
					
					ax.fill_between([start_pos_x-spacer_size*n-spacing, start_pos_x-spacer_size*n-spacing-spacer_size+spacing], start_pos_y-sp_width, start_pos_y+sp_width, color=spcolour[0], edgecolor=spcolour[1], linewidth=outline, joinstyle='miter')
				if fade_ancestral:
					ax.fill_between([start_pos_x-spacer_size*spacer_count-spacing, -label_pad], start_pos_y-sp_width-1, start_pos_y+sp_width+0.6, color="#ffffff", alpha=0.4, joinstyle='miter', zorder=10)


		if not emphasize_diffs:
			# Then add spacers
			spacers = array_dict[array].spacers
			if not no_align_cartoons:
				if location[0] != 0:
					ax.plot([-dash_shift, location[0]-dash_shift], [location[1], location[1]], linestyle='--', color='black', linewidth = 1, dashes=(10, 2), alpha=0.5)
				start_pos_x = -label_pad
			else:
				start_pos_x = location[0]-label_pad # Start a bit to the left to leave room for the label
			start_pos_y = location[1] 
			for n, spacer in enumerate(reversed(spacers)): # work backwards through the array plotting from right to left
				if spacer in spacer_cols_dict.keys():
					spcolour = spacer_cols_dict[spacer]
					sp_width = spacer_width
				else:
					spcolour = ("#000000", "#000000") #black
					sp_width = 0.25*spacer_width
				
				ax.fill_between([start_pos_x-spacer_size*n-spacing, start_pos_x-spacer_size*n-spacing-spacer_size+spacing], start_pos_y-sp_width, start_pos_y+sp_width, color=spcolour[0], edgecolor=spcolour[1], linewidth=outline, joinstyle='miter')


	ymin, ymax = ax.get_ylim()
	xmin, xmax = ax.get_xlim()
	
	fig.set_size_inches(0.15*(xmax-xmin), 0.4*(ymax-ymin))

	plt.axis('off')
	plt.tight_layout()
	plt.savefig(filename, dpi=dpi)
