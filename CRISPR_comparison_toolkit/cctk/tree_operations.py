from string import ascii_lowercase
from math import ceil, log, sin, pi
from itertools import product
from copy import copy
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np



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

	node_list = [n for n in yield_nodes(tree.seed_node)]

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


def draw_branches(tree, node_locs, ax, branch_lengths=True, brlen_scale=0.5,
	font_scale=1, line_scale=1, label_text_size=False, annot_text_size=False):

	for name, location in node_locs.items():
		x, y = location
		node = tree.find_node_with_taxon_label(name)
		
		# Draw branch from this node to its parent
		
		# identify x location for parent depth if the parent is not seed
		if node != tree.seed_node and node.parent_node.taxon != None:
			x2 = node_locs[node.parent_node.taxon.label][0]
		else:
			x2 = x + brlen_scale

		ax.plot([x, x2], [y, y], color='black', linewidth=line_scale,
			solid_capstyle="round")

		# If internal node, draw line connecting lines to children
		if node.is_internal():
			child_locations = [
				node_locs[i.taxon.label] for i in node.child_nodes()]
			# Draw line between the highest and smallest child y value
			y1 = max([i[1] for i in child_locations])
			y2 = min([i[1] for i in child_locations])

			ax.plot([x, x], [y1, y2], color='black', linewidth=line_scale,
				solid_capstyle="round")

	return ax


def add_labels(tree, node_locs, ax, branch_lengths=True,font_scale=1,
	line_scale=1, label_text_size=False, annot_text_size=False,
	no_align_labels=False, brlen_scale=0.5,):

	if annot_text_size:
		brlen_font_size = annot_text_size
	else:
		brlen_font_size = 8*font_scale
	
	for name, location in node_locs.items():
		# Extract relevant node from tree
		node = tree.find_node_with_taxon_label(name)
		
		# Add node label
		x, y = location
		label_color = "#000000"
		if not no_align_labels:
			# Draw labels at the most extreme x position (negativen number)
			label_x = min([int(i[0]) for i in node_locs.values()])
			# Draw dashed line from node to label for visual aid
			ax.plot([label_x, x], [y,y], linestyle='--', color='black',
				linewidth=line_scale*0.5, dashes=(10, 2), alpha=0.3)
		else:
			label_x = x

		ax.text(label_x-0.1, y, name, ha='right', va='center_baseline', 
			fontsize=label_text_size, color=label_color)

		# add branch lengths if user desires and if branch has a length
		if not branch_lengths:
			continue	
		if node.edge_length == 0:
			continue
		
		ax.text(x + (node.edge_length/2)*brlen_scale, y-0.6,node.edge_length,
			ha='center', va='top', fontsize=brlen_font_size)
	
	return ax


def add_cartoons(node_locs, ax, array_dict, spacer_cols_dict, spacer_width=0.3, 
	spacer_outline=3, spacer_spacing=0.45, spacer_size=2.5, label_pad=5,
	fade_ancestral=False, no_align_cartoons=False, emphasize_diffs=False,
	v_scaling_factor=1, annot_text_size=False):

	# Add cartoon arrays to show hypothetical ancestral states
	# plot each array using the coordinates of the array label on the plotted
	# tree.
	
	# If repeat_indel found with multiple arrays as possible partners,
	# annotate one on the tree and print the rest to stdout if user wants
	# emphasis.
	rep_indel_report_count = 1 
	
	# Print a help message the first time a list of possible rep indel partners
	# are found.
	rep_indel_message_printed = False 
	
	for array, location in node_locs.items():
		child = array_dict[array]
		x, y = location
		if not no_align_cartoons:
			x = min([int(i[0]) for i in node_locs.values()])

		start_pos_x = x-label_pad

		spacer_count = 0 # How many spacers have been plotted?
		reshift_loc = 1000 # init to number that won't be seen

		for n, diff_type in reversed(child.module_lookup.items()):
			spacer = child.aligned[n]

			nspacers = len([
				child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])

			# Add change info
			if n == reshift_loc:
				start_pos_x-=0.5 # Shift future spacers to make space for line
	
			if n == diff_type.indices[-1]:
				# If this is the start of a new module (working from the
				# back of the module) then process this module
				if diff_type == "shared":
					# If these spacers are shared then no annotation.
					continue
			
				elif diff_type.type == "no_ident":
					print('noident')
					# Plot blue box around spacers
					ax = plot_box(
						ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, "#02a8b1",
						box_thickness=0.3, hpad=1, vpad=1,
						v_scaling_factor=v_scaling_factor)

					# Shift future spacers to make spacer for this line.
					start_pos_x-=0.5
					# Shift again after the indel region for next line
					reshift_loc = diff_type.indices[0]-1
					
				elif diff_type.type == "repeated_indel":
					print('rep')
					# Plot a red box around repeated indels
					ax = plot_box(
						ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, "#cc3300",
						box_thickness=0.5, hpad=1, vpad=0.2,
							v_scaling_factor=v_scaling_factor)

					# Shift future spacers to make spacer for this line.
					start_pos_x-=0.5
					# Shift again after the indel region for next line
					reshift_loc = diff_type.indices[0]-1

					# Add Array ID of the array in which the spacers of this predicted repeated_indel can be found
					num_partners = len(diff_type.partner)
					
					ax = plot_rep_indel_text(
						ax, nspacers, start_pos_x, y, spacer_size, spacer_count,
						spacer_spacing,	spacer_width, num_partners,
						rep_indel_report_count,	rep_indel_message_printed,
						annot_text_size=annot_text_size)

					rep_indel_report_count += 1
					rep_indel_message_printed = True
					
				elif diff_type.type in ['indel_gap', 'indel_mm', 'indel']:
					if spacer == '-':
						# Draw an X as this is a deletion
						ax = plot_x_or_slash(ax, start_pos_x, y,
							spacer_size, spacer_count, spacer_spacing, spacer_width,
							'x')
						# Shift future spacers a bit to make spacer for
						# this line.
						spacer_count+=1
					else:
						# Insertion so plot a grey box
						# hpad = box_thickness + start_pos_x shift
						ax = plot_box(
							ax, nspacers-1, start_pos_x, y, spacer_size,
							spacer_count, spacer_spacing, spacer_width,
							"#02a8b1", box_thickness=0.5, hpad=1, vpad=0.2,
							v_scaling_factor=v_scaling_factor)

						# Shift future spacers to make spacer for this line.
						start_pos_x-=0.5
						# Shift again after the indel region for next line
						reshift_loc = diff_type.indices[0]-1

				elif diff_type.type == "acquisition":
					ax = plot_line(ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, wavy=True)

				elif diff_type.type == "trailer_loss":
					if spacer == '-': # Draw a single sloped line
						ax = plot_x_or_slash(ax, start_pos_x, y, 
							spacer_size, spacer_count, spacer_spacing, spacer_width,
							'/')
						 # Shift future spacers a bit to make spacer for this
						 # line.
						spacer_count+=1
				elif diff_type.type == "duplication":
					ax = plot_line(ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, wavy=False)

			# Plot spacer cartoon
			if spacer != '-':
				if spacer in spacer_cols_dict.keys():
					spacer_colour = spacer_cols_dict[spacer]
					sp_width = spacer_width
				else:
					spacer_colour = ("#000000", "#000000") #black
					sp_width = 0.25*spacer_width
				ax = plot_spacer(ax, start_pos_x, y, spacer_size, sp_width,
					spacer_colour, spacer_outline, spacer_count, spacer_spacing,
					v_scaling_factor)
				spacer_count+=1

	return ax


def plot_spacer(ax, x, y, spacer_size, spacer_width, spacer_colour,
	spacer_outline, spacer_count, spacer_spacing, v_scaling_factor):

	ax.fill_between([x-spacer_size*spacer_count-spacer_spacing,
		x-spacer_size*spacer_count-spacer_size],
		y-spacer_width/2, y+spacer_width/2,	color=spacer_colour[0],
		edgecolor='none',
		joinstyle='miter') # linewidth=spacer_outline, edgecolor=spacer_colour[1],
	
	# Plot outline
	ax = plot_box(ax, 0, x, y, spacer_size, spacer_count, spacer_spacing,
		spacer_width, spacer_colour[1], box_thickness=spacer_outline, hpad=0, vpad=0,
		v_scaling_factor=v_scaling_factor)

	return ax


def plot_box(ax, nspacers, x, y, spacer_size, spacer_count, spacer_spacing,
	spacer_width, box_colour, box_thickness, hpad, vpad, v_scaling_factor):
	
	box_thickness_v = box_thickness * v_scaling_factor


	# Right bar
	ax.fill_between(
		[x-spacer_size*spacer_count-spacer_spacing-box_thickness/2,
			x-spacer_size*spacer_count-spacer_spacing+box_thickness/2],
		y+spacer_width/2+vpad+box_thickness_v/2, y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none')
	# Left bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size-box_thickness/2-hpad,
			x-spacer_size*(spacer_count+nspacers)-spacer_size+box_thickness/2-hpad],
		y+spacer_width/2+vpad+box_thickness_v/2, y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none')

	# Top bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size-box_thickness/2-hpad,
			x-spacer_size*spacer_count-spacer_spacing+box_thickness/2],
		y+spacer_width/2+vpad+box_thickness_v/2, 
		y+spacer_width/2+vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none')

	# Bottom bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size-box_thickness/2-hpad,
			x-spacer_size*spacer_count-spacer_spacing+box_thickness/2],
		y-spacer_width/2-vpad+box_thickness_v/2, 
		y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none')

	return ax


def plot_x_or_slash(ax, x, y, spacer_size, spacer_count,
	spacer_spacing, spacer_width, x_or_slash):
		
	shift = 0.3*(spacer_width/spacer_size)

	xs = np.linspace(x-spacer_size*spacer_count-spacer_spacing,
		x-spacer_size*spacer_count-spacer_size, 100)

	ys = np.linspace(y-spacer_width/2-0.1,
		y+spacer_width/2+0.2, 100)
	ys2 = [i-shift for i in ys]

	ax.fill_between(xs, ys, ys2, color="#666666", joinstyle='round')
	
	if x_or_slash == 'x':
		# Reverse order of y values for second line
		ys = [i for i in reversed(ys)]
		ys2.reverse()
		ax.fill_between(xs, ys, ys2, color="#666666")

	return ax


def plot_line(ax, nspacers, x, y, spacer_size, spacer_count, spacer_spacing, 
	spacer_width, wavy=False):
		
	padded_y = y-spacer_width/2-0.4
	shift = 0.1

	xs = np.linspace(x-spacer_size*(spacer_count+nspacers),
		x-spacer_size*spacer_count-spacer_spacing, nspacers*100)

	if wavy:
		ys = [padded_y+0.1*sin(2*pi*i) for i in xs]
		ys2 = [(padded_y-shift)+0.1*sin(2*pi*i) for i in xs]
		

	else:
		ys = [padded_y for _ in xs]
		ys2 = [(padded_y-shift) for _ in xs]

	ax.fill_between(xs, ys, ys2, color="#666666")

	return ax


def plot_rep_indel_text(ax, nspacers, x, y, spacer_size, spacer_count,
	spacer_width, partners, rep_indel_report_count, rep_indel_message_printed,
	font_scale=1, annot_text_size=False):
	
	if annot_text_size:
		font_size = annot_text_size
	else:
		font_size = 1.5*font_scale

	if len(partners) > 2:
		vpad = 1.3
		message = "\n".join(
				[partners[0],
				"event {}".format(rep_indel_report_count)])
		
		if not rep_indel_message_printed:
			print("Repeated indels were identified with multiple possible "
				+ "partners. Those cases will be annotated in the tree png "
				+ "file with the one of the arrays identified as a partner "
				+ "followed by an event number corresponding to one of the "
				+ "lists of partner arrays below:\n\n")
		print("Event {}: {}\n\n".format(
			rep_indel_report_count,
			" ".join(diff_type.partner)))
	else:
		if len(partners) == 2:
			vpad = 1.3
			message = "\n".join(partners)
		else:
			vpad = 1
			message = partners[0]

		ax.text(x-spacer_size*(spacer_count+nspacers/2)-0.5,
			y-spacer_width/2-vpad, message,
			color="#cc3300", ha='center', fontsize=font_size)

	return ax
			

def calc_vh_ratio_and_label_pad(tree, array_dict, spacing, spacer_size, 
	branch_spacing, brlen_scale, fig_w, fig_h, text_size):
	
	# Find approx max horizontal axis range
	max_h = 0 # Initialize at 0
	for leaf in tree.leaf_node_iter():
		name = leaf.taxon.label
		brlen = leaf.root_distance
		nspacers = len(array_dict[name].spacers)
		total_h = (brlen*brlen_scale + nspacers*(spacing+spacer_size))
		if total_h > max_h:
			max_h = total_h
	
	max_v = len(tree.nodes())*branch_spacing

	# Find label padding size
	longest_label = max([i for i in array_dict.keys()], key=len)

	# Find approx asymptote point as label pad will increase with axis
	# limits, but increasing label pad will increase axis limits.
	last_max = copy(max_h)-1 # Init at higher value to init loop
	new_max_h = copy(max_h) # Store max_h updated with each pad
	
	while new_max_h > last_max+0.0001: # Stop once label pad stable
		last_max = copy(new_max_h)
		label_pad = calc_label_pad_size(
			longest_label, text_size, fig_w, fig_h,
			xlim=(0, new_max_h), ylim=(0, max_v))
		new_max_h = max_h+label_pad

	v_scaling_factor = max_v/new_max_h

	return v_scaling_factor, label_pad


def calc_label_pad_size(label, text_size, fig_w, fig_h, xlim, ylim):

	f,ax = plt.subplots(figsize=(fig_w,fig_h))
	plt.ylim(ylim)
	plt.xlim(xlim)
	r = f.canvas.get_renderer()
	t = plt.text(5, 5, label, fontsize=text_size)

	bb = t.get_window_extent(renderer=r)
	width = Bbox(ax.transData.inverted().transform(bb)).width

	# Add a bit of space so the cartoons aren't right up against labels
	width += 0.5

	return width


def plot_tree(tree, array_dict, filename, spacer_cols_dict, 
	branch_lengths=False, emphasize_diffs=False, dpi=600, line_scale=1,
	brlen_scale=0.5, branch_spacing=2.3, font_scale=1, fig_h=1, fig_w=1,
	no_align_cartoons=False, no_align_labels=False, fade_ancestral=False,
	label_text_size=False, annot_text_size=False):

	outline = 0.3 #line_widths = 3 # Thickness of lines
	spacing = 0.5
	spacer_size = 2

	if not label_text_size:
		label_text_size = 10 * font_scale

	node_locs = find_node_locs(tree, branch_spacing=branch_spacing,
		brlen_scale=brlen_scale)

	v_scaling_factor, label_pad = calc_vh_ratio_and_label_pad(tree, array_dict,
		spacing, spacer_size, branch_spacing, brlen_scale, fig_w, fig_h,
		label_text_size)

	# Set spacer width based on vertical vs horizontal ratio. Max 1.3
	spacer_width = min([spacer_size*v_scaling_factor, 1.3])

	fig, ax = plt.subplots(figsize=(fig_w,fig_h))

	ax = draw_branches(tree, node_locs, ax, font_scale=font_scale,
		line_scale=line_scale, branch_lengths=branch_lengths,
		label_text_size=label_text_size, annot_text_size=annot_text_size,
		brlen_scale=brlen_scale)

	ax = add_labels(tree, node_locs, ax, font_scale=font_scale,
		branch_lengths=branch_lengths, label_text_size=label_text_size,
		annot_text_size=annot_text_size, no_align_labels=no_align_labels,
		brlen_scale=brlen_scale)

	ax = add_cartoons(node_locs, ax, array_dict, spacer_cols_dict,
		label_pad=label_pad, spacer_size=spacer_size, spacer_outline=outline,
		spacer_spacing=spacing, spacer_width=spacer_width,
		no_align_cartoons=no_align_cartoons, emphasize_diffs=emphasize_diffs,
		v_scaling_factor=v_scaling_factor, annot_text_size=annot_text_size)

	plt.axis('off')
	plt.tight_layout()

	plt.savefig(filename, dpi=dpi, bbox_inches='tight')

