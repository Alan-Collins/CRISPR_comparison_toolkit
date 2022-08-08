from math import ceil, log, sin, pi
from copy import copy
import itertools
import sys

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np

from . import (
	tree_operations,
	sequence_operations)

def draw_branches(tree, node_locs, ax, branch_lengths=True, brlen_scale=1,
	font_scale=1, line_scale=1, label_text_size=False, annot_text_size=False,
	branch_support=False):

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

			# add branch support if appropriate
			if branch_support and node != tree.seed_node:
				support = node.node_support
				ax.scatter(
					x, y,
					c=support,
					cmap='cividis_r',
					vmin=0, vmax=100)


def add_labels(tree, node_locs, ax, branch_lengths=True,font_scale=1,
	line_scale=1, label_text_size=False, annot_text_size=False,
	no_align_labels=False, brlen_scale=1, no_fade_ancestral=False):

	if annot_text_size:
		brlen_font_size = annot_text_size
	else:
		brlen_font_size = 6*font_scale

	leaves = [i.taxon.label for i in tree.leaf_nodes()]
	
	for name, location in node_locs.items():
		# Extract relevant node from tree
		node = tree.find_node_with_taxon_label(name)
		
		# Add node label
		x, y = location
		label_color = "#000000"

		# Fade if not leaf node if no_fade_ancestral == False
		if not no_fade_ancestral and node.taxon.label not in leaves:
			label_color += "80"
		if not no_align_labels:
			# Draw labels at the most extreme x position (negativen number)
			label_x = min([float(i[0]) for i in node_locs.values()])
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


def add_cartoons(tree, node_locs, ax, array_dict, non_singleton_spacers, spacer_cols_dict,
	spacer_width=0.3, spacer_outline=0.3, spacer_spacing=0.45, spacer_size=2.5,
	label_pad=5, no_fade_ancestral=False, no_align_cartoons=False,
	emphasize_diffs=False, v_scaling_factor=1, annot_text_size=False,
	label_pad_dict={}):

	# Add cartoon arrays to show hypothetical ancestral states
	# plot each array using the coordinates of the array label on the plotted
	# tree.
	
	# List of leaf nodes for identification of internal nodes used if
	# no_fade_ancestral == False

	leaves = [i.taxon.label for i in tree.leaf_nodes()]

	# If repeat_indel found with multiple arrays as possible partners,
	# annotate one on the tree and print the rest to stdout if user wants
	# emphasis.
	rep_indel_report_count = 1 
	
	# Print a help message the first time a list of possible rep indel partners
	# are found.
	rep_indel_message_printed = False 
	
	for array, location in node_locs.items():
		if array not in leaves:
			ancestral = True
		else:
			ancestral = False

		child = array_dict[array]
		x, y = location
		if not no_align_cartoons:
			x = min([int(i[0]) for i in node_locs.values()])
			start_pos_x = x-label_pad

		else:
			start_pos_x = x - label_pad_dict[array]

		spacer_count = 0 # How many spacers have been plotted?
		reshift_loc = 1000 # init to number that won't be seen

		for n, diff_type in reversed(child.module_lookup.items()):
			spacer = child.aligned[n]

			nspacers = len([
				child.aligned[i] for i in diff_type.indices if child.aligned[i] != '-'])

			# Add change info
			if n == reshift_loc:
				start_pos_x-=0.5 # Shift future spacers to make space for line
	
			if n == diff_type.indices[-1] and emphasize_diffs:
				# If this is the start of a new module (working from the
				# back of the module) then process this module
				if diff_type == "shared":
					# If these spacers are shared then no annotation.
					continue
			
				elif diff_type.type == "no_ident":
					colour = "#02a8b1"
					if not no_fade_ancestral and ancestral:
						colour += "80"
					# Plot blue box around spacers
					plot_box(
						ax, nspacers-1, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, colour,
						box_thickness=0.5, hpad=1, vpad=0.2,
						v_scaling_factor=v_scaling_factor)

					# Shift future spacers to make spacer for this line.
					start_pos_x-=0.5
					# Shift again after the indel region for next line
					reshift_loc = diff_type.indices[0]-1
					
				elif diff_type.type == "repeated_indel":
					colour = "#cc3300"
					if not no_fade_ancestral and ancestral:
						colour += "80"
					# Plot a red box around repeated indels
					plot_box(
						ax, nspacers-1, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, colour,
						box_thickness=0.5, hpad=1, vpad=0.2,
						v_scaling_factor=v_scaling_factor)

					# Shift future spacers to make spacer for this line.
					start_pos_x-=0.5
					# Shift again after the indel region for next line
					reshift_loc = diff_type.indices[0]-1

					# Add Array ID of the array in which the spacers of this predicted repeated_indel can be found
					
					plot_rep_indel_text(
						ax, nspacers, start_pos_x, y, spacer_size, spacer_count,
						spacer_width, diff_type.partner,
						rep_indel_report_count,	rep_indel_message_printed,
						annot_text_size=annot_text_size)
					if len(diff_type.partner) > 1:
						rep_indel_report_count += 1
					rep_indel_message_printed = True
					
				elif diff_type.type in ['indel_gap', 'indel_mm', 'indel']:
					if spacer == '-':
						colour = "#666666"
						if not no_fade_ancestral and ancestral:
							colour += "80"
						# Draw an X as this is a deletion
						plot_x_or_slash(ax, start_pos_x, y,
							spacer_size, spacer_count, spacer_spacing, 
							spacer_width, 'x', colour)
						# Shift future spacers a bit to make spacer for
						# this line.
						spacer_count+=1
					else:
						colour = "#666666"
						if not no_fade_ancestral and ancestral:
							colour += "80"
						# Insertion so plot a grey box
						# hpad = box_thickness + start_pos_x shift
						plot_box(
							ax, nspacers-1, start_pos_x, y, spacer_size,
							spacer_count, spacer_spacing, spacer_width,
							colour, box_thickness=0.5, hpad=1, vpad=0.2,
							v_scaling_factor=v_scaling_factor)

						# Shift future spacers to make spacer for this line.
						start_pos_x-=0.5
						# Shift again after the indel region for next line
						reshift_loc = diff_type.indices[0]-1

				elif diff_type.type == "acquisition":
					colour = "#666666"
					if not no_fade_ancestral and ancestral:
						colour += "80"
					plot_line(ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, wavy=True,
						colour=colour)

				elif diff_type.type == "trailer_loss":
					if spacer == '-': # Draw a single sloped line
						colour = "#666666"
						if not no_fade_ancestral and ancestral:
							colour += "80"
						plot_x_or_slash(ax, start_pos_x, y, 
							spacer_size, spacer_count, spacer_spacing, spacer_width,
							'/', colour)
						 # Shift future spacers a bit to make spacer for this
						 # line.
						spacer_count+=1
				elif diff_type.type == "duplication":
					colour = "#666666"
					if not no_fade_ancestral and ancestral:
						colour += "80"
					plot_line(ax, nspacers, start_pos_x, y, spacer_size,
						spacer_count, spacer_spacing, spacer_width, wavy=False,
						colour=colour)

			# Plot spacer cartoon
			if spacer == '-':
				continue
			if spacer in non_singleton_spacers:
				spacer_colour = spacer_cols_dict[spacer]
				sp_width = spacer_width
			else:
				if spacer in spacer_cols_dict:
					spacer_colour = spacer_cols_dict[spacer]
					spacer_outline = 0.2
				else:
					spacer_colour = ("#000000", "#000000") #black
				sp_width = 0.25*spacer_width
			plot_spacer(ax, start_pos_x, y, spacer_size, sp_width,
				spacer_colour, spacer_outline, spacer_count, spacer_spacing,
				v_scaling_factor)
			spacer_outline = 0.3
			if not no_fade_ancestral and ancestral:
				# Plot translucent box over spacers to fade them
				ax.fill_between(
					[(start_pos_x-spacer_size*spacer_count-spacer_spacing
						+spacer_outline/2),
					(start_pos_x-spacer_size*spacer_count-spacer_size
						-spacer_outline/2)],
					y-spacer_width/2-(spacer_outline*v_scaling_factor)/2,
					y+spacer_width/2+(spacer_outline*v_scaling_factor)/2,
					edgecolor='none', linewidth=0, color="#ffffff80")
			spacer_count+=1


def plot_spacer(ax, x, y, spacer_size, spacer_width, spacer_colour,
	spacer_outline, spacer_count, spacer_spacing, v_scaling_factor):

	ax.fill_between([x-spacer_size*spacer_count-spacer_spacing,
		x-spacer_size*spacer_count-spacer_size],
		y-spacer_width/2, y+spacer_width/2,	color=spacer_colour[0],
		edgecolor='none', linewidth=0)
	
	# Plot outline
	plot_box(ax, 0, x, y, spacer_size, spacer_count, spacer_spacing,
		spacer_width, spacer_colour[1], box_thickness=spacer_outline, hpad=0, vpad=0,
		v_scaling_factor=v_scaling_factor)


def plot_box(ax, nspacers, x, y, spacer_size, spacer_count, spacer_spacing,
	spacer_width, box_colour, box_thickness, hpad, vpad, v_scaling_factor):
	
	box_thickness_v = box_thickness * v_scaling_factor


	# Right bar
	ax.fill_between(
		[x-spacer_size*spacer_count-spacer_spacing-box_thickness/2,
			x-spacer_size*spacer_count-spacer_spacing+box_thickness/2],
		y+spacer_width/2+vpad+box_thickness_v/2, y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none', linewidth=0)
	# Left bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size-box_thickness/2-hpad,
			x-spacer_size*(spacer_count+nspacers)-spacer_size+box_thickness/2-hpad],
		y+spacer_width/2+vpad+box_thickness_v/2, y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none', linewidth=0)

	# Top bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size+box_thickness/2-hpad,
			x-spacer_size*spacer_count-spacer_spacing-box_thickness/2],
		y+spacer_width/2+vpad+box_thickness_v/2, 
		y+spacer_width/2+vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none', linewidth=0)

	# Bottom bar
	ax.fill_between(
		[x-spacer_size*(spacer_count+nspacers)-spacer_size+box_thickness/2-hpad,
			x-spacer_size*spacer_count-spacer_spacing-box_thickness/2],
		y-spacer_width/2-vpad+box_thickness_v/2, 
		y-spacer_width/2-vpad-box_thickness_v/2,
		color=box_colour, edgecolor='none', linewidth=0)


def plot_x_or_slash(ax, x, y, spacer_size, spacer_count,
	spacer_spacing, spacer_width, x_or_slash, colour):
		
	shift = 0.1#*(spacer_width/spacer_size)

	xs = np.linspace(x-spacer_size*spacer_count-spacer_spacing,
		x-spacer_size*spacer_count-spacer_size, 100)

	ys = np.linspace(y+shift/2-spacer_width/2,
		y+shift/2+spacer_width/2, 100)
	ys2 = [i-shift for i in ys]

	ax.fill_between(xs, ys, ys2, color=colour, edgecolor='none', linewidth=0)
	
	if x_or_slash == 'x':
		# Reverse order of y values for second line
		ys = [i for i in reversed(ys)]
		ys2.reverse()
		ax.fill_between(xs, ys, ys2, color=colour, edgecolor='none',
			linewidth=0)


def plot_line(ax, nspacers, x, y, spacer_size, spacer_count, spacer_spacing, 
	spacer_width, wavy=False, colour="#666666"):
		
	padded_y = y-spacer_width/2-0.2
	shift = 0.1

	xs = np.linspace(x-spacer_size*(spacer_count+nspacers),
		x-spacer_size*spacer_count-spacer_spacing, nspacers*100)

	if wavy:
		ys = [padded_y+0.05*sin(2*pi*i) for i in xs]
		ys2 = [(padded_y-shift)+0.05*sin(2*pi*i) for i in xs]
		

	else:
		ys = [padded_y for _ in xs]
		ys2 = [(padded_y-shift) for _ in xs]

	ax.fill_between(xs, ys, ys2, color=colour, edgecolor='none', linewidth=0)


def plot_rep_indel_text(ax, nspacers, x, y, spacer_size, spacer_count,
	spacer_width, partners, rep_indel_report_count, rep_indel_message_printed,
	font_scale=1, annot_text_size=False):
	
	if annot_text_size:
		font_size = annot_text_size
	else:
		font_size = 6*font_scale

	vpad = 0.8

	if len(partners) > 1:
		message = "event {}".format(rep_indel_report_count)
		
		if not rep_indel_message_printed:
			sys.stderr.write("Repeated indels were identified with multiple possible "
				"partners. Those cases will be annotated in the tree png "
				"file with the one of the arrays identified as a partner "
				"followed by an event number corresponding to one of the "
				"lists of partner arrays below:\n\n")
		sys.stderr.write("Event {}: {}\n\n".format(
			rep_indel_report_count,
			" ".join(partners)))
	else:
		message = partners[0]

	ax.text(x-spacer_size*(spacer_count+nspacers/2)-0.5,
		y-spacer_width/2-vpad, message,
		color="#cc3300", ha='center', fontsize=font_size)


def calc_vh_ratio_and_label_pad_tree(tree, array_dict, spacing, spacer_size, 
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


def calc_font_scale(array_dict, fig_w, fig_h, text_size, axis_height,
	max_text_height):
	# find longest label

	longest_label = ""
	longest_len = 0

	for label in array_dict.keys():
		f,ax = plt.subplots(figsize=(fig_w,fig_h))
		plt.ylim(100)
		plt.xlim(100)
		r = f.canvas.get_renderer()
		t = plt.text(1, 1, label, fontsize=text_size, ha='right',
			va='bottom')

		bb = t.get_window_extent(renderer=r)
		width = abs(Bbox(ax.transData.inverted().transform(bb)).width)
		plt.close()

		if width > longest_len:
			longest_len = width
			longest_label = label
	scale = 1

	# if needed, shrink text until smaller than 1/4 of plot axes

	width = longest_len
	while width > 25:
		scale*=0.9
		scaled_text = text_size*scale

		f,ax = plt.subplots(figsize=(fig_w,fig_h))
		plt.ylim(100)
		plt.xlim(100)
		r = f.canvas.get_renderer()
		t = plt.text(1, 1, longest_label, fontsize=scaled_text, ha='right',
			va='bottom')

		bb = t.get_window_extent(renderer=r)
		width = abs(Bbox(ax.transData.inverted().transform(bb)).width)
		plt.close()

	# Now check height. Shrink font if too tall

	f,ax = plt.subplots(figsize=(fig_w,fig_h))
	plt.ylim(axis_height)
	plt.xlim(10)
	r = f.canvas.get_renderer()
	t = plt.text(1, 1, label, fontsize=text_size, ha='right',
		va='bottom')

	bb = t.get_window_extent(renderer=r)
	height = abs(Bbox(ax.transData.inverted().transform(bb)).height)
	plt.close()

	height_scale = max_text_height/height

	if height_scale < scale:
		scale = height_scale

	return scale	


def calc_vh_ratio_and_label_pad_diffplot(array_dict, spacer_size,
	spacer_spacing, ylim, text_size=10, fig_w=3, fig_h=3):

	max_h = 0
	for array, spacers in array_dict.items():
		spacer_h = len(spacers)*(spacer_size+spacer_spacing)
		last_max = spacer_h-1 # Init at higher value to init loop
		new_max_h = spacer_h # Store max_h updated with each pad
		
		while new_max_h > last_max+0.0001: # Stop once label pad stable
			last_max = copy(new_max_h)
			pad = calc_label_pad_size(array, text_size, fig_w, fig_h,
				xlim=(0,new_max_h), ylim=ylim)
			new_max_h = spacer_h-pad
		if new_max_h > max_h:
			max_h = copy(new_max_h)
	v_scaling_factor = (ylim[1]-ylim[0])/max_h

	return v_scaling_factor


def calc_label_pad_size(label, text_size, fig_w, fig_h, xlim, ylim):

	f,ax = plt.subplots(figsize=(fig_w,fig_h))
	plt.ylim(ylim)
	plt.xlim(xlim)
	r = f.canvas.get_renderer()
	t = plt.text(1, 1, label, fontsize=text_size, ha='right',
		va='center_baseline')

	bb = t.get_window_extent(renderer=r)
	width = Bbox(ax.transData.inverted().transform(bb)).width
	plt.close()

	# Add a bit of space so the cartoons aren't right up against labels
	width += 0.5

	return width


def draw_diffplot_lines(ax, y, array_order, array_dict, spacer_colours,
				spacer_size, spacer_width, spacer_spacing, spacer_outline,
				vh_ratio, line_width, connection_outline):
	array1_id = array_order[y]
	array2_id = array_order[y+1]
	array1_spacers_rev = [i for i in reversed(array_dict[array1_id])]
	array2_spacers_rev = [i for i in reversed(array_dict[array2_id])]
	for n, spacer in enumerate(array1_spacers_rev):
		if spacer not in spacer_colours:
			continue
		spcolour = spacer_colours[spacer]
		array1_indices = sequence_operations.find_indices(
			array1_spacers_rev, spacer)
		array2_indices = sequence_operations.find_indices(
			array2_spacers_rev, spacer)
		for a, b in itertools.product(array1_indices, array2_indices):
			sp_x1 = -a*spacer_size-(spacer_size/2)-(spacer_spacing/2)
			sp_y1 = y+1+(spacer_width/2)+(spacer_outline*vh_ratio)/4
			sp_x2 = -b*spacer_size-(spacer_size/2)-(spacer_spacing/2)
			sp_y2 = y+2-(spacer_width/2)-(spacer_outline*vh_ratio)/4
			
			if connection_outline:
				ax.plot(
					[sp_x1, sp_x2],
					[sp_y1, sp_y2],
					color=spcolour[1],
					linewidth=2*line_width,
					solid_capstyle="round",
					label=spacer, 
					zorder=0
					)

			ax.plot(
				[sp_x1, sp_x2],
				[sp_y1, sp_y2],
				color=spcolour[0],
				linewidth=1*line_width,
				solid_capstyle="round",
				label=spacer,
				zorder=0.5
				)


def plot_tree(tree, array_dict, filename, non_singleton_spacers, spacer_cols_dict, 
	branch_lengths=False, emphasize_diffs=False, dpi=600, line_scale=1,
	brlen_scale=1, branch_spacing=2, font_scale=1, fig_h=1, fig_w=1,
	no_align_cartoons=False, no_align_labels=False, no_fade_ancestral=False,
	label_text_size=False, annot_text_size=False, branch_support=False):

	outline = 0.3 # Thickness of lines
	spacing = 0.5
	spacer_size = 2

	# scale branches to take up 25% of x axis
	brlen_scale *= 25/tree.length()
	
	tree_height = branch_spacing*len(tree.nodes())

	if not label_text_size:
		font_scale = calc_font_scale(array_dict, fig_w, fig_h, 10, tree_height, 1.8)
		label_text_size = 10 * font_scale

	if not annot_text_size:
		annot_text_size = 6 * font_scale

	node_locs = tree_operations.find_node_locs(tree, branch_spacing=branch_spacing,
		brlen_scale=brlen_scale)

	v_scaling_factor, label_pad = calc_vh_ratio_and_label_pad_tree(
		tree, array_dict, spacing, spacer_size, branch_spacing, brlen_scale,
		fig_w, fig_h, label_text_size)

	# Calculate label size if not aligning cartoons
	label_pad_dict = {}
	if no_align_cartoons:
		for array in array_dict.values():
			max_v = len(tree.nodes())*branch_spacing
			max_h = max_v/v_scaling_factor
			label_pad_dict[array.id] = calc_label_pad_size(
			array.id, label_text_size, fig_w, fig_h,
			xlim=(0,max_h), ylim=(0,max_v))

	# Set spacer width based on vertical vs horizontal ratio. Max 1.3
	spacer_width = min([spacer_size*v_scaling_factor, 1.3])

	fig, ax = plt.subplots(figsize=(fig_w,fig_h))

	draw_branches(tree, node_locs, ax, font_scale=font_scale,
		line_scale=line_scale, branch_lengths=branch_lengths,
		label_text_size=label_text_size, annot_text_size=annot_text_size,
		brlen_scale=brlen_scale, branch_support=branch_support)

	add_labels(tree, node_locs, ax, font_scale=font_scale,
		branch_lengths=branch_lengths, label_text_size=label_text_size,
		annot_text_size=annot_text_size, no_align_labels=no_align_labels,
		brlen_scale=brlen_scale, no_fade_ancestral=no_fade_ancestral)

	add_cartoons(tree, node_locs, ax, array_dict, non_singleton_spacers, spacer_cols_dict,
		label_pad=label_pad, spacer_size=spacer_size, spacer_outline=outline,
		spacer_spacing=spacing, spacer_width=spacer_width,
		no_align_cartoons=no_align_cartoons, emphasize_diffs=emphasize_diffs,
		v_scaling_factor=v_scaling_factor, annot_text_size=annot_text_size,
		no_fade_ancestral=no_fade_ancestral, label_pad_dict=label_pad_dict)

	plt.axis('off')
	if branch_support:
		import matplotlib
		plt.rcParams.update({'font.size': label_text_size})
		fig.colorbar(
			matplotlib.cm.ScalarMappable(
				matplotlib.colors.Normalize(
					vmin=0,vmax=100),
				cmap="cividis_r"),
			ax=ax,
			location = "bottom",
			pad=0.02,
			label="Node support (%)",
			shrink=0.5)
	plt.tight_layout()
	plt.margins(0,0)
	fig.set_size_inches(fig_w,fig_h)

	plt.savefig(filename, dpi=dpi)


def plot_diffplot(array_dict, array_order, imp_spacers, spacer_colours,
	text_size=10, plot_width=3, plot_height=3, dpi=300,
	outfile="diffplot.png", line_width=1, connection_outline=False):

	outline = 0.3 # Thickness of lines
	spacer_spacing= 0.5
	spacer_size = 2

	vh_ratio = calc_vh_ratio_and_label_pad_diffplot(array_dict, spacer_size+spacer_spacing,
	spacer_spacing, (0.5, len(array_dict)+0.5), text_size=text_size,
	fig_w=plot_width, fig_h=plot_height)

	w_h_ratio = plot_width/plot_height
	vh_ratio*=w_h_ratio
	spacer_width = min(spacer_size*vh_ratio, 0.5)

	fig, ax = plt.subplots(figsize=(plot_width,plot_height))

	for y, array in enumerate(array_order):
		# Add array name
		ax.text(0.5, y+1, array, ha='left', va='center_baseline', 
			fontsize=text_size, color='black')

		# plot spacers
		for n, spacer in enumerate(reversed(array_dict[array])):
			if spacer in imp_spacers:
				spacer_w = spacer_width
				spcolour = spacer_colours[spacer]
			else:
				if spacer in spacer_colours:
					spcolour = spacer_colours[spacer]
				else:	
					spcolour = ("#000000", "#000000") #black
				spacer_w = 0.25*spacer_width
				outline = 0.2

			plot_spacer(ax, 0, y+1, spacer_size, spacer_w,
				spcolour, outline, n, spacer_spacing, vh_ratio)
			outline = 0.3
		if array != array_order[-1]:
			# plot lines between shared spacers
			draw_diffplot_lines(ax, y, array_order, array_dict, spacer_colours,
				spacer_size, spacer_width, spacer_spacing, outline, vh_ratio,
				line_width=line_width, connection_outline=connection_outline)


	plt.ylim(0.5,len(array_dict)+0.5)
	plt.axis('off')
	plt.tight_layout()
	# plt.show()
	plt.savefig(outfile, dpi=dpi)

