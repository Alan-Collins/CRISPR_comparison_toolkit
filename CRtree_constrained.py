import dendropy
import sys
from CRISPR_mp import *

intree = sys.argv[1]
array_file = sys.argv[2]
genome_array_file = sys.argv[3]


array_spacers_dict = {}
with open(array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		array_spacers_dict[bits[0]] = bits[2:]

genome_array_dict = {}
with open(genome_array_file, 'r') as fin:
	for line in fin.readlines():
		bits = line.split()
		genome_array_dict[bits[1]] = bits[0]
		if bits[0].lower() == "outgroup":
			outgroup_taxon = bits[1]

tree = dendropy.Tree.get(path=intree, schema="newick")

# If tree not rooted using outgroup node then root it
outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
tree.to_outgroup_position(outgroup_node, update_bipartitions=False)

# Now remove outgroup as it isn't used for array analysis
tree.prune_taxa_with_labels([outgroup_taxon])

print(tree.as_ascii_plot(plot_metric='length'))
