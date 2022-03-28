CRISPRtree
==========

Introduction
------------

``cctk CRISPRtree`` uses a maximum parsimony approach to create a tree representation of the relationship between a set of arrays. It produces a visualization of this hypothesis that indicates how arrays are related to one another and which events may have occurred in each array since an ancestral state. Branch lengths in the tree correspond to the weighted parsimony cost between an array and its hypothetical ancestor. All arrays (both input, extant arrays, and hypothetical ancestral arrays) are depicted in the same style as ``cctk CRISPRdiff``. In addition, any events that are hypothesized to have occurred in an array since its ancestor are indicated using symbols corresponding to the :ref:`tree-key`. The process by which ``cctk CRISPRtree`` constructs trees is discussed in the :ref:`tree-process` section.

.. _tree-before-you-run:

Before you run
--------------

``cctk CRISPRtree`` requirs only an :ref:`array-ids` or :ref:`array-seqs` file as input. By default, all arrays present in the input file will be analyzed. However, ``cctk CRISPRtree`` requires that all input arrays are related to one another (i.e., all arrays are part of a single cluster in a network representation of their relationships) if not all of the arrays are related to one another the resulting plot will be harder to interpret (busy plot, harder to assign visually distinct colours to spacers). It is therefore recommended that you run ``cctk CRISPRtree`` only on smaller batches of your arrays that share spacers.

If you identified CRISPR arrays using ``cctk minced`` or ``cctk blast``, you will have a :ref:`array-network` file among the output of those tools. This file can be visualized using a network visualization tool such as `cytoscape <https://cytoscape.org/download.html>`_ and clusters of related arrays can be selected easily. See the section :ref:`network-tutorial` for an example of how this workflow may look.

.. _tree-basic:

Basic Usage
-----------

``cctk CRISPRtree`` requires two command line inputs: an :ref:`array-ids` (or :ref:`array-seqs`) file using ``-a``.

.. code-block:: shell
	
	cctk CRISPRtree -a <Array_IDs.txt>

The above command will produce a newick string for the most parsimonious tree(s), but will not produce a graphical represenation of the tree. To save a graphical representation of the most parsimonious tree(s) to a file, you must provide a destination filename using ``-o``. e.g.

.. code-block:: shell

	cctk CRISPRtree -a <Array_IDs.txt> -o <output plot with desired extension>

**N.B.** ``cctk CRISPRtree`` uses `matplotlib <https://matplotlib.org/>`_ to perform all plotting functions. You can specify the format of the output file by providing a filename with an extension corresponding to the desired file format. E.g. out_file.png will produce a PNG format file, while out_file.svg will produce an SVG format file. Any file format compatible with `matplotlib.pyplot.savefig() <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html>`_ should work.

Outputs
-------

``cctk CRISPRtree`` produces outputs both to specified file locations (depending on command line options specified) and to the ``stdout`` and ``stderr`` as described below.

.. _tree-plot:

Tree plot
"""""""""

If ``cctk CRISPRtree`` is run with ``-o``, a graphical representation of the most parsimonious tree(s) will be saved to the specified file path. If more than one equally parsimonious tree is found to have the best parsimony score, then they will be saved to separate files which are derived from the specified filename through the addition of '_1', '_2' etc before the file extension. E.g., given the output file 'output/tree_image.png', two equally parsimonious trees would be saved to 'output/tree_image_1.png' and 'output/tree_image_2.png'.

An example tree plot is show below. This tree was generated using the same arrays as were used for the :ref:`diff-output` example plots. i.e.,

.. code-block:: shell

	20      16 13 9 5 4 3 2 1
	19      21 15 12 6 5 4 3 2 1
	16      18 15 12 6 5 4 3 2 1
	7       11 6 5 4 3 2 1
	15      17 8 5 4 1

.. image:: images/CRISPRtree_eg_tree.png

The tree topology on the right of the plot indicates the hypothetical relationship between the plotted arrays. Each input array as well as each hypothetical ancestral array is represented in the same style as that used in the :ref:`diff-output`. Key visual elements of the tree plot, corresponding to the numbers in the above image, are described below:

1. The ID of the array being plotted. Input arrays use the ID that was present in the input file provided using ``-a``.

2. The ID of a hypothetical ancestral array. Ancestral arrays are assigned identifiers beginning with "Anc" (abbreviated from ancestral) and a letter identifier. If more than 26 ancestral states are predicted, then two letter idetnfiers will be used. Ancestral array leaf labels and array cartoons are plotted with slight transparency to provide visual contrast with input arrays.

3. Spacers that are only present in a single input array are represented using a thin, black rectangle. Importantly, the presence of a spacer in an input array AND a hypothetical ancestral array is not sufficient for that spacer to be assigned a colour. The presence of a spacer in an ancestral array is not considered when choosing how to represent the spacer.

4. Spacers that are present in more than one input array are represented using thick rectangles. Each combination of fill and outline colour corresponds to each unique spacer. Note that hypothetical ancestral arrays are depicted using slight transparency to distinguish them from input arrays.

5. The relationships between the arrays being analyzed is represented as a tree. Branch lengths correspond to the weighted parsimony cost of all events that are predicted to have occurred between each array and its hypothetical ancestral array.

6. The root of the tree corresponds to the array that ``cctk CRISPRtree`` hypothesizes is the last common ancestor of all the arrays being considered. The spacers shown in this array are not necessarily all the spacers that would have been present in this array, but are all the spacers for which ``cctk CRISPRtree`` has evidence.

7. Events that ``cctk CRISPRtree`` hypothesizes may have occurred between each array and its hypothetical ancestor are indicated using symbols that correspond to the :ref:`tree-key`. Acquisition and duplication events are indicated with symbols that are placed below the acquired or duplicated spacers. Deletions and trailer loss events are indicated with symbols in-place of the lost spacers. insertions, no identity, and redundant gain events are indicated with a box that surrounds all spacers that are hypothesized to have been involved in the inicated event.


.. _tree-key:

CRISPRtree event key
""""""""""""""""""""

.. image:: images/tree_key.png


stdout
""""""

Newick format tree strings are written to the stdout. If multiple, equally parsimonious trees are found, they will be written to stdout separated with newlines.


stderr
""""""



.. _tree-advanced:

Advanced Usage
--------------



.. _tree-process:

CRISPRtree tree-making process
------------------------------
