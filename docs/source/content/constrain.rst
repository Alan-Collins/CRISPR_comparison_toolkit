#########
constrain
#########

************
Introduction
************

``cctk constrain`` can be used to assess how CRISPR arrays may have evolved given a certain tree topology (i.e., how would CRISPR arrays have changed if their evolution were "constrained" by a certain topology).

This analysis can be useful when ``cctk crisprtree`` produces a topology that differs from that produced using other genomic data. In such a case, ``cctk constrain`` can provide information about whether a tree is still a reasonable explanation of CRISPR array relationships, but simply a less parsimonious topology than that inferred by ``cctk crisprtree``.

``cctk constrain`` can also indicate that horizontal gene transfer may have occurred. ``cctk constrain`` can hypothesize when a topology would require that arrays at different points in the tree independently acquire the same spacers (either by insertion or by leader-end acquisition). Such an event is highlighted in the produced plot and information concerning the event is sent to the ``stderr``.

.. _constrain-before-you-run:

Before you run
==============

In addition to files produced by other CCTK tools such as ``cctk minced`` or ``cctk blast``, you will also need a tree of the relationships between your isolates. Information about which leaf labels in the tree you wish to analyze correspond to which arrays must also be provided. An example dataset and workflow are shown in the `tutorial <tutorial.html>`_.


.. _constrain-basic:

Basic Usage
===========

The basic command to run ``cctk constrain`` requires three input files and the path to an output file:

.. code-block:: shell

	cctk constrain -a <Array_IDs.txt> -t <newick format tree file> -g <array-leaf file> -o <output tree plot>

Further description of the required arguments:

``-a`` File describing the spacers in each array being analyzed e.g. :ref:`array-ids` produced by ``cctk minced`` and ``cctk blast``.

``-t`` Tree in newick format. Tree will be treated as rooted unless ``-u`` is used (see)

``-g`` File corresponding leaf IDs in tree file to array IDs in file provided with ``-a``. Tab-delimited, 1 line per array-leaf ID pair. e.g., (Note array 1 present associated with two leaves on different lines)

.. code-block:: shell
	
	array1	leaf1
	array1	leaf2
	array2	leaf3
	...

.. _constrain-advanced:

Advanced Usage
==============


Output files
============

.. _constrain-treeplot:

Tree plot
---------

The tree plot produced by ``cctk constrain`` is much like that produced by ``cctk crisprtree`` (see :ref:`tree-plot`). However, there are some differences which will be discussed below using this example image (generated using the dataset provided in the `tutorial <tutorial.html>`_).

.. image:: images/large_cluster_constrain.png

.. _constrain-tree-key:

Constrain event key
-------------------

Events that Constrain can hypothesize are shown in the below key. They are the same as those for the :ref:`tree-key`

.. image:: images/tree_key.png


stdout
------

Newick format tree strings and ascii representation of the tree are written to the stdout.


stderr
------

If redundant acquisition events were found involving more than 2 arrays, the Event number and corresponding arrays will be written to stderr.
