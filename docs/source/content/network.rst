network
=======

Introduction
------------

``cctk network`` performs pairwise comparisons of CRISPR arrays and represents their spacer-sharing relationships in as a network. This functionality is included in both ``cctk minced`` and ``cctk blast`` which output :ref:`array-network`. This network creation functionality is provided separately in ``cctk network`` to allow the creation of networks in, for example, the following situations:
	
#. When analyzing CRISPR arrays that were not identified by CCTK tools.

#. To generate a network of a the relationships between a subset of arrays.

#. To generate a network of arrays that have been modified after their identification by ``cctk minced`` or ``cctk blast`` (e.g. following manual curation).

Before you run
--------------

CCTK is built around the :ref:`array-ids` file (or :ref:`array-seqs` if you prefer) and that file is the only required input for ``cctk network``. If you identified CRISPR arrays using a non-CCTK tool then you will need to format your arrays into the format of the :ref:`array-ids` file. That format is composed of two tab-delimited columns containing the following information:

#. A unique identifier for the array. This can be any string of characters including numbers, letters, and symbols.

#. A space-delimited list of identifiers for the spacers in this array. Spacer identifiers can also be any string of characters including numbers, letters, and symbols. The only requirement is that the identifier used for a spacer is the same between arrays in which it is contained.

**N.B.** While ``cctk minced`` and ``cctk blast`` output :ref:`array-ids` in the format 

.. code-block:: shell

	<Array ID>\t<spacer> <spacer>...

with a tab betewen the two columns, this file can be successfully read with any whitespace delimiters. Array ID must come before the first whitespace character and the spacers must be separated by whitespace, but tabs and spaces will work equally well in all locations. You can use all tabs or all spaces when generating this file if it is easier.

.. _network-basic:

Basic usage
-----------

The minimal command for ``cctk network`` requirs only an :ref:`array-ids` or :ref:`array-seqs` file as input. ``cctk network`` will write a file called :ref:`array-network` to the directory specified using ``-o`` which is your currently directory (``./``) by default.

.. code-block:: shell

	cctk network -i <Array file>

This produces a version of the :ref:`array-network` file produced by ``ctk minced`` and ``cctk blast`` that contains only the first 4 columns. 

.. code-block:: shell

	Array_A	Array_B	Shared_spacers	Jaccard_similarity
	1	2	2	0.2
	1	3	5	0.5
	1	4	8	0.8

See :ref:`network-advanced` for information about producing a file including the remaining columns.

An :ref:`array-clusters` file like that produced by ``ctk minced`` and ``cctk blast`` is also output.

.. _network-advanced:

Advanced Usage
--------------

With the command described in the :ref:`network-basic`, the output file will contain only the first 4 columns contained in the :ref:`array-network` files produced by ``cctk minced`` and ``cctk blast``. This is because ``cctk network`` does not attempt to infer CRISPR type information from your array file. If you want columns describing the CRISPR type of arrays in the network, you need to provide that information as a separate file using ``-t``.

.. _network-array-types:

Array types file
^^^^^^^^^^^^^^^^

The array types file contains two columns. The first is the Array ID used in :ref:`array-ids` or :ref:`array-seqs`, the second is the CRISPR type that you would like written to your :ref:`array-network` file and can be any string of characters you wish (it does not need to correspond to any other information and could, in fact, be any kind of annotation as long as it contains no whitespace characters).

.. code-block:: shell

	1	1F
	2	Type_I_Subtype_F
	3	I-E
	4	Other_annotation

If you are simply making a network from a subset of arrays in an existing :ref:`array-ids` file then you can quickly and easily produce the array types file using ``grep`` and ``sed``:

.. code-block:: shell

	$ cat Array_IDs.txt

	1	1F_42 1F_18 1F_153 1F_53 1F_82
	2	1F_90 1F_56 1F_166 1F_26 
	3	1F_56 1F_166 1F_26 1F_141 
	4	1F_156 1F_19 1F_26 1F_141
	...
	
	$ sed 's/_.*//' Array_IDs.txt | grep -Ew "^1|^2|^3|^4|..." > Array_types.txt

	$ cat Array_types.txt

	1	1F
	2	1F
	3	1F
	4	1F
	...


Complete command
^^^^^^^^^^^^^^^^

We can then use the :ref:`network-array-types` to annotate edges in the output network with information about the connected arrays with the following command (also including output directory specification for completeness):

.. code-block:: shell

	cctk network -i <Array file> -t <array types file> -o <output dir>

If we used the example file described in the :ref`network-array-types` section, then this would result in the following output file:

.. code-block:: shell

	Array_A	Array_B	Shared_spacers	Jaccard_similarity	Array_A_type	Array_B_type
	1	2	2	0.2	1F	Type_I_Subtype_F
	1	3	5	0.5	1F	I-E
	1	4	8	0.8	1F	Other_annotation
