minced
======

Introduction
------------

CCTK minced identifies CRISPR arrays in assemblies using the tool `minced <https://github.com/ctSkennerton/minced>`_. The output produced by minced is then processes into several output files containing information about identified CRISPR arrays and spacers.

.. _minced-basic:

Basic Usage
-----------

``cctk minced`` will identify CRISPRs in all assemblies contained in a directory and process the output using the following minimal command:

.. code-block:: shell

	cctk minced -i <directory of assemblies> -o <output directory> -m -p

**N.B.** Basic usage may include ``-l`` and ``-r`` options as well depending on your situation. See the :ref:`minced-advanced` section for details. 

.. _minced-output:

Output files
------------

Both ``cctk minced`` and ``cctk blast`` produce the same output files. However, ``cctk minced`` also produces output files made my minced that are retained for the user's reference.

``cctk minced`` creates two directories in the user-specified output directory (if they did not already exist): MINCED_OUT/ and PROCESSED/. MINCED_OUT/ contains all of the output files produced by minced. PROCESSED/ contains files summarizing various features of the CRISPR arrays identified by minced in more useful formats than the raw minced output. In addition to creating output files, a summary of the run is written to stderr stating the number of unique spacers that were identified and the number of distinct arrays in which they were identified.

.. _crispr-spacers:

CRISPR_spacers.fna
^^^^^^^^^^^^^^^^^^

**Summary**

Sequences of all of the unique spacers that were identified in the provided assemblies. Fasta headers are constructed from the best matching repeat ID and an integer that counts the number of spacers found associated with that repeat. See :ref:`minced-advanced` for details about how the repeat ID to assign to each spacer is chosen and how orientation of the spacer sequence is chosen by ``cctk minced``. ``cctk blast`` gets the repeat ID portion of the spacer fasta header from the fasta header of the repeat used to identify the CRISPR array by BLAST.

**Format**

Fasta nucleotide sequence

**Example**

.. code-block:: shell

	>1F_1
	GGTACGTGGTTTCGACCAACAGCACTGCCCAA
	>1F_2
	AGGCTGCCAAGTCGGTGCGCGAGGCCGGCTTT
	>1F_3
	TGCAGCGATTGCACCTTGGCCTGCTGCCGATC
	>1E_1
	CATCTGGCCGGGGCTCGGGTCTGGTTCTACGA
	>1E_2
	GATGGCAACCGGCGTTTGTCCGCGCTGAACTG

.. _array-ids:

Array_IDs.txt
^^^^^^^^^^^^^

**Summary**

CRISPR arrays are defined as distinct by ``cctk`` if they contain a single different spacer to any other arrays. Each distinct array is assigned a numerical ID based on the order in which they are found in input assemblies. This file contains the IDs of the spacers in each array.

**Format**

2 columns, tab-delimited.

Column 1: ID of array
Column 2: Space-delimited list of IDs (fasta headers) of spacers in this array

**Example**

.. code-block:: shell

	1	1F_42 1F_18 1F_153 1F_53 1F_82 1F_148
	2	1E_90 1E_56 1E_166 1E_26 1E_141 1E_13 1E_77

.. _array-seqs:

Array_seqs.txt
^^^^^^^^^^^^^^

**Summary**

This file contains the sequence of the spacers in each array.

**Format**

2 columns, tab-delimited.

Column 1: ID of array
Column 2: Space-delimited list of sequence of spacers in this array

**Example**

.. code-block:: shell

	1	GGTACGTGGTTTCGACCAACAGCACTGCCCAA AGGCTGCCAAGTCGGTGCGCGAGGCCGGCTTT 
	2	CATCTGGCCGGGGCTCGGGTCTGGTTCTACGA GATGGCAACCGGCGTTTGTCCGCGCTGAACTG

.. _array-locations:

Array_locations.bed
^^^^^^^^^^^^^^^^^^^

**Summary**

Contig names and contig locations in which CRISPR arrays were identified.

**Format**

BED format.

First line is a "#" character followed by tab-delimited column names.

Name column contains the ID of the array at the indicated location. This ID corresponds to the IDs in :ref:`array-ids` and :ref:`array-seqs`

**Example**

N.B. when viewing this file in a text editor, the headings and column contents will usually not line up, visually. If you wish to view this file for manual inspection, it will read into excel with proper column assignments or can be viewed in the terminal using ``column -t Array_locations.bed | less``

.. code-block:: shell

	#contig             contigStart  contigEnd   name   score   strand
	Assembly1_contig2   208444       209013      6      0       -
	Assembly1_contig6   19991        20559       7      0       +
	Assembly2_contig1   29424        30050       11     0       -

.. _array-reps:

Array_representatives.txt
^^^^^^^^^^^^^^^^^^^^^^^^^

**Summary**

This file indicates which assemblies each array was found in. e.g.,

**format**

2 columns, tab-delimited.

Column 1: ID of array
Column 2: Space-delimited list of sequences in which this array was identified

**Example**

.. code-block:: shell

	array1	assembly1
	array2	assembly2 assembly3
	...


.. _array-network:

Array_network.txt
^^^^^^^^^^^^^^^^^

**Summary**

Network representation of the number and proportion of spacers that arrays have in common with one another. Each pair of arrays that share one or more spacers are respresented by an edge in the network. The similarity between arrays is represented as both the number of spacers in common, and the Jaccard similarity index of the two arrays. The repeat ID associated with each array is also included.

This file can be easily read into a network visualization software such as cytoscape, as demonstrated in the `tutorial <tutorial.html>`_.

Jaccard similarity between two arrays is defined as the number of unique spacers in common between the two arrays, divided by the combined number of unique spacers present in the two arrays. 

e.g. for the following 2 arrays (as they would be represented in Array_IDs.txt):

.. code-block:: shell

	Array	Spacers
	1	1F_1 1F_2 1F_3
	2	1F_4 1F_2 1F_3

The array both contain spacers 1F_2 and 1F_3, while each array also contains one spacer that is not present in the other array. Therefore, the 2 shared spacers are 1F_2 and 1F_3, while the list of 4 total unique spacers in the two arrays is 1F_1, 1F_2, 1F_3, and 1F_4. This results in a Jaccard similarity index of 2/4 = 0.5

Jaccard is an effective similarity measure for comparing CRISPR arrays as it takes into account both the number of spacers in common between two arrays, and the spacers present in each array that are not shared.

**Format**

Tab-delimited.

First line is header information

**Example**

.. code-block:: shell

	Array_A	Array_B	Shared_spacers	Jaccard_similarity	Array_A_type	Array_B_type
	6	4	9	0.75	1F	1F
	11	1	10	0.5263157894736842	1F	1F
	13	8	1	0.02127659574468085	1F	1F
	2	9	12	0.3333333333333333	1F	1F

.. _array-clusters:

Array_clusters.txt
^^^^^^^^^^^^^^^^^^

**Summary**

Clusters of arrays are identified within the network represented in :ref:`array-network`. A cluster of arrays is defined as a set of arrays in which each array shared at least N spacers with one or more other members of the set, where N is the number provided by the user with the option ``--min-shared``.

This file is provided so that you can easily analyse a cluster of arrays using `CCTK CRISPRdiff <CRISPRdiff.html>`_ or `CCTK CRISPRtree <CRISPRtree.html>`_. Instead of typing out the list of arrays, you can copy it from this file or iterate over the lines of this file to analyze every identified cluster.

**Format**

Space-delimited

Each line is a different cluster. Arrays only appear in one cluster each. Arrays that are not present in any cluster of two or more arrays are not present in the file.

**Example**

.. code-block:: shell

	49 67 71 70 13 24 73 4
	77 15
	14 47 57 76 50 58 56

.. _crispr-sum-csv:

CRISPR_summary_table.csv
^^^^^^^^^^^^^^^^^^^^^^^^

**Summary**

Summary of CRISPR arrays found in each assembly with information about each array. This file is designed to be read into Microsoft Excel or a similar program to view.

**Format**

comma-delimited (csv) table

Columns:

#. Sequence_ID: Name of assembly (extracted from input file name)
#. Has_CRISPR: Boolean whether and CRISPR arrays were found
#. Array_count: Number of CRISPR arrays found. No further columns are populated if no arrays were found.
#. Spacers: List of spacer sequences found in each array. Corresponds to sequences in :ref:`crispr-spacers`
#. Spacer_IDs: List of spacer IDs found in each array. Corresponds to IDs in :ref:`crispr-spacers`
#. Array_IDs: List of array IDs corresponding to :ref:`array-ids` and :ref:`array-seqs` files
#. Array_locations: List of array locations (contig name, start, stop)
#. Repeat_sequences: Sequence of the most common repeat in each array
#. Array_CRISPR_types: Most similar repeat type found
#. Array_repeats: Array repeat ID corresponding to sequences in :ref:`array-repeats`
#. Array_repeat_score: Most similar CRISPR repeat to the repeat in this array and the number of mismatches. Format: repeat(mismatches)

In columns 4-11, arrays are numbered according to the order in which they were found in the input assembly file. These numbers correspond between columns in a given row such that the spacer IDs for array 1 correspond to the spacer sequences of array 1 etc.


**Example**

.. image:: images/cr_sum_tab.png

.. _crispr-sum-txt:

CRISPR_summary_table.txt
^^^^^^^^^^^^^^^^^^^^^^^^

**Summary**

Summary of CRISPR arrays found in each assembly with information about each array. This file is easier to interact with programatically.

**Format**

Tab-delimited table with "|" (pipe)-delimited lists of arrays in columns 4-11 within each array, elements are space-delimited.

**Example**

.. code-block:: shell

	Sequence_ID	Has_CRISPR	Array_count	Spacers Spacer_IDs	Array_IDs	Array_locations Repeat_sequences	Array_CRISPR_types	Array_repeats	Array_repeat_score
	Assembly1	True	3	TAGCTGATCAGCAGGCCGACAGTCAGGCCTGC TACCCGAATACGACTTGCGCGAGGAAGACGGT AGCATCGCATCAAATCGTGCAGAACACGATAA TGGTCGAGCAGTTCGGCAAAGGGGCCGTGGTT TTCACCTGGTCGCCGGCCAGGCTGATCACTGC TACAAGGTCATGGCGCTCGGCAACGTGGTGGAA GCTGTGCGTCGCCGTGGTCTGACGGTCGAATC AGCAGATACCCGAACCACTGGAGGTACATGCA TTCATCAGGATGCCGCCAAGGGTCCGCATAAT|AGGTCGAGGTGGGCTCGGCGGCGATGATCGAT GGTACGTGGTTTCGACCAACAGCACTGCCCAA TAAAGGAGATTGCCATGCTGATCAAACTTCCC GTCAGGGTCGTGCATGACTCCGATGTGGTGGC CGTCCAGAACGTCACACGCTCGCCGTCGATGT AACCGGAGCCTTCGGGCCGCGTTGGGATCCAC TTGACTGCTGGGGCCTGACGCTCATCGCGCGG GCGACCCTGGCCAGGGCGGCGTCGCGCTCTGC TTGAGCACAACCGGCTGAGCCAGCTGGTTGTC|CAGCAGCGGCTCCAGGAAGAGGGGCGCTGCCT AAGAGTCGCGGCGACAACTACCAGACGTCCGC GTATGGCTCTCTCCATTGGGGTGGCGATACTC GATCTGGGGCGGCATCATCACAGCAGAATCTA ACAACATCAATCGCCTGATGCTGGGGCACCTG AGCTTCGGCACCCTGATGCGCGCCGTCGAGGG AATGCGGTCCTGCGCATCCGAACTGGTAAGTG GACCCCCGGAGGACCAACCGTGGACAACGACA TCCTTCGGCTCCGCCGGCCGGATCGCTGCAT GTCGCGAAGTTCATAAGCGGGCTTAGGGCGA	1F_156 1F_19 1F_46 1F_123 1F_59 1F_64 1F_34 1F_93 1F_33|1F_99 1F_1 1F_45 1F_83 1F_124 1F_126 1F_30 1F_39 1F_49|1F_134 1F_81 1F_55 1F_84 1F_16 1F_5 1F_51 1F_100 1F_106 1F_145 6|7|11	Assembly1_contig2 209013-208445|Assembly1_contig4 19992-20559|Assembly1_contig4 30050-29425 GTTCACTGCCGTATAGGCAGCTAAGAAA|GTTCACTGCCGTGTAGGCAGCTAAGAAA|GTTCACTGCCGTATAGGCAGCTAAGAAA	1F|1F|1F	1_a|2_a|3_a	1F(2)|1F(1)|1F(0)

.. _array-repeats:

Array_repeats.txt
^^^^^^^^^^^^^^^^^

**Summary**

Sequences of repeats in each identified CRISPR array. If multiple repeat sequence variants were found associated with the same Array ID, then these are indicated using letters (e.g., 1a, 1b when repeat variants are identified in array 1; in the below example the second repeat in 2_a starts GTGTTCCCT... instead of GTGTTCCCC...).

These repeats can be visualized using `CCTK CRISPRdiff <CRISPRdiff.html>`_ just like with spacers (although you may want to set ``--line-width`` to 0).

**Format**

2 columns, tab-delimited.

Column 1: ID of array repeats
Column 2: Space-delimited list repeat sequences in this array

**Example**

.. code-block:: shell
	
	1_a GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG
	2_a	GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCTACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG
	2_b	GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG GTGTTCCCCACATGCGTGGGGATGAACCG 

.. _spacer-cluster-reps:

Spacer_cluster_members.txt
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Summary**

When running ``cctk minced`` with ``-s`` to cluster similar spacers, this file is produced to provide details of which spacers were identified as similar to one another.

**Format**

Tab-delimited table with two columns. Each line represents a distinc cluster of spacers. Column 1 is the ID of the spacer chosen as the representative of the cluster. The ID (or its corresponding sequence - see CRISPR_spacers.fna) is used to represent all cluster members in any files in which they are described. Column 2 is a space-delimited list of the sequences of spacers that are members of the cluster (not including the sequence of the spacer chosen as the representative.)

**Example**

.. code-block:: shell

	CRtype_1	GCCCAGGCACGTTTGCTCGCGCTTTGATCTCA
	CRtype_13	TGTCCCGAAGTTCATAAGCGGGCTTCGGGCGA GTCGCGAAGTTCATAAGCGGGCTTCGGGCGA
	CRtype_42	AGCCGATGGCCCGCAGTAGTACCCCGATCAGT

.. _minced-advanced:

Advanced Usage
--------------

The usage of ``cctk minced`` described in the :ref:`minced-basic` is sufficient to identify CRISPR arrays in assemblies. The most likely situations in which you will need more complex usage of ``cctk minced`` are:

Specifying the path of your MinCED installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you installed minced manually and it is not in your path you can specify the path to MinCED in your ``cctk minced`` command using the ``-l`` option. 

**N.B.** This is not a problem if you install using conda.

e.g.
.. code-block:: shell

	cctk minced -i <directory of assemblies> -o <output directory> -l <path to minced> -m -p

Specifying the CRISPR types of repeats in your assemblies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cctk minced`` has a default database of CRISPR subtype repeats representing variants of the following subtypes: I-A, I-B, I-C, I-D, I-E, I-F, I-G, II-A, II-B, II-C, III-A, and III-B. If you are analyzing assemblies that have CRISPR systems of other subtypes (or variants not in the default database), you will want to specify the repeat (in the correct orientation relative to the leader end) here to ensure that your CRISPR arrays are correctly oriented and categorized.

**N.B.** The repeats in the default database may be added to in the future. (Please send characterized CRISPR repeats with known orientation via email (crisprtoolkit@gmail.com) or as an issue on the `CCTK github <https://github.com/Alan-Collins/CRISPR_comparison_toolkit>`_ and I will be happy to add them. If you have a citable reference for the repeat and its correct orientation all the better!)

``cctk minced`` uses repeats to add CRISPR type information to spacer fasta headers, but also (and more importantly) to figure out the correct orientation of CRISPR arrays with regards to their leader and trailer ends as minced does not check array orientation itself. This information is essential if you wish to analyze your CRISPR arrays using ``cctk crisprtree``.

Relying on the built-in repeat sequences will result in consistent orientation of CRISPR arrays with the same repeat sequence. However, there is a roughly 50% chance your arrays will be output in the reverse orientation.

If you wish to provide your own repeat sequences in order to properly characterize repeat type and orient your arrays correctly, you can provide any number of repeats in fasta format using the ``-r`` option. It is important that your repeat sequences be oriented so that the leader end of the array is 5' of the repeat.

Repeats are only used during processing steps so you do not need to run minced again if you have already done so (i.e. you do not need ``-i`` or ``-m`` inputs). An example command to include user-defined repeats is 

.. code-block:: shell

	cctk minced -o <output directory> -r <repeats file> -p

.. _minced-append:

Appending to an existing dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have previously run ``cctk blast`` on some assemblies and wish to add CRISPR arrays from additional assemblies to the existing dataset, you can do this using ``--append``.

``cctk minced`` and ``cctk blast`` both have an ``--append`` option which can be set to activate the appending mode. In this mode, existing CRISPR information files are read and then their data are added to. **N.B** The files are overwritten in the process so make sure to duplicate your files to another location if you wish to preserve them.

Append mode expects to find existing files within the directory structure created by CCTK (i.e., in the "PROCESSED" directory at the path specified using ``-o``). You do not need to provide all files at this location, but must at least provide a CRISPR_spacers.fna file with your spacers in fasta format. Any existing files will be read to initialize the dataset (for example, spacer ID and array ID assignments). Any files that are absent will simply not be used to initialize a dataset to be added to and will instead be created as if you were not using append mode.

Manually curating MinCED output upstream of CCTK processing steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may wish to manually curate the output of ``cctk minced``. For example, you may find that arrays in your output files seem like they are not actually CRISPR arrays. You may also find that an array has the same one or two bases on the end of every spacer (see :ref:`minced-limitations` for an explanation of how this may occur).

``cctk minced`` does not include functionality for fine control over outputs or how arrays are identified. Instead you must laboriously modify the minced output files. However, while ``cctk minced`` won't help you with this process, it does retain all the minced output files in the MINCED_OUT/ directory in your specified output directory. Furthermore it will allow you to process the modified minced output files without rerunning minced by omitting the ``-m`` flag in your command as in the example below.

.. code-block:: shell

	cctk minced -o <output directory containing MINCED_OUT/> -p

When running only processing steps ``cctk minced`` will read and process all files in the MINCED_OUT/ directory in your specified output directory (Crucially not the input directory specified with ``-i``, but instead the output directory specified with ``-o``). The only requirement is that the format of the minced output files is not changed. You can delete whole arrays from these files and can modify the sequence of spacers and repeats and ``cctk minced`` should process them without issue.

Consider the following example minced output file. 

.. code-block:: shell
	
	Sequence 'Assembly1_contig1' (209122 bp)

	CRISPR 1   Range: 208445 - 208593
	POSITION	REPEAT				SPACER
	--------	----------------------------	--------------------------------
	208445		AAAAAAAAAAAAAAAAAAAAAAAAAAAA	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	[ 28, 32 ]
	208505		AAAAAAAAAAAAAAAAAAAAAAAAAAAA	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	[ 28, 32 ]
	208565		AAAAAAAAAAAAAAAAAAAAAAAAAAAA	

	--------	----------------------------	--------------------------------
	Repeats: 3	Average Length: 28		Average Length: 32

	Time to find repeats: 8 ms


	Sequence 'Assembly1_contig2' (84619 bp)

	CRISPR 2   Range: 19992 - 20319
	POSITION	REPEAT				SPACER
	--------	---------------------------	--------------------------------
	19992		TTCACTGCCGTGTAGGCAGCTAAGAAA	AGGTCGAGGTGGGCTCGGCGGCGATGATCGAT	[ 27, 32 ]
	20052		TTCACTGCCGTGTAGGCAGCTAAGAAA	GGTACGTGGTTTCGACCAACAGCACTGCCCAAG	[ 27, 33 ]
	20112		TTCACTGCCGTGTAGGCAGCTAAGAAA	TAAAGGAGATTGCCATGCTGATCAAACTTCCCG	[ 27, 33 ]
	20172		TTCACTGCCGTGTAGGCAGCTAAGAAA	GTCAGGGTCGTGCATGACTCCGATGTGGTGGCG	[ 27, 33 ]
	20232		TTCACTGCCGTGTAGGCAGCTAAGAAA	CGTCCAGAACGTCACACGCTCGCCGTCGATGTG	[ 27, 33 ]
	20292		TTCACTGCCGTGTAGGCAGCTAAGAAA	
	--------	---------------------------	--------------------------------
	Repeats: 6	Average Length: 27		Average Length: 33

In this example file the first array is clearly nonsense, while the second array has what looks like a type I-F repeat missing the first G and most of the spacers have a G on one end. It seems like the first array should be removed, while the second array should be modified to correct the mischaracterization of the repeat boundaries.

In minced output files, the information about a CRISPR array begins on the line starting with the word "CRISPR" and ends on the line starting with the word "Repeats". In addition, If multiple arrays are identified in the same contig, they will have a single line starting with the word "Sequence" that identifies all of the subsequent arrays as being found in the names contig.

If you wish to delete an array, remove all lines describing that CRISPR array. If it is the only array found in that contig, remove the line above it starting with "Sequence" as well.

Modifying repeat and spacer sequences is easier. Just make the desired changes. You do not need to change the length information on the right of each line as ``cctk minced`` does not use that information. Additionally, you do not need to worry about the number of blank lines.

Making the above changes would result in the following modified file:

.. code-block:: shell

	Sequence 'Assembly1_contig2' (84619 bp)

	CRISPR 2   Range: 19992 - 20319
	POSITION	REPEAT				SPACER
	--------	---------------------------	--------------------------------
	19992		GTTCACTGCCGTGTAGGCAGCTAAGAAA	AGGTCGAGGTGGGCTCGGCGGCGATGATCGAT	[ 27, 32 ]
	20052		GTTCACTGCCGTGTAGGCAGCTAAGAAA	GGTACGTGGTTTCGACCAACAGCACTGCCCAA	[ 27, 33 ]
	20112		GTTCACTGCCGTGTAGGCAGCTAAGAAA	TAAAGGAGATTGCCATGCTGATCAAACTTCCC	[ 27, 33 ]
	20172		GTTCACTGCCGTGTAGGCAGCTAAGAAA	GTCAGGGTCGTGCATGACTCCGATGTGGTGGC	[ 27, 33 ]
	20232		GTTCACTGCCGTGTAGGCAGCTAAGAAA	CGTCCAGAACGTCACACGCTCGCCGTCGATGT	[ 27, 33 ]
	20292		GTTCACTGCCGTGTAGGCAGCTAAGAAA	
	--------	---------------------------	--------------------------------
	Repeats: 6	Average Length: 27		Average Length: 33


.. _minced-limitations:

Limitations and considerations
------------------------------

Minced uses a sliding window to detect regions containing more than two (roughly) equally spaced (approximately) repeated sequences. The first two repeated sequences that are found (as the window slides 5' to 3' along the sequence) are used to define the repeat sequence. Additional windows are then added, the same distance apart until no more repeats are found. See the `CRT publication <https://doi.org/10.1186/1471-2105-8-209>`_ for further description. This approach results in a few behaviours that a user should bear in mind:

* Not all regions containing 3 or more (approximate) repeats are CRISPRs. Manual curation is important to confirm that predicted CRISPR arrays are to be believed.

* By only comparing a few, short sequences (i.e. the contents of the sliding windows), minced tolerates relatively large numbers of differences between repeats while still being confident the sequences are related. This can result in the inclusion of spacers flanked by fairly degenerate repeats.

* By determining the repeat sequence using the first repeats encountered, minced is vulnerable to mischaracterizing the repeat sequence in the rest of the array if these first repeats are degenerate. CRISPR array trailer repeats often contain mutations not present in more leader-proximal repeats. If minced finds an array encoded on the minus strand (i.e. it finds the trailer end first while scanning the plus strand 5' to 3') and the array has differences in it's trailer-most repeats end-most bases, this can result in minced miscalling the boundaries of the repeat and including one or two repeat bases in all (or most) spacers in the array.
