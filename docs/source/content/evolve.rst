######
evolve
######

************
Introduction
************

``cctk evolve`` is a tool designed to perform *in-silico* evolution of CRISPR arrays. It was developed to provide test datasets to evaluate the performance of ``cctk crisprtree`` and is packaged as part of CCTK to allow the generation of test or example datasets for use with CCTK tools.

.. _evolve-process:

How does cctk evolve work
=========================

``cctk evolve`` performs *in silico* evolution of CRISPR arrays according to the following process in which spacers and arrays are represented simply by integers:

1. An initial array of length set using ``-i`` is created with spacers numbered 1 through whatever number is set (default: 5). This array is assigned the number 0.

2. Based on the relative frequencies of acquisition, deletion, and trailer loss set uing ``-a``, ``-d``, and ``-t``, respectively, an event is chosen at random (defaults are 75, 10, and 15, respectively).

3. The array is duplicated and the chosen event is applied to the copy of the array using the following process for the different event types:
	
  * Acquisition: A new spacer is created and assigned the next unused integer. This spacer is added to the leader end of the array.
  
  * Trailer-loss: the spacer at the trailer end of the array is removed.
  
  * Deletion: Deletions are performed using a multi-step process.
  	
     #. First the start and stop of the deletion region are selected by chosing a random number according to a probability distribution that follows a normal distribution with a mean equal to the length of the array divided by 2 and a standard deviation equal to the length of the array divided by 4. Numbers are chosen in this way until two different numbers are selected that would not result in the deletion of every spacer in the array.
    
     #. The spacers between the indices described by the selected numbers are removed from the array.

4. It is next decided whether to delete the original array according to a percent probability set using ``-l`` (default: 50%).

5. The newly created array is set as a descendent of the parent array in a tree object. This tree preserves information about which array gave rise to which other arrays

6. Of the remaining arrays (i.e., those that have not been removed in step 4), an array is selected at random and the process is repreated from step 2 until the number of events specified with ``-n`` has been performed.

7. Once all events have been performed, the tree object is processed to remove any nodes in the tree that only have a single descendant. Thus, each node in the tree is either a leaf of an ancestral array that gave rise to more than one descendant.

.. _evolve-basic:

***********
Basic Usage
***********

``cctk evolve`` requires a single input: the number of events (``-n``) for which the evolution should be allowed to continue. e.g.:

.. code-block:: shell

	cctk evolve -n 20

In addition, reproducable simulations can be performed by using the seed option ``-s``. All examples described on this page are run with the following command using a seed value of 1. e.g.:

.. code-block:: shell

	cctk evolve -n 20 -s 1

************
Output files
************

All ``cctk evolve`` output files are named according to a fixed convention and can not be set by the user at run time. File names are composed of the data type followed by the settings used in that run. 

For example, for the file "evolved_arrays_5_20_75_10_15_50_1.txt" the parts of the filename are as follows.

"evolved_arrays" - indicates that this file contains the arrays generated during the *in silico* evolution process. This file is equivalent to the :ref:`array-ids` file used by other CCTK tools.

"5_20_75_10_15_50_1" - these numbers correspond to the following settings respectively:
	
* ``-i`` initial array length
* ``-n`` number of events
* ``-a`` acquisition rate
* ``-d`` deletion rate
* ``-t`` trailer-loss rate
* ``-l`` array-loss rate
* ``-s`` seed used for random processes

.. _evolve-treeplot:

Tree plot
=========

The tree plot produced by ``cctk evolve`` is much like that produced by ``cctk crisprtree`` and ``cctk constrain``. (See  :ref:`tree-plot`.) However, as shown in the below image, the plots are not identical.

.. image:: images/eg_evolve.png

In the above image, ancestral arrays are slightly transparent just like in other CCTK tree plots. However, instead of being assigned IDs such as "Anc a", ancestral arrays in tree plots produced by ``cctk evolve`` are assigned numbers.

The reason for this difference is that ancestral arrays in tree plots produced by ``cctk evolve`` are not inferred based on extant arrays. Instead these arrays are the true ancestral arrays of their descendents in the tree (as described in the :ref:`evolve-process` section).

The numbers assigned to each array in the tree correspond to the event number that created them. For example, in the above image, array 0 was the initial array that is the last common ancestor of all others. Array 1 was created during the first event in the simulation.

.. _evolve-arrayfile:

Evolved arrays
==============

``cctk evolve`` produces a file with a name beginning "evolved_arrays". This file is the equivalent of the :ref:`array-ids` file produced by other CCTK tools. It contains the arrays and spacers produced during the simulation.

.. _evolve-colourfile:

color_scheme
============

``cctk evolve`` produces a "color_scheme" file for every run. This file is the same as that described for other CCTK tools such in the :ref:`CRISPRtree-json` section of the `CRISPRtree <CRISPRtree.html>`_ documentation page.

.. _evolve-treefile:

Tree file
=========

``cctk evolve`` writes a newick representation of the produced tree to a file. This tree corresponds to the tree shown in the :ref:`evolve-treeplot` file.

.. _evolve-advanced:

**************
Advanced Usage
**************

Changing evolution parameters
=============================

The following evolution parameters can be changed using command-line options to control the simulation of CRISPR array evolution:

* ``-i``, or ``--initial-length``  length of the starting array. Default = 5
* ``-a``, or ``--acquisition``     relative frequency of spacer acquisitions. Default = 75
* ``-t``, or ``--trailer-loss``    relative frequency of trailer spacer decay. Default = 15
* ``-d``, or ``--deletion``        relative frequency of deletions . Default = 10
* ``-l``, or ``--loss-rate``       rate arrays are lost after spawning descendant. Default = 50%

The impact of the above parameters is explained in the :ref:`evolve-process` section.

Controlling plot elements and size
==================================

Plot element control
--------------------

Several visual elements of the plot produced by ``cctk evolve`` can be controlled using command line options.

The default behaviour of ``cctk evolve`` is to de-emphasize ancestral arrays by applying transparency to their node labels and array cartoons. This can be disabled using the ``--no-fade-anc`` option.

The default behaviour of ``cctk evolve`` is to annotate hypothetical events onto arrays. This can be disabled using the ``--no-emphasize-diffs`` option.

The inclusion of branch length annotations can be controlled using ``-b``. Branch lengths correspond to the weighted parsimony cost of events between an array and its ancestor. Branch length labels are added at the midpoint of the corresponding branch.

Branch lengths can be scaled by a (floating point number) factor provided using ``--brlen-scale``. This can be used to increase or decrease all branch lengths. Horizontal space taken up by branches in the tree reduces the space available for CRISPR array cartoons so this option can be used to control the amount of space in the plot used by those two components.

The default behaviour of ``cctk evolve`` is to align node labels and array cartoons.The alignment of both cartoons and labels can be deactivated using ``--no-align``.

Plot size and resolution
------------------------

The size and resolution of plots produced by ``cctk evolve`` can be controlled using command line options. These options can be used to generate images of the exact specification required for a figure, or may be necessary to create a sensibly scaled image (see :ref:`tree-limitations`).

Plot height and width can be set using the options ``--plot-width`` and ``--plot-height`` and providing the desired size in inches.

pixel density (DPI) can be set using ``--dpi``. The images on this page were generated at 600 DPI. **N.B.** DPI settings are only relevant for images generated by ``cctk evolve`` in raster formats such as PNG. SVG outputs are unaffected by DPI settings.

``--font-override-labels`` and ``--font-override-annotations`` can be used to control the size of text in the plot (default value is 10pt).
