######
evolve
######

************
Introduction
************

``cctk evolve`` is a tool designed to perform *in-silico* evolution of CRISPR arrays. It was developed to provide test datasets to evaluate the performance of ``cctk crisprtree`` and is packaged as part of CCTK to allow the generation of test or example datasets for use with CCTK tools. 

.. _evolve-basic:

Basic Usage
===========

``cctk evolve`` requires a single input: the number of events (``-n``) for which the evolution should be allowed to continue. e.g.:

.. code-block:: shell

	cctk evolve -n 20

Output files
============

All ``cctk evolve`` output files are named according to a fixed convention and can not be set by the user at run time. File names are composed of the data type followed by the settings used in that run. 

For example, for the file "evolved_arrays_5_20_75_10_15_50_1.txt" the parts of the filename are as follows.

"evolved_arrays" - indicates that this file contains the arrays generated during the *in silico* evolution process. This file is equivalent to the :ref:`array-ids` file used by other CCTK tools.

"5_20_75_10_15_50_1" - these numbers correspond to the following settings respectively:
	
``-i`` initial array length
``-n`` number of events
``-a`` acquisition rate
``-d`` deletion rate
``-t`` trailer-loss rate
``-l`` array-loss rate
``-s`` seed used for random processes

.. _evolve-treeplot:

Tree plot
---------

.. _evolve-arrayfile:

Tree plot
---------

.. _evolve-colourfile:

Tree plot
---------

.. _evolve-treefile:

Tree plot
---------


.. _evolve-advanced:

Advanced Usage
==============
