################################
CRISPR Comparison Toolkit (CCTK)
################################

**CCTK** is a Python toolkit for the identification and comparison of CRISPR arrays.
It includes tools to go from assemblies to publication-quality images and is built around simple file formats to allow users to easily fit CCTK tools into existing workflows.

See the :doc:`usage` section for further information, including
:ref:`Installation` instructions.

.. _tools:

**********
CCTK Tools
**********

* `blast <content/blast.html>`_ - Find CRISPR arrays in assemblies using BLASTn
* `minced <content/minced.html>`_ - Find CRISPR arrays in assemblies using minced 
* `crisprdiff <content/CRISPRdiff.html>`_ - Produce a CRISPRdiff plot comparing CRISPR arrays
* `crisprtree <content/CRISPRtree.html>`_ - Perform a maximum parsimony analysis on CRISPR arrays
* `constrain <content/constrain.html>`_ - Predict array relationships constrained by a tree
* `network <content/network.html>`_ - Produce a network representation of spacer sharing among arrays
* `evolve <content/evolve.html>`_ - Perform in silico evolution of CRISPR arrays
* `spacerblast <content/spacerblast.html>`_ - BLAST spacers against a BLASTdb, process output & check for PAMs


.. toctree::
   :caption: Contents:
   :glob:
   :maxdepth: 2

   usage
   Contact Us <contact>

   content/*
