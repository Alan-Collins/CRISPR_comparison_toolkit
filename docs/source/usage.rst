Usage
=====

.. _installation:

Installation
------------

To use CCTK, first install it. The easiest installation method is using ``conda``. CCTK is currently available on my personal Anaconda channel. In the future I intend to host it through Bioconda.

**Conda - Recommended**

.. code-block:: console

   $ conda create -n cctk -c alan-collins cctk
   $ conda activate cctk
   (cctk) $ 

Conda has the benefit of handling the installation of the correct version of all depencies and adds the ``cctk`` (and dependency) executables to your PATH, which makes usage simpler.

**Git clone**

You can also install dependencies seperately and download CCTK from github.

To download CCTK using ``git``:

.. code-block:: console

  $ git clone https://github.com/Alan-Collins/CRISPR_comparison_toolkit.git

You can find the ``cctk`` executable at CRISPR_comparison_toolkit/cctk

**Dependencies**

* `python3 <https://www.python.org/downloads/>`_ (Tested with Python >= 3.8)
   * `Matplotlib <https://matplotlib.org/3.1.1/users/installing.html>`_
   * `DendroPy <https://dendropy.org/downloading.html>`_
   * `NumPy <https://numpy.org/install/>`_

* `minced <https://github.com/ctSkennerton/minced>`_
* `BLAST+ <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ (tested with BLAST 2.9.0)


Basic Usage
-----------

During all examples, we will use CCTK as if it were installed using conda and is in our PATH.

All tools available in CCTK can be accessed through the ``cctk`` executable. Each tool can be executed by calling the corresponding command within ``cctk``. A list of available commands can be found in the help message returned by the ``cctk`` executable when run with the ``-h`` or ``--help`` options.

.. code-block:: console

   (cctk) $ cctk -h
   usage: cctk [-h] [--version]

   optional arguments:
     -h, --help  show this help message and exit
     --version   show program's version number and exit

   Available commands in the CRISPR comparison toolkit:

     Call any command followed by -h or --help to show help for that command

     Find CRISPR arrays in assemblies:
       blast        find CRISPR arrays with user-provided repeat(s) using BLASTn
       minced       find CRISPR arrays using minced

     Analyze the differences between CRISPR arrays:
       CRISPRdiff   produce a CRISPRdiff plot comparing CRISPR arrays
       CRISPRtree   perform a maximum parsimony analysis on CRISPR arrays
       constrain    predict array relationships constrained by a tree ## TO ADD ##
       network      produce a network representation of spacer sharing among arrays

     Other:
       evolve       perform in silico evolution of CRISPR arrays
       spacerblast  BLAST spacers against a BLASTdb, process output & check for PAMs

Usage of each command can be found by calling that command with ``-h`` or ``--help``. e.g. ``cctk blast -h``

Details of the specific usage for each tool in CCTK can be found in the :ref:`Tools` section.
