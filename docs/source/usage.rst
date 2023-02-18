Usage
=====

.. _installation:

Installation
------------

To use CCTK, first install it. The easiest installation method is using ``conda``. CCTK is currently available on my personal Anaconda channel. In the future I intend to host it through Bioconda.

Conda - Recommended
^^^^^^^^^^^^^^^^^^^

**N.B.** Some dependencies of CCTK are distributed through the bioconda and conda-forge channels. If you do not have those in your conda config you can add them as follows:

.. code-block:: shell

  conda config --append channels conda-forge
  conda config --append channels bioconda

You can then create a conda environment and install CCTK and all dependencies into it with the following:

.. code-block:: shell

  conda create -n cctk cctk
  conda activate cctk

If you do not wish to modify your conda config you will need to specify the above channels in your installation command.

.. code-block:: shell

  conda create -n cctk -c conda-forge -c bioconda cctk

Conda has the benefit of handling the installation of the correct version of all depencies and adds the ``cctk`` (and dependency) executables to your PATH, which makes usage simpler.

Git
^^^

You can also install dependencies seperately and download CCTK from github.

To download CCTK using ``git``:

.. code-block:: shell

  $ git clone https://github.com/Alan-Collins/CRISPR_comparison_toolkit.git

You can find the ``cctk`` executable at CRISPR_comparison_toolkit/cctk

Dependencies
^^^^^^^^^^^^

* `python3 <https://www.python.org/downloads/>`_ (Tested with version >= 3.8)
   * `Matplotlib <https://matplotlib.org/3.1.1/users/installing.html>`_ (Tested with version 3.5.0)
   * `DendroPy <https://dendropy.org/downloading.html>`_ (Tested with version 4.5.2)
   * `NumPy <https://numpy.org/install/>`_ (Tested with version 1.21.2)

* `MinCED <https://github.com/ctSkennerton/minced>`_
* `BLAST+ <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ (tested with version >= 2.9.0)

Quick run
---------

CCTK includes a ``quickrun`` command to run some CCTK tools with default settings to allow you to get a quick look at your data. ``cctk quickrun`` identifies CRISPRs using `CCTK Minced <content/minced.html>`_ and then visualizes all identified groups of related arrays (groups that contain between 3 and 15 arrays by default) using `CCTK CRISPRdiff <content/CRISPRdiff.html>`_ and `CCTK CRISPRtree <content/CRISPRtree.html>`_.

``cctk quickrun`` is a good way to check if your dataset contains any related arrays and to get a sense of their relationships. If you would like to see what a more involved analysis using CCTK would look like then have a look at the `tutorial <content/tutorial.html>`_ page.

Example Workflow
----------------

An example workflow of using all of the tools in CCTK is given in the `tutorial <content/tutorial.html>`_ page.

Usage Help
----------

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
       crisprdiff   produce a CRISPRdiff plot comparing CRISPR arrays
       crisprtree   perform a maximum parsimony analysis on CRISPR arrays
       constrain    predict array relationships constrained by a tree ## TO ADD ##
       network      produce a network representation of spacer sharing among arrays

     Other:
       evolve       perform in silico evolution of CRISPR arrays
       spacerblast  BLAST spacers against a BLASTdb, process output & check for PAMs

Usage of each command can be found by calling that command with ``-h`` or ``--help``. e.g. ``cctk blast -h``

Details of the specific usage for each tool in CCTK can be found in the :ref:`Tools` section.

Adapting your data from a non-CCTK pipeline
-------------------------------------------

CCTK is primarily built around two files: 

#. a fasta file of your spacer sequences, and 
#. :ref:`array-ids`

If you identified CRISPR arrays using a method other than CCTK, then in order to perform analyses using CCTK you will need to generate one or both of these two files. For the tools `CCTK Network <content/network.html>`_, `CCTK CRISPRdiff <content/CRISPRdiff.html>`_, `CCTK CRISPRtree <content/CRISPRtree.html>`_, and `CCTK Constrain <content/constrain.html>`_, you only need the :ref:`array-ids` file. For `CCTK Spacerblast <content/spacerblast.html>`_, you will need a fasta file of your spacers.

How you can convert your data to these CCTK data formats will depend on the method you used to identify CRISPR arrays. However, CCTK does provide one simple way to generate all of the CCTK files if you already have CRISPR spacers in fasta format. Specifically, `CCTK Minced <content/minced.html>`_ and `CCTK Blast <content/blast.html>`_ have an option that allows :ref:`blast-append`. You can use that option using the following steps to run `CCTK Minced <content/minced.html>`_ or `CCTK Blast <content/blast.html>`_:

N.B., We recommend using `CCTK Blast <content/blast.html>`_ for this, assuming you know the sequence of the CRISPR repeats in your assemblies. Minced can make mistakes about repeat boundaries and can identify false-positive arrays. If you know your repeat sequences, `CCTK Blast <content/blast.html>`_ will produce a better output.

#. Make an output directory to store the CCTK output files. (If running `CCTK Minced <content/minced.html>`_, make another directory inside of this one called "PROCESSED")
#. Copy your CRISPR spacers fasta file into the newly created directory (the PROCESSED directory if using `CCTK Minced <content/minced.html>`_ )
#. Run `CCTK Minced <content/minced.html>`_ or `CCTK Blast <content/blast.html>`_ with ``--append`` and whichever other options you wish

You should now have a directory containing all of the normal CCTK output files. The spacer fasta file you provided was read by CCTK and used during the process of identifying CRISPR arrays. In addtion, any identified CRISPR spacers that were not in the original set you provided will have been added to the end of the file. If you look at the produced files, you should see that the spacers you provided and their fasta headers have been used by CCTK.
