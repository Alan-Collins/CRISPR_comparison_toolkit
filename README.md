[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Conda](https://anaconda.org/bioconda/cctk/badges/installer/conda.svg)](https://anaconda.org/bioconda/cctk)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cctk/badges/downloads.svg)](https://anaconda.org/bioconda/cctk)
[![Documentation Status](https://readthedocs.org/projects/crispr-comparison-toolkit/badge/?version=latest)](https://crispr-comparison-toolkit.readthedocs.io/en/latest/?badge=latest)


# CRISPR Comparison ToolKit (CCTK)

Tools to identify and compare CRISPR arrays.

Manuscript now available on bioRxiv (https://doi.org/10.1101/2022.07.31.502198)

## What is it?

CCTK is a collection of tools written in Python3 that are focused on the comparison of CRISPR arrays that share spacers with one another. CCTK includes two main tools: CRISPRdiff and CRISPRtree which can be used to visualize and analyze the relationships between CRISPR arrays. In addition, CCTK includes two scripts to identify CRISPR arrays genome assemblies and a script to generate a network representation of array relationships. The flow of data between CCTK tools is shown below.

![Flow of data between CCTK tools](https://github.com/Alan-Collins/CRISPR_comparison_toolkit/blob/main/images/CCTK_flowchart.png)

CCTK tools:
1. `cctk blast` and `cctk minced` Identify CRISPR arrays in assemblies using a method based on [MinCED](https://github.com/ctSkennerton/minced) or [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
2. `cctk network` Represent spacers shared between arrays as a network.
3. `cctk crisprdiff` Visualize alignment of homologous arrays.
4. `cctk crisprtree` Infer maximum parsimony tree explaining array relationships.
5. `cctk constrain` Assess whether CRISPR array relationships support phylogenetic relationships inferred based on other data (e.g. MLST)
6. `cctk spacerblast` Identify targets of CRISPR spacers in a BLAST database: assess whole length of imperfect spacer-protospacer matches and check for presence of protospacer adjacent motif (PAM).

## Documentation and tutorials

Documentation can be found at [Read the docs](https://crispr-comparison-toolkit.readthedocs.io/en/latest/)

## Installation

CCTK is a collection of scripts that designed to be run from the command line. Installation is as simple as downloading this repository to your computer and running the scripts in the terminal.

### Conda - recommended

CCTK has a number of dependencies to function. These are listed in the [Dependencies](#dependencies) section below.

Conda can be used to install CCTK and all of its dependencies easily.

The latest version of CCTK can be installed from Anaconda using:

`conda install -c bioconda cctk`

### Git clone

To download this repository using git:

`git clone https://github.com/Alan-Collins/CRISPR_comparison_toolkit.git`

## Dependencies

N.B. some dependencies are only required for certain tools. If you are only using CRISPRtree and CRISPRdiff then you only need python3 with matplotlib and dendropy.

- [python3](https://www.python.org/downloads/) (Tested with Python >= 3.8)
  - [Matplotlib](https://matplotlib.org/3.1.1/users/installing.html)
  - [DendroPy](https://dendropy.org/downloading.html)
  - [NumPy](https://numpy.org/install/)

- [MinCED](https://github.com/ctSkennerton/minced)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with BLAST 2.9.0)




