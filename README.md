[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Conda](https://anaconda.org/alan-collins/cctk/badges/installer/conda.svg)](https://anaconda.org/Alan-Collins/cctk)

# CRISPR Comparison ToolKit (CCTK)

Tools to identify and compare CRISPR arrays.

**Currently under construction.** 

N.B. All included tools work as intended, but are being actively improved. Command line options, defaults, and project organization will change over the coming weeks. 

## What is it?

CCTK is a collection of tools written in Python3 that are focused on the comparison of CRISPR arrays that share spacers with one another. CCTK includes two main tools. CRISPRdiff and CRISPRtree which can be used to visualize and analyze the relationships between arrays. In addition, CCTK includes several scripts to identify CRISPR arrays in illumina sequencing reads or assemblies and a script to generate a network representation of array relationships.

## Installation

CCTK is a collection of scripts that designed to be run from the command line. Installation is as simple as downloading this repository to your computer and running the scripts in the terminal.

### Conda - recommended

CCTK has a number of dependencies to function. These are listed in the [Dependencies](#dependencies) section below.

Conda can be used to install CCTK and all of its dependencies easily.

The latest version of CCTK can be installed from Anaconda using:

`conda install -c alan-collins cctk`

### Git clone

To download this repository using git:

`git clone https://github.com/Alan-Collins/CRISPR_comparison_toolkit.git`

## Dependencies

N.B. some dependencies are only required for certain tools. If you are only using CRISPRtree and CRISPRdiff then you only need python3 with matplotlib and dendropy.

- [python3](https://www.python.org/downloads/) (Tested with Python >= 3.8)
  - [matplotlib](https://matplotlib.org/3.1.1/users/installing.html)
  - [dendropy](https://dendropy.org/downloading.html)

- [minced](https://github.com/ctSkennerton/minced)
- [CRISPRCasFinder](https://crisprcas.i2bc.paris-saclay.fr/Home/Download)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with BLAST 2.9.0)
- [fuzznuc](http://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html)
- [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/download.html)



