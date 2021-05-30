# CRISPR_scripts

There are two groups of scripts in this repository. Scripts to find CRISPRs in assemblies or reads, and scripts to compare the CRISPRs found in those assemblies or reads. The CRISPR finding scripts are found in the [Find_CRISPR_arrays](https://github.com/Alan-Collins/CRISPR_scripts/tree/main/Find_CRISPR_arrays) directory. The CRISPR analysis scripts are in the [Visualizing_array_relatedness](https://github.com/Alan-Collins/CRISPR_scripts/tree/main/Visualizing_array_relatedness) directory.

## Find_CRISPR_arrays directory

Scripts to identify and store information about CRISPR arrays in reads or assemblies.

### What is each script useful for?

The choice of which of these scripts to use depends on 2 things:
  1. Are you starting with reads or assemblies? 
  2. Do you know do you know the repeat sequence of the CRISPR type you are looking for?

If you have assemblies, but you don't know which CRISPR types you are looking for, you can use either minced2arrays.py or run_CRISPR_cas_finder_MP.py (followed by Process_CRISPR_cas_finder_out.py).
If you have assemblies and you know the repeat sequence(s) of the CRISPR type(s) you are looking for then you can use reps2spacers.py.
If you have reads then you will also need to know the repeat sequence(s) of the CRISPR type(s) you are looking for and should use narbl_pipeline.py and the associated scripts.



### minced2arrays.py
#### Dependencies
[minced](https://github.com/ctSkennerton/minced)

#### Description



### run_CRISPR_cas_finder_MP.py
#### Dependencies
[CRISPRCasFinder](https://crisprcas.i2bc.paris-saclay.fr/Home/Download)

[perl](https://www.perl.org/)

### reps2spacers.py
#### Dependencies
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with blast 2.9.0)



### narbl_pipeline.py
#### Dependencies
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with blast 2.9.0)

[fuzznuc](http://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html)

[FASTX](http://hannonlab.cshl.edu/fastx_toolkit/download.html)

## Visualizing_array_relatedness directory

Scripts to compare and visualize the different CRISPR arrays you have.

### What is each script useful for?

### CRISPRspacers2network.py
#### Dependencies
None

#### Description

In order to understand how the CRISPR arrays in your dataset are related to one another, you can assess similarity by measuring how many spacers a given pair of arrays share. This can be acheived using the script CRISPRspacers2network.py. It outputs a network representation of your arrays as well as a file detailing which arrays are found in which of your assemblies. The network file can be visualized with software such as [cytoscape](https://cytoscape.org/). See image below for an example of a network representation of CRISPR arrays.

![Example array network in cytoscape](https://github.com/Alan-Collins/CRISPR_scripts/blob/main/Example_images/array_network_viridis.png)

### CRISPR_alignment_plot.py
#### Dependencies
[matplotlib](https://matplotlib.org/3.1.1/users/installing.html)

#### Description

If you want to understand more about how the arrays in a cluster of your network are related to one another you can use one of the two CRISPR alignment plot scripts to visualize the CRISPR arrays with identical spacers colour-coded to highlight similarities. 

CRISPR_alignment_plot.py can be run on a user defined selection of array IDs and will plot a cartoon representation of those arrays in which shared spacers are coloured the same and spacers that are unique to an array are represented with a thin black line. In order to more clearly highlight identical spacers shared between arrays, a line is drawn connecting identical spacers between adjacently plotted arrays. See example below.

The order of arrays in the plot can be either provided by the user with the -preordered option, an approximate order can be provided with the -approxordered option, or the best order will be found by a search algorithm to try to put the most related arrays next to one another in the plot (default behaviour).

![Example array network in cytoscape](https://github.com/Alan-Collins/CRISPR_scripts/blob/main/Example_images/related_array_cluster.png)


### All_CRISPR_alignment_plot.py
#### Dependencies
[matplotlib](https://matplotlib.org/3.1.1/users/installing.html)

#### Description

If you just want to plot alignments of every array in your network that shares spacers then this is the script for you! It builds the network of related arrays and plots all the clusters of related arrays (connected by at least a user-provided threshold of shared spacers; default = 2) using the CRISPR_alignment_plot.py script in this repository.


