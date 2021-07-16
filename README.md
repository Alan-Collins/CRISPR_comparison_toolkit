# CRISPR_tree

Generate tree representations of the relationships between CRISPR arrays. Two scripts are included which use different approaches to construct a tree.

CRISPR_nj.py uses a basic neighbour-joining algorithm. This works well for simple array relationships where the only differences are new leader-end spacers. However, it produces negative branch lengths and unintuitive trees when given more complicated data.

CRISPR_mp.py uses a branch-and-bound maximum parsimony approach to find an approximately best tree. If the number of leaves are small then all trees can feasibly be checked and the best tree can be found. Otherwise a number of trees to test can be set by the user and the best tree or trees found during that search will be returned. As part of the maximum parsimony search algorithm, ancestral states are inferred that are intermediate between each input array. The accuracy of these intermediate states is important when assessing the quality of the final tree. This script can produce a png output (using `-o`) showing the topology of an arbitrarily rooted tree in which all leaf and internal CRISPR arrays are cartooned next to their node labels. This output can be used to manually assess how you feel about the parismony of the output tree.

## Input file format

N.B. All arrays must be properly oriented with leader-end spacers at the start of the list in order for the parsimony model to function properly. If you do not know the correct orientation of your arrays then the output tree may be meaningless.

Both CRISPR_nj.py and CRISPR_mp.py take the same input: A tab-delimited text file description of Array IDs followed by a column that is not used and can contain descriptive information if you like, followed by a tab-delimited list of either spacer IDs or sequences. The arrays are compared by string comparisons of each spacer so it doesn't matter what the spacers are represented by as long as the representation is consistent between arrays. An example of an input file is:

```
1231    spacer_IDs: 2698 274 3191 18 3300 3280 7180 2034 6045 1737
1245    spacer_IDs: 961 1677 4809 694 1073 5324 6045 1737
1401    spacer_IDs: 961 1677 4809 694 1073 5324 4913 2014 2034 6045 1737 3917
1902    spacer_IDs: 6569 130 117 830 42 1015 1677 4809 2034 6045 1737 3917
```

The above example inputs were used to generate all below example figures. The arrays to be aligned are specified in the command as a space-separated list at the end of the command.

The relationship of the above example arrays are illustrated in the below image, where:
* Position on the x-axis is the position in each CRISPR array relative to the leader end
* Colour indicates identical spacers.
* Lines drawn between identical colours in arrays plotted adjacent help to highlight identical spacers
* Thin black line in arrays indicate spacers unique to that array


![Alignment of CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/alignment.png)

## CRISPR_nj.py

Use the neighbour joining algorithm to generate a tree based on a pairwise distance matrix of distances between CRISPR arrays.

Uses a Needleman-Wunsch alignment algorithm to align pairs of CRISPR arrays based on their spacer content and then uses one of two distance measures to score pairwise differences. These measures can be selected using the `-d` option.

* `hamming`: Use hamming distance to count the number of positions in the aligned sequences that differ.
* `evol`: Use a simple evolutionary model to count differences in the form of events. This model takes into account some simple elements of CRISPR evolution. E.g. each leader-end difference is an event (acquisition) while a gap in the middle is a single event (indel) no matter how big the gap.

## Example output

Using `python3 CRISPR_nj.py -a array_reps.txt -p -d evol 1902 1245 1401 1231`, the CRISPR_nj.py script produces the following unrooted tree in newick format as well as an optional graphical representation produced by dendopy:

```
(((1902:8.0, 1231:7.0):7.0, 1245:1.0):0, 1401:-3.5);

        /--- 1902
    /---+
/---+   \--- 1231
|   |
+   \------- 1245
|
\----------- 1401
```

N.B. While the representation produced by dendropy appears rooted, the trees produced by this approach are unrooted. The newick representation of this tree can be given to software such as iToL to examine it in an unrooted format such as that in the below image:

![NJ Tree representation of relationships between CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/nj_tree.png)

---

## CRISPR_mp.py

Use a branch-and-bound process to explore possible tree topologies and keep the most parsimonious tree or trees found during the search. This process is carried out using the following steps:

* The order of the list arrays is randomized. (this can be overruled using the `-x` option which allows the user to define the addition order to produce the same tree every time.)
* The first two arrays in the list are selected.
* The arrays are aligned using a needleman-wunsch algorithm.
* The aligned arrays are compared and an intermediate state between the two arrays is inferred.
* The input arrays (referred to hereafter as "extant" arrays) are added to a tree object as leaves and the inferred array (referred to hereafter as the "ancestor" is set as the parent node of the extant arrays)
* The number of events required to turn the ancestor into each extant array is counted and set as the branch length for each leaf.
* The next array in the list is selected and compared to each array already in the tree (extant an ancestor). Each comparison is scored by the number of events that would be required to turn the existing array in the tree into the array being considered.
* The most similar array in the tree is identified and an ancestor inferred between it and the array to be added.
* The new array, and newly inferred ancestor are added to the tree by joining them to the most similar array already in the tree.
* This process is repeated until all arrays have been added to the tree.
* Every time a new array is added, the total length of all branches in the tree is calculated and compared to the total branch length in the best tree found yet. If the total branch length exceeds that of the best tree found so far then the tree is abandoned and the process restarts with a new array order.
* If the tree is finished and has a lower total branch length than the best tree found before it, the new tree replaces the best tree and the best score becomes the total branch length of the newly identified best tree.
* If a tree is finished and has a score equal to the best tree found before it, both are kept and will be printed as output for the user after the script finishes running.

## Example output

Using `python3 CRISPR_mp.py -a array_reps.txt -o images/mp_tree_annotations.png -r 100 -p 1902 1245 1401 1231`, the CRISPR_mp.py  script produces the following unrooted tree in newick format as well as an optional graphical representation produced by dendopy:

```
             /------ 1401
       /-----Int_b
/------Int_a \------ 1245
|      |
Int_c  \------------ 1902
|
\------------------- 1231


(((1401:1,1245:1)'Int_b':2,1902:6)'Int_a':2,1231:8)'Int_c':0;
```

The produced newick string can then be given to tree viewing software such as iToL where it can be viewed in its real, unrooted form. The following image shows the tree produced by an unweighted parimony analysis (i.e. `-q`, `-i`, `-d` are not given and so default costs of 1 are used):

![Unrooted MP tree representation of relationships between CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/mp_tree.png)

In addition to viewing the produced tree, the user is able to assess the information that went into the choice of the output tree as the most parsimonious. Using the `-o` option and providing the path to an output .png file will result in the production of one image per best tree illustrating the inferred ancestral states that were used to calculate parsimony scores across the tree. In the event that there is more than one best tree, the last 4 characters of your given output file (e.g. ".png") will be replaced with "\_1.png" where 1 is an increasing count of the best tree the image refers to (i.e. tree number 1, tree number 2 etc in the order they appear in the terminal output.)

The tree image produced by the script is drawn in a rooted form with the same structure as the ascii tree printed as terminal output and the same representation as the rooted form the iToL will produce by default when given the output newick string. An example of the rooted form of the tree produced by iToL is shown below.

![Rooted MP tree representation of relationships between CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/mp_tree_rooted.png)

The cartoon representations of the CRISPR spacers in both the extant and inferred ancestral arrays are the same as those in the CRISPR array representations shown above. i.e. 
* Spacers are colour-coded so that identical spacers are coloured the same.
* Arrays are oriented with the leader end on the left and trailer end on the right.
* Thinner, black lines indicate spacers that are unique to the array in which they are drawn.

Inferred ancestral arrays are assigned node names prefixed with "Int_" to indicate that they are internal nodes (i.e. not leaves). Any node names "Int_..." was inferred by the CRISPR_mp.py script.

The below image shows the representation of the arrays used to infer the "best" tree for the example data

![Rooted MP tree with ancestral and extant array cartoons](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/mp_tree_annotations.png)

In addition to showing the cartooned CRISPR arrays of extant and hypothetical ancestor arrays, the events that were predicted to have occured between these states can also be shown. The below image shows the same tree with each event that has occured in an array since its hypothetical ancestor (the array of the parent node in the tree). 3 predicted event types are shown:

1. A large "X" indicates a region has been deleted.
2. A box around a run of spacers indicates these spacers have been acquired through extopic acquisition or recombination
3. A squiggly line under a run of spacers indicated these spacers have been added to the leader end of the array.

![Rooted MP tree with ancestral and extant array cartoons plus highlighted events](https://github.com/Alan-Collins/CRISPR_tree/blob/master/images/mp_tree_annotations_highlight.png)
