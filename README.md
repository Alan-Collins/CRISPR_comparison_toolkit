# CRISPR_tree

Generate tree representations of the relationships between CRISPR arrays.

## CRISPR_nj.py

Use the neighbour joining algorithm to generate a tree based on a pairwise distance matrix of distances between CRISPR arrays.

Uses a Needleman-Wunsch alignment algorithm to align pairs of CRISPR arrays based on their spacer content and then uses hamming distance to count the number of differences between those aligned sequences.

## Example output

Given the below cartooned alignments of CRISPR arrays, where:
* Position on the x-axis is the position in each CRISPR array relative to the leader end
* Colour indicates identical spacers.
* Lines drawn between identical colours in arrays plotted adjacent help to highlight identical spacers
* Thin black line in arrays indicate spacers unique to that array

![Alignment of CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/images/alignment.png)

The CRISPR_nj.py script produces the following tree in newick format as well as an optional graphical representation produced by ete3:

```
(((A:7.0, B:1.0):2.0, C:3.0), D);

         /-A
      /-|
   /-|   \-B
  |  |
--|   \-C
  |
   \-D
```

N.B. While the representation produced by ete3 appears rooted, the trees produced by this approach are unrooted. The newick representation of this tree can  be given to software such as iToL to examine it in an unrooted format such as that in the below image:

![Tree representation of relationships between CRISPR arrays](https://github.com/Alan-Collins/CRISPR_tree/images/tree.png)
