# CRISPR_scripts
Scripts to identify and store information about CRISPR arrays in reads or assemblies.

## What is each script useful for?

The choice of which of these scripts to use depends on 2 things:
  1. Are you starting with reads or assemblies? 
  2. Do you know do you know the repeat sequence of the CRISPR type you are looking for?

If you have assemblies, but you don't know which CRISPR types you are looking for, you can use either minced2arrays.py or run_CRISPR_cas_finder_MP.py (followed by Process_CRISPR_cas_finder_out.py).
If you have assemblies and you know the repeat sequence(s) of the CRISPR type(s) you are looking for then you can use reps2spacers.py.
If you have reads then you will also need to know the repeat sequence(s) of the CRISPR type(s) you are looking for and should use narbl_pipeline.py and the associated scripts.



## minced2arrays.py
### Dependencies
[minced](https://github.com/ctSkennerton/minced)


## run_CRISPR_cas_finder_MP.py
### Dependencies
[CRISPRCasFinder](https://crisprcas.i2bc.paris-saclay.fr/Home/Download)

[perl](https://www.perl.org/)

## reps2spacers.py
### Dependencies
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with blast 2.9.0)



## narbl_pipeline.py
### Dependencies
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested with blast 2.9.0)

[fuzznuc](http://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html)

[FASTX](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
