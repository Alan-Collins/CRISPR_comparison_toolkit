spacerblast
===========

.. _spacerblast-intro:

Introduction
------------

``cctk spacerblast`` is essentially a wrapper script for ``blastn`` that performs additional functions that are useful when searching for the targets of CRISPR spacers. The additional functions performed by ``cctk spacerblast`` are:

1. Extends sequence matches to cover the entire length of the query sequence. For short query sequences, a small number of mismatches can result in the BLAST algorithm not extending the match through the mismatches even if additional sequence identity would be found. Users can specify a percent identity requirement that takes into account the full sequence length.

2. (Optional) Retrieves flanking sequence either side (or both sides) of the match location in the subject sequence.

3. (Optional) Identifies the presence of user-defined sequence motifs in sequence flanking the match and sorts blast hits based on the presence or absence of the motif.

4. (Optional) Masks user-specified regions in the subject sequence and ignores blast hits in those regions.

5. (Optional) Specify spacer seed region to exclude hits with mismatches within that region.

Before you run
--------------

``cctk spacerblast`` requires you to provide the sequences you wish to search in the form of a blast database. This can be acheived simply using the ``makeblastdb`` command from NCBI BLAST+ which is included in the conda installation of CCTK. However, before making your blastdb, please confirm that your sequences meet the following requirements:

1. No pipe symbols ("|") in any of your fasta headers.
2. None of the fasta headers in the sequences are the same. 

Making a blastdb from a directory containing your sequences to search can be acheived in two steps (if you are only searching a single sequence file skip the first step). Assuming your sequences are in a directory called sequences/:

.. code-block:: shell

	cat sequences/* > all_seqs.fna
	makeblastdb -in all_seqs.fna -out seq_blastdb -dbtype nucl -parse_seqids

**N.B.** It is essential to include the ``-parse_seqids`` option when creating the blastdb. You can be check if your blastdb was made using ``-parse_seqids`` by looking for files ending in .nog and .nos. If those files exist associated with your blastdb then it was.

Basic Usage
-----------

Simple BLAST with match extension
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most basic usage requires only 2 inputs: the blastdb name (i.e. without file extensions) to search and the query sequences.

.. code-block:: shell

	cctk spacerblast -d <blastdb path> -s <query sequence file>

N.B. If you are on a computer with multiple threads, you can speed up ``cctk spacerblast`` by specifying a number of threads to use with the ``-t`` option.

Not specifying an output file results in the output being directed to the stdout for easy downstream filtering with command-line tools. The above command results in an output like the example below:

.. code-block:: shell

	Spacer_ID	Target_contig	Protospacer_start	Protospacer_end	Percent_identity	mismatches	protospacer_sequence	mismatch_locations	target_strand
	1F_1	Assembly1_contig1	66234	66265	100.0	0	GGTACGTGGTTTCGACCAACAGCACTGCCCAA	GGTACGTGGTTTCGACCAACAGCACTGCCCAA	minus
	1F_1	Assembly5_contig2	91689	91720	78.125	7	GGTACGTGGTTTCGACAAGCAACGTGGCCCAG	GGTACGTGGTTTCGAC.A.CA.C...GCCCA.	plus
	1F_1	Assembly6_contig1	15636	15667	59.375	13	TGGTGCACGTTTCGACCAACAGCCTAGCGCCC	.G......GTTTCGACCAACAGC...GC.C..	plus

Checking flanking sequence for PAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also retrieve just the hits that have an adjacent PAM sequence. For example, to check for the *Pseudomonas aeruginosa* type I-F PAM, CC, which occurs 5' (or upstream) of the protospacer:

.. code-block:: shell

	cctk spacerblast -d <blastdb path> -s <query sequence file> -P CC -l up

You can specify your PAM in this way using any characters in the `IUPAC nucleotide alphabet <https://www.bioinformatics.org/sms/iupac.html>`_.

The above command returns something like the following output:

.. code-block:: shell
	
	Your specified PAM is at least 2 bases, but you only requested 0 upstream bases. 2 bases will now be retrieved on the upstream side.
	Spacer_ID	Target_contig	Protospacer_start	Protospacer_end	Percent_identity	mismatches	protospacer_sequence	mismatch_locations	upstream_bases	target_strand
	1F_1	Assembly1_contig1	66234	66265	100.0	0	GGTACGTGGTTTCGACCAACAGCACTGCCCAA	GGTACGTGGTTTCGACCAACAGCACTGCCCAA	CC	minus
	1F_1	Assembly6_contig1	15636	15667	59.375	13	TGGTGCACGTTTCGACCAACAGCCTAGCGCCC	.G......GTTTCGACCAACAGC...GC.C..	CC	plus

There are two differences to note here between the first output and the output when specifying a PAM sequence:

1. The output includes an additional column "upstream_bases" which contains the 2 bases flanking the protospacer on the 5' side. If the PAM had been specified as being on the other side (i.e. ``-l down``) then this column would be "downstream_bases".

2. A warning message was written to stderr stating that our PAM is two bases long while we didn't ask for any flanking bases to be checked. This message is intended to make clear why a different number of bases are returned if you specifically request fewer bases than the provided PAM requires. ``cctk spacerblast`` will automatically determine the shortes length of sequence that your PAM could match and will return at least that much sequence.


Output files
------------

The default behaviour of ``cctk spacerblast`` is to direct ouputs to the stdout and information and error messages to stderr. However, two output files can be produced if requested by the user using the ``-o`` and ``-q`` options.

``-o`` Main output file / hits with PAMs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option directs any output that would have been sent to stdout to the specified file instead. You can name this file and specify its location by providing the path to a file (i.e. ``-o <path to file>``)

If no PAM information is provided then this output file contains all BLAST hits that meet the percent identity and evalue thresholds. If PAM information is provided, this file will contain just the hits that were found to have an adjacent PAM.

``-q`` Hits without PAMs
^^^^^^^^^^^^^^^^^^^^^^^^

This file will only be generated if PAM information is provided. You can name this file and specify its location by providing the path to a file (i.e. ``-q <path to file>``).

If PAM information is provided, this file will contain all BLAST hits that were not found to have an adjacent PAM. Only hits that exceed the percent identity and evalue thresholds will be stored.

Advanced Usage
--------------

Advanced usage of ``cctk spacerblast`` is not much more complicated than the basic usage described above. There are three cases in which a more complicated usage is required:

Control amount and location of flanking sequence retrieved
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The amount of sequence retrieved from each side of BLAST hits can be controlled using command line input with the optione ``-n``, ``-u``, and ``-w``. If also specifying a PAM, at least enough sequence to match the PAM will be retrieved. If you request less sequence than is required to match the provided PAM, the length of sequence retrieved will be adjusted and an informative message will be written to stderr informing you.

``-n`` can be used to retrieve the same length of flanking sequence on both sides of BLAST hits,

Specify a PAM using a regex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would prefer to define your PAM as a regex rather than using IUPAC nucleotide codes, you can do that using the ``-R`` option. Regex PAM definition is useful when the number of bases is flexible or if you prefer to specify e.g. A, T, or G with "[ATG]" rather than using the IUPAC "D".

Mask regions of sequences in you blastdb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to ignore hits in certain regions of your subject sequences you can maks regions by providing a `BED format <https://en.wikipedia.org/wiki/BED_(file_format)#Format>`_ file with the ``-r`` option. Only the first 3 columns of the .bed file will be read so all other columns are optional.

This can be useful when extracting spacers and searching for CRISPR targets in the same set of sequences. It will allow you to ignore hits against CRISPR arrays as each spacer will return a perfect match against its location in the genome in which it was found. Both `cctk blast <blast.html>`_ and `cctk minced <minced.html>`_ return a .bed file of CRISPR array locations that can be used for this purpose.

Specify seed region
^^^^^^^^^^^^^^^^^^^

CRISPR targeting can often occur even without 100% identity between the spacer and protospacer. However, some CRISPR types have a strict requirement for 100% sequence identity within a "seed region".

You can specify a region to be treated as a seed region using ``-E``. Only hits with 100% identity in this region will be returned.

The format to specify a region if of the format start_index:end_index and is a 0-base range. i.e., the range starts at the first number and ends 1 before the last number. For example, to specify a seed region covering the first 6 bases of the spacer you would use ``-E 0:6``. 0:6 specifies the bases at index 0, 1, 2, 3, 4, and 5.

Note that if your range starts at the first base of the spacer or ends at the last base of the spacer, you can omit that index. E.g., ``-E :6`` is equivalent to ``-E 0:6`` shown above. ``-E 6:`` specifies a range from the seventh base (index 6 with 0-base counting) to the end of the spacer

If you wish to specify a range relative to the end of the spacer. you can use negative indices (counting back from -1, which is the index of the last base). Note that if the end of your desired range includes the last base, you must omit the end index from your range. E.g., to specify a seed region of the last 6 bases of a spacer you can use ``-E -6:``. Using negative indices can be useful if your seed region is a fixed length at the end of the spacer, but spacers can differ in length.

The option to specify a seed region is only relevant when also specifying a percent identity threshold (``-p``) below 100%. Percent identity between the spacer and protospacer will still be calculated as a percent of the total spacer length. E.g., if you are looking for at least 90% identity matches for a 30 base spacer and specify that the seed region is the first 10 bases, up to 3 mismatches will be tolerated as long as they do not occur in the first 10 bases.
