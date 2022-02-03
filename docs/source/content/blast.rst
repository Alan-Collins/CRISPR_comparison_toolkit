blast
=====

Before you run
--------------

``cctk blast`` requires you to provide the sequences you wish to search in the form of a blast database. This can be acheived simply using the ``makeblastdb`` from NCBI BLAST+ command which is included in the conda installation of CCTK. However, before making your blastdb, please confirm that your sequences meet the following requirements:

1. No pipe symbols ("|") in any of your fasta headers.
2. None of the fasta headers in the sequences are the same. 

Making a blastdb from a directory containing your sequences to search can be acheived in two steps (if you are only searching a single sequence file skip the first step). Assuming your sequences are in a directory called sequences/:

.. code-block:: shell

	cat sequences/* > all_seqs.fna
	makeblastdb -in all_seqs.fna -out seq_blastdb -dbtype nucl -parse_seqids

**N.B.** It is essential to include the ``-parse_seqids`` option when creating the blastdb. You can be check if your blastdb was made using ``-parse_seqids`` by looking for files ending in .nog and .nos. If those files exist associated with your blastdb then it was.

Introduction
------------

