#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-8-3
# DESCRIPTION :  Process BLAST output of spacers against a blastdb. For results that have cut off due to mismatches, extend the hit to the full length and report mismatches. Report up- and down-stream bases for PAM analysis.

import sys
import argparse

parser = argparse.ArgumentParser(
	description="Process BLAST output of spacers against a blastdb. For results that have cut off due to mismatches, extend the hit to the full length and report mismatches. Report up- and down-stream bases for PAM analysis."
	)
parser.add_argument(
    "-d", dest="blast_db_path", required = True,
    help="path to blast db files (not including the file extensions). The blastdb must have been made with the option '-parse-seqids' for this script to function."
    )
parser.add_argument(
    "-s", dest="spacer_file", required = True,
    help="The file with your spacers in fasta format."
    )
parser.add_argument(
    "-o", dest="outfile", required = True,
    help="path to output file."
    )
parser.add_argument(
    "-e", dest="evalue", required = False, default='10',
    help="DEFAULT: 10. set the evalue cutoff below which blastn will keep blast hits when looking for CRISPR repeats in your blast database. Useful for reducing inclusion of low quality blast hits with big databases in combination with the -m option."
    )
parser.add_argument(
    "-m", dest="max_target_seqs", required = False, default='10000',
    help="DEFAULT: 10000. Set the max_target_seqs option for blastn when looking for CRISPR repeats in your blast database. Blast stops looking for hits after finding and internal limit (N_i) sequences for each query sequence, where N_i=2*N+50. These are just the first N_i sequences with better evalue scores than the cutoff, not the best N_i hits. Because of the nature of the blast used here (small number of queries with many expected hits) it may be necessary to increase the max_target_seqs value to avoid blast ceasing to search for repeats before all have been found. The blast default value is 500. The default used here is 10,000. You may want to reduce it to increase speed or increase it to make sure every repeat is being found. If increasing this value (e.g. doubling it) finds no new spacers then you can be confident that this is not an issue with your dataset."
    )
parser.add_argument(
    "-t", dest="num_threads", required = False, default=1, type=int,
    help="DEFAULT: 1. Number of threads you want to use for the blastn step of this script."
    )
parser.add_argument(
    "-x", dest="other_blast_options", required = False, default='',
    help="DEFAULT: none. If you want to include any other options to control the blastn command, you can add them here. Options you should not provide here are: blastn -query -db -task -outfmt -num_threads -max_target_seqs -evalue"
    )


args = parser.parse_args(sys.argv[1:])
