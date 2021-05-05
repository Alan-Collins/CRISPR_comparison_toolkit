#!/usr/bin/env python3


# AUTHOR        :    ALAN COLLINS
# VERSION       :    v2
# DATE          :    2021-5-5
# DESCRIPTION   :    Given fasta CRISPR repeats and a blast db of genomes, pulls out spacer arrays of >= 2 spacers
# CHANGELOG
# V1.1 
#   Added contig information to array location column in .csv output file. 
#   Contig information is taken from the sseqid of the blast result.
# V1.1
#   Removed print commands that I had forgotten to get rid of.
# v1.2
#   Fixed mistake where the wrong array_class attribute was being used for storing array ends points for 3'-5' arrays.
#   added trailing newlines to the end of output files.
# v1.3
#   Fixed issue where the end of an array was not found properly if there was only one array in one genome.
#   Removed btop from blast as it was unused
#   Added stderr catching and error printing for blastn and blastdbcmd run in case it fails.
# v1.3.1
#   Fixed issue handling large blastdbs that are split over multiple .nhr files.
# v2
#   Added e-value cutoff and max_target_seqs options to the blastn step so that blast doesn't stop looking for repeats prematurely with large databases.
#   Added argparse options to control thread count, evalue, and max_target_seqs blast options
#   Added multithreading to the spacer finding function to the blastdbcmd step that finds spacer sequences as this is the slowest part.
#   Added user provided option to control type of regex used for grep against database files.

import sys
import os
import argparse
import subprocess
import multiprocessing
import re
from collections import defaultdict



class blast_result():

    def __init__(self, blast_line):
        bits = blast_line.split('\t')
        self.qseqid = bits[0]

        self.sseqid = bits[1]
        self.pident = float(bits[2])
        self.length = int(bits[3])
        self.mismatch = int(bits[4])
        self.gapopen = bits[5]
        self.qstart = int(bits[6])
        self.qend = int(bits[7])
        self.sstart = int(bits[8])
        self.send = int(bits[9])
        self.evalue = bits[10]
        self.bitscore = bits[11]
        self.qlen = int(bits[12])
        self.slen = bits[13]
        self.strand = 'plus' if self.sstart < self.send else 'minus'
        self.trunc = False # Store whether repeat might be truncated (i.e. blast result wasn't full qlen)
        self.sstart_mod = False # If repeat appears tuncated, store revised genome locations here
        self.send_mod = False

class blast_entry():
    def __init__(self, blast_result, args):
        self.qseqid = blast_result.qseqid
        self.sseqid = blast_result.sseqid
        #self.regex_pattern = args.regex_pattern
        #self.blastdb = args.blast_db_path
        if blast_result.sstart < blast_result.send:
            self.strand = 'plus'
            if blast_result.trunc:
                self.sstart = blast_result.sstart_mod
                self.send = blast_result.send_mod
            else:
                self.sstart = blast_result.sstart
                self.send = blast_result.send
        else:
            self.strand = 'minus'
            if blast_result.trunc:
                self.sstart = blast_result.send_mod
                self.send = blast_result.sstart_mod
            else:
                self.sstart = blast_result.send
                self.send = blast_result.sstart

class array_class():
    def __init__(self):
        self.genome = ''
        self.contig = ''
        self.start = 0
        self.stop = 0
        self.repeats = []
        self.repeat_id = ''
        self.spacers = []
        self.spacer_ids = []
        self.array_id = 0

def hamming(string1, string2):
    dist = 0
    for i in range(max(len(string1), len(string2))):
        if i < len(string1) and i < len(string2):
            if string1[i] != string2[i]:
                dist += 1
        else:
            dist += 1

    return dist

def rev_comp(string):
    rev_str = ''
    rev_comp_lookup = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c"}
    for i in reversed(string):
        if i in "ATCGatcg":
            rev_str += rev_comp_lookup[i]
        else:
            rev_str += i
    return rev_str


def run_blastcmd(db, seqid, start, stop, strand):
    x = subprocess.run("blastdbcmd -db {} -entry {} -range {}-{} -strand {}".format(db, seqid, start, stop, strand), shell=True, universal_newlines = True, capture_output=True) 
    if x.stderr:
        print("ERROR running blastdbcmd on {}:\n{}".format(db, x.stderr))
        sys.exit()
    else:
        return x.stdout.split('\n')[1]


def pool_MP_spacer_finder(array_entries, threads, chunksize):

    pool = multiprocessing.Pool(processes=int(threads))

    pool.imap_unordered(build_arrays_MP, array_entries, chunksize)
    pool.close()
    pool.join()


def build_arrays_MP(array_entry):
    try:
        array = array_class()
        array.genome = re.match(args.regex_pattern, array_entry[0].sseqid)[0]
        n_reps = len(array_entry)
        array.repeat_id = array_entry[0].qseqid
        array.contig = array_entry[0].sseqid
        if array_entry[0].strand == 'plus':
            array.start = array_entry[-1].send
            array.stop = array_entry[0].sstart
            for i, entry in enumerate(array_entry):

                rep = run_blastcmd(args.blast_db_path, entry.sseqid, entry.sstart, entry.send, entry.strand)
                array.repeats.append(rep)
                if i+1 != n_reps:
                    spacer = run_blastcmd(args.blast_db_path, entry.sseqid, entry.send+1, array_entry[i+1].sstart-1, entry.strand)
                    array.spacers.append(spacer)
        
        else:
            array.start = array_entry[0].send
            array.stop = array_entry[-1].sstart
            for i, entry in enumerate(reversed(array_entry)):
                rep = run_blastcmd(args.blast_db_path, entry.sseqid, entry.sstart, entry.send, entry.strand)
                array.repeats.append(rep)
                if i+1 != n_reps:
                    spacer = run_blastcmd(args.blast_db_path, entry.sseqid, array_entry[-(i+2)].send+1, entry.sstart-1, entry.strand)
                    array.spacers.append(spacer)
        return array


    except Exception as e:
        print(e)
        sys.exit()


def blastn_to_arrays(args):
    
    blastn_command = "blastn -query {} -db {} -task blastn-short -outfmt '6 std qlen slen' -num_threads {} -max_target_seqs {} -evalue {} {}".format(args.repeats_file, args.blast_db_path, args.num_threads, args.max_target_seqs, args.evalue, args.other_blast_options)
    blast_run = subprocess.run(blastn_command, shell=True, universal_newlines = True, capture_output=True)
    # blast_lines = [blast_result(i) for i in subprocess.run(blastn_command, shell=True, universal_newlines = True, capture_output=True).stdout.split('\n') if len(i) > 0]
    if blast_run.stderr:
        print("ERROR running blast on {}:\n{}".format(args.blast_db_path, blast_run.stderr))
        sys.exit()
    blast_lines = [blast_result(i) for i in blast_run.stdout.split('\n') if len(i) > 0]
    

    good_hits = []

    for line in blast_lines:
        if line.length > 0.9 * line.qlen: # Keep any blast hits where the match length is over 90% of query
            if line.pident > 90:
                if line.length != line.qlen: # If match length isn't 100%, extend repeat to expected length
                    line.trunc = True
                    if line.sstart < line.send: # if the blast hit is on the positive strand
                        if line.qstart != 1: # if the 5' end of repeat is missing
                            line.sstart_mod = line.sstart - (line.qstart - 1) # adjust sstart by the number of bases missing from repeat 5' end
                        if line.qend != line.length: # if the 3' end of the repeat is missing
                            line.send_mod = line.send + (line.qend - line.length) # adjust send by the number of bases missing from repeat 3' end
                    elif line.sstart > line.send: # if the blast hit is on the negative strand
                        if line.qstart != 1: # if the 5' end of repeat is missing
                            line.send_mod = line.send - (line.qstart - 1) # adjust sstart by the number of bases missing from repeat 5' end
                        if line.qend != line.length: # if the 3' end of the repeat is missing
                            line.sstart_mod = line.sstart + (line.qend - line.length) # adjust send by the number of bases missing from repeat 3' end
                good_hits.append(blast_entry(line, args))
    good_hits_sorted = sorted(good_hits, key = lambda x: (x.sseqid, x.strand, x.sstart))
    array_entries = identify_same_array_hits(good_hits_sorted, args)
    # print(array_entries)
    

    if len(array_entries) > 0:
        arrays = pool_MP_spacer_finder(array_entries, args.num_threads, int(len(array_entries)//int(args.num_threads)))

    return arrays


def identify_same_array_hits(blast_entries, args): 
    # Using position and orientation of sorted (by contig and start position) repeats, 
    # group them into arrays if there are at least 3 repeats in a group. 
    # Otherwise discard them
    arrays = []
    n_spacers = 0 # count spacers in the current array
    last_contig = False
    last_loc = False # store start coord of last repeat
    last_strand = False # store orientation of last repeat
    for entry in blast_entries:
        if last_loc:
            if last_loc > entry.sstart - 150 and last_strand == entry.strand and last_contig == entry.sseqid:
                this_array.append(entry)
                last_loc = entry.sstart
                if entry == blast_entries[-1]:
                    if len(this_array) > 2:
                        arrays.append(this_array)
                        last_loc = entry.sstart
                        last_strand = entry.strand
                        last_contig = entry.sseqid
                        this_array = [entry]

                    else:
                        last_loc = entry.sstart
                        last_strand = entry.strand
                        last_contig = entry.sseqid
                        this_array = [entry]
                
            else:
                if len(this_array) > 2:
                    arrays.append(this_array)
                    last_loc = entry.sstart
                    last_strand = entry.strand
                    last_contig = entry.sseqid
                    this_array = [entry]

                else:
                    last_loc = entry.sstart
                    last_strand = entry.strand
                    last_contig = entry.sseqid
                    this_array = [entry]
        else:
            last_loc = entry.sstart
            last_strand = entry.strand
            last_contig = entry.sseqid
            this_array = [entry]

    return arrays

parser = argparse.ArgumentParser(
    description="Given fasta CRISPR repeats and a blast db of genomes, pulls out spacer arrays of >= 2 spacers")
parser.add_argument(
    "-r", dest="repeats_file", required = True,
    help=""
    )
parser.add_argument(
    "-d", dest="blast_db_path", required = True,
    help="path to blast db files (not including the file extensions)"
    )
parser.add_argument(
    "-o", dest="outdir", required = True,
    help="path to directory you want output files saved in"
    )
parser.add_argument(
    "-e", dest="evalue", required = False, default='1',
    help="set the evalue cutoff below which blastn will keep blast hits when looking for CRISPR repeats in your blast database. Useful for reducing inclusion of low quality blast hits with big databases in combination with the -m option."
    )
parser.add_argument(
    "-m", dest="max_target_seqs", required = False, default='10000',
    help="Set the max_target_seqs option for blastn when looking for CRISPR repeats in your blast database. Blast stops looking for hits after finding and internal limit (N_i) sequences for each query sequence, where N_i=2*N+50. These are just the first N_i sequences with better evalue scores than the cutoff, not the best N_i hits. Because of the nature of the blast used here (small number of queries with many expected hits) it may be necessary to increase the max_target_seqs value to avoid blast ceasing to search for repeats before all have been found. The blast default value is 500. The default used here is 10,000. You may want to reduce it to increase speed or increase it to make sure every repeat is being found. If increasing this value (e.g. doubling it) finds no new spacers then you can be confident that this is not an issue with your dataset."
    )
parser.add_argument(
    "-t", dest="num_threads", required = False, default='1',
    help="Number of threads you want to use for the blastn step of this script."
    )
parser.add_argument(
    "-p", dest="regex_pattern", required = True,
    help="In order to identify which genome a spacer was found in, put the genome id in your genome files fasta headers before making the blast db and then provide a regex pattern that can extract that id from the fasta header here. e.g. for the fasta header: >animaloris_GCA_900637855.1_LR134440.1 the genome id (animaloris_GCA_900637855.1) can be extracted with [A-Za-z]+_GCA_[0-9]+\.[0-9]"
    )
parser.add_argument(
    "-q", dest="regex_type", required = True, default='P',
    help="Regex type declaration option for grep (e.g. E, P, etc). Your regex pattern will be used with grep to find genome names in the .nhr file of your blast database. If you want to use regex patterns like \d+ you need to use the P option for Perl regex. Test your regex using the function 'grep -haoP' or 'grep -haoE' etc with your regex pattern against the nhr file in your blastdb."
    )
parser.add_argument(
    "-x", dest="other_blast_options", required = False, default='',
    help="If you want to include any other options to control the blastn command, you can add them here. Options you should not provide here are: blastn -query -db -task -outfmt -num_threads -max_target_seqs -evalue"
    )


args = parser.parse_args(sys.argv[1:])

outdir = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir

# Find the names of all the genomes being searched by looking for the user-provided regex in the blast-db .nhr files.
# for all genomes, start their entry in the genome_CRISPR_dict with the placeholder False. This indicated no CRISPRs found.
# If a CRISPR array is found later this entry will be overwritten with infor about the arrays.
genome_CRISPR_dict = { k : ['False'] for k in subprocess.run("blastdbcmd -db {} -entry all -outfmt '\%a' | grep -o{} '{}' | sort | uniq".format(args.blast_db_path, args.regex_type, args.regex_pattern), shell=True, universal_newlines = True, capture_output=True).stdout.split('\n') if len(k) > 0} 

all_arrays = blastn_to_arrays(args)


spacer_dict = {}
spacer_count_dict = defaultdict(int)
array_dict = {}
array_count = 0

for i, array in enumerate(all_arrays):
    try:
        print(array.spacers)
        sys.exit()
    except:
        pass
#     for spacer in array.spacers:
#         if spacer in spacer_dict.keys():
#             all_arrays[i].spacer_ids.append(spacer_dict[spacer])
#         else:
#             spacer_count_dict[array.repeat_id] +=1
#             spacer_dict[spacer] = 'rep_' + array.repeat_id + '_sp_' + str(spacer_count_dict[array.repeat_id])
#             all_arrays[i].spacer_ids.append(spacer_dict[spacer])
            
#     if tuple(array.spacer_ids) in array_dict.keys():
#         all_arrays[i].array_id = array_dict[tuple(array.spacer_ids)]
#     else:
#         array_count += 1
#         all_arrays[i].array_id = array_count
#         array_dict[tuple(array.spacer_ids)] = array_count

# genome_arrays = defaultdict(list) # list of arrays present in each genome
# array_genomes = defaultdict(list) # list of genomes with this array
# spacer_genomes = defaultdict(list) # list of genomes with this spacers

# for array in all_arrays:
#     genome_arrays[array.genome].append(array)
#     array_genomes[array.array_id].append(array.genome)
#     for spacer in array.spacer_ids:
#         spacer_genomes[spacer].append(array.genome)


# for genome, arrays in genome_arrays.items():
#     n_arrays = str(len(arrays))
#     spacers_list = "\t".join(["{}: {}".format(i+1, " ".join(array.spacers)) for i, array in enumerate(arrays)])
#     spacer_id_list = "\t".join(["{}: {}".format(i+1, " ".join(array.spacer_ids)) for i, array in enumerate(arrays)])
#     array_id_list = "\t".join(["{}: {}".format(i+1, array.array_id) for i, array in enumerate(arrays)])
#     array_locs = "\t".join(["{}: {} {} {}".format(i+1, array.contig, array.start, array.stop) for i, array in enumerate(arrays)])
#     genome_CRISPR_dict[genome] = ['True', n_arrays, spacers_list, spacer_id_list, array_id_list, array_locs]

# outcontents = ["Genome,Has_CRISPR,Array_count,Spacers,Spacer_IDs,Array_IDs,Array_locations"]

# for k,v in genome_CRISPR_dict.items():
#     outcontents.append(",".join([k]+v))

# with open(outdir + "CRISPR_summary_table.csv", 'w') as fout:
#     fout.write('\n'.join(outcontents)+ '\n')

# outcontents = ["Array_ID\tgenomes_with_array"]

# for k,v in array_genomes.items():
#     outcontents.append('{}\t{}'.format(k, ' '.join(v)))

# with open(outdir + "Array_representatives.txt", 'w') as fout:
#     fout.write('\n'.join(outcontents) + '\n')

# outcontents = []

# for k,v in spacer_dict.items():
#     outcontents.append(">{}\n{}".format(v,k))

# with open(outdir + "CRISPR_spacers.fna", 'w') as fout:
#     fout.write('\n'.join(outcontents) + '\n')

# outcontents = ["Spacer_ID\tGenomes_with_spacer"]

# for k,v in spacer_genomes.items():
#     v = list(set(v))
#     outcontents.append("{}\t{}".format(k, " ".join(v)))

# with open(outdir + "Spacer_representatives.txt", 'w') as fout:
#     fout.write('\n'.join(outcontents) + '\n')
