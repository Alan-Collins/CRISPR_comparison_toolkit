#!/usr/bin/env python3


import sys
import os
import argparse
import subprocess
import multiprocessing
import re
from collections import defaultdict

from . import sequence_operations, file_handling

description = """
usage: cctk blast [-h] -r <path> -d <path (no extension)> -o <path> \
[-a <path>] [-p <string>] [-q <string>] [-t <int>] [-i <int>] [-c <float>] \
[-s <int>] [--min-shared <int>] [--append] [-e <float>] [-m <int>] [-b <int>] \
[-x <string>]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -r, --repeats         CRISPR repeats in FASTA format
  -d, --blastdb         path to blast db (excluding file extension)
  -o, --outdir          directory to store output files

other inputs:
  -a, --assemblies      file containing a list of your assembly names
  -p, --regex-pattern   pattern describing your assembly names
  -q, --regex-type      {E, P} regex type describing assembly names. Default: P
  -t, --threads         number of threads to use. Default: 1
  -i, --repeat-interval maximum interval between repeats. Default: 80
  -c, --percent-id      minumum percent ID of repeat BLAST hits. Default: 80
  -s, --snp-thresh      number of SNPs to consider spacers the same. Default: 0
  --min-shared          minimum spacers shared to draw an edge in network
  --append              add new CRISPR data to a previous run

BLASTn settings:
  -e, --evalue          blastn evalue. Default: 10
  -m, --max-target-seqs max_target_seqs option for blastn. Default: 10000
  -b, --batch-size      batch size for blastdbcmd. Default: 1000
  -x, --blast-options   input additional blastn options. Forbidden options: \
blastn -query -db -task -outfmt -num_threads -max_target_seqs -evalue
"""

def pool_MP_spacer_finder(array_entries, args, threads, chunksize):
    """
    Manages the multiprocess worker pool and returns to arrays found by
    build_arrays_MP to whatever called it.
    Args:
        array_entries (list): A list of lists of blast_entry class
          instances.
        args (argparse class): All of the argparse options given by the
          user.
        threads (int): Number of threads to use for multiprocessing.
        Chunksize (int):    Approximate size of batches into which pool
          should split array_entries 
    
    Returns:
        (list) List of FoundArray instances representing distinct arrays
    """

    pool = multiprocessing.Pool(processes=threads)

    options = [(array_entry, args) for array_entry in array_entries]
    output = pool.starmap(build_arrays_MP, options, chunksize)
    output = [i for i in output]
    pool.close()
    pool.join()

    return output


def build_arrays_MP(array_entry, args):
    """
    Builds the fstring and batch_locations objects to give to
    run_blastdbcmd and then calls that and returns the output.
    Args:
        array_entry (list): list of blast_entry class instances.
        args (argparse class): All of the argparse options given by the
          user.
    
    Returns:
        (FoundArray) FoundArray instance representing the contiguous
          array identified.
    
    Raises:
        Exception: Raises an exception and prints the exception if one
          occurs.
    """
    try:
        # retrieve spacer seqs in batch for quicker blastdbcmd usage.
        sp_fstring = '' 
        sp_batch_locations = ''
        # same for repeats
        rep_fstring = '' 
        rep_batch_locations = ''
       
        array = file_handling.FoundArray()
        if args.regex_pattern:
            array.genome = re.match(
                args.regex_pattern, array_entry[0].sseqid)[0]
        elif args.assemblies:
            array.genome = which_substring_in_string(
                args.assemblies,
                array_entry[0].sseqid)
        else:
            array.genome = sseqid
        n_reps = len(array_entry)
        array.repeats = [a.sseq for a in array_entry] 

        ### Retrieve repeat sequence using blastdbcmd and fstring

        array.repeat_id = array_entry[0].qseqid
        array.contig = array_entry[0].sseqid
        if array_entry[0].strand == 'plus':
            array.start = array_entry[0].sstart
            array.stop = array_entry[-1].send
            for i, entry in enumerate(array_entry):
                if i+1 != n_reps:
                    sp_fstring += '%s %s %s\n'
                    sp_batch_locations += '{} {}-{} {} '.format(
                        entry.sseqid,
                        entry.send+1,
                        array_entry[i+1].sstart-1,
                        entry.strand
                        )
                rep_fstring += '%s %s %s\n'
                rep_batch_locations += '{} {}-{} {} '.format(
                    entry.sseqid,
                    entry.sstart,
                    entry.send,
                    entry.strand
                    )
        else:
            array.reverse = True
            array.start = array_entry[-1].sstart
            array.stop = array_entry[0].send
            rev_array_entry = [a for a in reversed(array_entry)]
            for i, entry in enumerate(rev_array_entry):
                if i+1 != n_reps:
                    sp_fstring += '%s %s %s\n'
                    sp_batch_locations += '{} {}-{} {} '.format(
                        entry.sseqid,
                        rev_array_entry[i+1].sstart+1,
                        entry.send-1,
                        entry.strand
                        )
                rep_fstring += '%s %s %s\n'
                rep_batch_locations += '{} {}-{} {} '.format(
                    entry.sseqid,
                    entry.send,
                    entry.sstart,
                    entry.strand
                    )
        _, array.spacers = sequence_operations.run_blastcmd(
            [], # only used for multithread running
            args.blast_db_path,
            sp_fstring,
            sp_batch_locations
            )
        _, array.repeats = sequence_operations.run_blastcmd(
            [],
            args.blast_db_path,
            rep_fstring,
            rep_batch_locations
            )

        return array


    except Exception as e:
        print(e)
        sys.exit()


def extend_hit(hit):
    if hit.sstart < hit.send:
        # if the blast hit is on the positive strand
        if hit.qstart != 1:
            # if the 5' end of repeat is missing
            # adjust sstart by the number of bases missing from
            # repeat 5' end
            hit.sstart = hit.sstart - (hit.qstart - 1)
        if hit.qlen != hit.qend:
            # if the 3' end of the repeat is missing
            # adjust send by the number of bases missing from
            # repeat 3' end
            hit.send = hit.send + (hit.qlen - hit.qend)
    elif hit.sstart > hit.send:
        # if the blast hit is on the negative strand
        if hit.qstart != 1:
            # if the 5' end of repeat is missing
            # adjust sstart by the number of bases missing from
            # repeat 5' end
            hit.sstart = hit.sstart + (hit.qstart - 1)
        if hit.qlen != hit.qend:
            # if the 3' end of the repeat is missing
            # adjust send by the number of bases missing from
            # repeat 3' end
            hit.send = hit.send - (hit.qlen - hit.qend)


def blastn_to_arrays(args):
    """
    Runs blastn of provided repeats against provided blastdb. Processes
    the blast output to find hits that are good enough to keep.
    Passes those good hits to identify_same_array_hits and stores the
    list of arrays returned to it. Passes that list of arrays to
    pool_MP_spacer_finder and stores the arrays returned to it that now
    have spacer information. If arrays were found, these are returned.
    If none are found an error is raised.
    Args:
        args (argparse class): All of the argparse options given by the user.

    
    Returns:
        (list): List of lists of FoundArray instances representing
          distinct arrays.
    
    Raises:
        No arrays found: If no arrays were found in the blast database
          then this raises an error stating that and exits.
    """ 
    blastn_command = "blastn -query {} -db {} -task blastn-short -outfmt \
    '6 std qlen slen sseq' -num_threads {} -max_target_seqs {} -evalue {} {}\
    ".format(
        args.repeats_file,
        args.blast_db_path,
        args.num_threads,
        args.max_target_seqs,
        args.evalue,
        args.other_blast_options
        )
    blast_run = subprocess.run(
        blastn_command,
        shell=True,
        universal_newlines=True,
        capture_output=True
        )

    if blast_run.stderr:
        print("ERROR running blast on {}:\n{}".format(
            args.blast_db_path, blast_run.stderr))
        sys.exit()
    blast_lines = [
        file_handling.BlastResult(
            i) for i in blast_run.stdout.split('\n') if len(i) > 0]

    queries_dict = file_handling.fasta_to_dict(args.repeats_file)
    # If blast result is < query length, extend it to check if it is
    # a good match and update pid by comparing to repeat sequence

    hits_to_update = []
    fstring = ''
    batch_locations = ''
    extended_hits = []
    n_args = 0

    for n, line in enumerate(blast_lines):
        if line.length == line.qlen:
            continue
        
        extend_hit(line)

        if ((int(line.sstart) < 0)
            or (int(line.send) < 0)
            or (int(line.sstart) > int(line.slen))
            or (int(line.send) > int(line.slen))
        ):
            # If extending the hit would leave the contig then discard
            line.pident = 0
            continue

        hits_to_update.append(n)
        n_args+=3
        if n_args > args.batch_size:
            x = len(extended_hits)
            extended_hits += sequence_operations.run_blastcmd(
                [], # Only used for multithread running
                args.blast_db_path,
                fstring,
                batch_locations
                )[1]
            fstring = ''
            batch_locations = ''
            n_args = 3

        fstring += '%s %s %s\n'
        if line.sstart < line.send:
            batch_locations += '{} {}-{} {} '.format(
                line.sseqid,
                line.sstart,
                line.send,
                line.strand
                )
        else:
            batch_locations += '{} {}-{} {} '.format(
                line.sseqid,
                line.send,
                line.sstart,
                line.strand
                )
    
    extended_hits += sequence_operations.run_blastcmd(
        [],
        args.blast_db_path,
        fstring,
        batch_locations
        )[1]

    for n, new_seq in zip(hits_to_update, extended_hits):
        hit = blast_lines[n]
        hit.length = hit.qlen
        hit.pident = sequence_operations.percent_id(
            new_seq,
            queries_dict[hit.qseqid]
            )


    good_hits = []

    for line in blast_lines:
        # Discard any blast hits where pident is < some percent
        if line.pident < args.percent_id:
            continue
        
        # # Adjust start and end for reverse orientation matches
        # if line.sstart > line.send:
        #     line.sstart, line.send = line.send, line.sstart
        good_hits.append(line)
    # check_blast_overlap(good_hits)
    good_hits_sorted = sorted(
        good_hits, key = lambda x: (x.sseqid, x.strand, x.sstart))
    array_entries = identify_same_array_hits(good_hits_sorted, args)
    

    if len(array_entries) > 0:
        chunksize = max(int(len(array_entries)//args.num_threads), 1)
        arrays = pool_MP_spacer_finder(
            array_entries,
            args,
            args.num_threads,
            chunksize)

        return arrays

    else:
        print("No arrays identified from the query file {} in the database {}\
            ".format(args.repeats_file, args.blast_db_path))
        sys.exit()


def identify_same_array_hits(blast_entries, args):
    """
    Given a list of blast results, identify which of them are from a
    contiguous spacer array by checking:
    Are they close enough?
    Are they on the same strand?
    Are they in the same genome, contig?
    If it finds 3 or more repeats close together on the same strand of
    the same contig, groups them into an array.
    Args:
        blast_entries (list):   List of blast entry classes sorted by
          contig, then start, then location.
        args (argparse class):  The arguments provided to argparse.
    
    Returns:
        (list) List of lists of blast_entry classes.
    """
    arrays = []
    n_spacers = 0 # count spacers in the current array
    last_contig = False
    last_loc = False # store start coord of last repeat
    last_strand = False # store orientation of last repeat
    for entry in blast_entries:
        if last_loc:
            if last_loc == entry.sstart:
                # If two parts of a repeat matched separately and the 
                # matches were extended, don't process both.
                continue
            if ((last_loc > entry.sstart - args.repeat_interval)
                and (last_strand == entry.strand)
                and (last_contig == entry.sseqid)):
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


def which_substring_in_string(substrings, string):
    present = []
    for s in substrings:
        if s in string:
            present.append(s)

    # If more than one present, longest will be best
    if len(present) > 1:
        return max(present, key=len)
    else:
        return present[0]


def check_repeat_similarity(repeats_file):
    blastn_command = "blastn -query {} -subject {} -task blastn-short -outfmt \
        '6 std qlen slen sseq'\
        ".format(
            repeats_file,
            repeats_file,
            )
    blast_run = subprocess.run(
        blastn_command,
        shell=True,
        universal_newlines=True,
        capture_output=True
        )

    blast_lines = [
        file_handling.BlastResult(
        i) for i in blast_run.stdout.split('\n') if len(i) > 0]
    similar_pairs = []
    for line in blast_lines:
        if (
            (line.length > 0.9*line.qlen or line.length > 0.9*line.slen)
            and (line.pident > 90)
            and (line.qseqid != line.sseqid)):
            
            if not (((line.qseqid, line.sseqid) in similar_pairs)
                or ((line.sseqid, line.qseqid) in similar_pairs)):
                similar_pairs.append((line.qseqid, line.sseqid))
    if len(similar_pairs) > 0:
        sys.stderr.write(
            "The following pairs of repeats in the provided file are "
            "over 90% identical across over 90% of their length. This "
            "program can only identify arrays associated with distinct "
            "repeats.\nSimilar pairs of arrays:\n{}\n".format(
                "\n".join(["\t".join(i) for i in similar_pairs])))
        sys.exit()


def build_parser(parser):
    req_options = parser.add_argument_group("Required arguments")
    req_options.add_argument(
        "-r", "--repeats",
        metavar=" ",
        dest="repeats_file",
        required=True,
        help="fasta format file of your CRISPR repeats"
        )
    req_options.add_argument(
        "-d", "--blastdb",
        metavar=" ",
        dest="blast_db_path",
        required=True,
        help="path to blast db files not including the file extensions"
        )
    req_options.add_argument(
        "-o", "--outdir",
        metavar=" ",
        required=True,
        help="path to directory you want output files saved in"
        )

    other_options = parser.add_argument_group("Other inputs")
    other_options.add_argument(
        "-a", "--assemblies",
        metavar=" ",
        required=False,
        help="file containing a list of your assembly names"
        )
    other_options.add_argument(
        "-p", "--regex-pattern",
        metavar=" ",
        dest="regex_pattern",
        required=False,
        help="pattern describing your assembly names"
        )
    other_options.add_argument(
        "-q", "--regex-type",
        metavar=" ",
        required=False,
        default='P',
        help="DEFAULT: P. {E, P} Regex type to be used in step that calls "
            "grep. See CCTK documentation for details."
        )
    other_options.add_argument(
        "-t", "--threads",
        metavar=" ",
        dest="num_threads",
        required=False,
        default=1,
        type=int,
        help="DEFAULT: 1. Number of threads to use"
        )
    other_options.add_argument(
        "-i", "--repeat-interval",
        metavar=" ",
        required=False,
        default=80,
        type=int,
        help="DEFAULT: 80. Set the expected interval between the start "
            "position of your repeats."
        )
    other_options.add_argument(
        "-c", "--percent-id",
        metavar=" ",
        required=False,
        default=80,
        type=float,
        help="DEFAULT: 80. Percent ID of blast hit to consider match "
        "a candidate repeat"
        )
    other_options.add_argument(
        "-s", "--snp-thresh",
        metavar=" ",
        type=int,
        default=0,
        required = False,
        help="Specify number of SNPs to consider spacers the same. Default: 0"
        )
    other_options.add_argument(
        "--min-shared",
        metavar=" ",
        default=1,
        type=int,
        help="minimum number of spacers shared to draw an edge in the output "
        "network"
        )
    other_options.add_argument(
        "--append",
        action='store_true',
        help="Indicate that you want to add new CRISPR data to a previous run"
        )

    blast_options = parser.add_argument_group("BLASTn settings")
    blast_options.add_argument(
        "-e", "--evalue",
        metavar=" ",
        required=False,
        default='10',
        help="DEFAULT: 10. set the evalue for blastn"
        )
    blast_options.add_argument(
        "-m", "--max-target-seqs",
        metavar=" ",
        required=False,
        default='10000',
        help="DEFAULT: 10000. Set the max_target_seqs option for blastn"
        )
    blast_options.add_argument(
        "-b", "--batch-size",
        metavar=" ",
        required=False,
        default=1000,
        type=int,
        help="DEFAULT: 1000. batch size for blastdbcmd. Higher values run "
        "quicker. Reduce if you get and error about argument list length"
        )
    blast_options.add_argument(
        "-x", "--blast-options",
        metavar=" ",
        dest="other_blast_options",
        required=False,
        default='',
        help="DEFAULT: none. Input additional blastn options. "
        "Forbidden options: blastn -query -db -task -outfmt -num_threads "
        "-max_target_seqs -evalue"
        )
    return parser


def main(args):

    outdir = args.outdir + '/' if args.outdir[-1] != '/' else args.outdir

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Make sure repeats aren't similar enough to return overlapping hits
    check_repeat_similarity(args.repeats_file)
        
    if args.regex_pattern:
        assembly_names = [a for a in subprocess.run(
            "blastdbcmd -db {} -entry all -outfmt '\%a' | grep -o{} '{}' | \
            sort | uniq".format(
                args.blast_db_path,
                args.regex_type,
                args.regex_pattern
                ),
            shell=True,
            universal_newlines=True,
            capture_output=True
            ).stdout.split('\n') if len(a) > 0]

    elif args.assemblies:
        # store assembly names in args.assemblies for easy delivery
        assembly_names = file_handling.read_assembly_list_file(args.assemblies)
        args.assemblies = assembly_names

    else:
        sys.stderr.write(
            "WARNING!!!\n\n"
            "In order to associate arrays found in different contigs from the "
            "same assembly, you must provide either a list of assembly names "
            "or a regex to identify your assemblies. Each sequence in your "
            "blastdb will be treated as a separate assembly.\n")
        assembly_names = [a for a in subprocess.run(
            "blastdbcmd -db {} -entry all -outfmt '%a'".format(
                args.blast_db_path
                ),
            shell=True,
            universal_newlines=True,
            capture_output=True
            ).stdout.split('\n') if len(a) > 0]
        args.regex_pattern = ".+"
    # Find the names of all the genomes being searched by looking for 
    # the user-provided regex in the blast-db sequence IDs.
    # for all genomes, start their entry in the AssemblyCRISPRs with
    # the placeholder False. This indicated no CRISPRs found.
    # If a CRISPR array is found later this entry will be overwritten
    # with infor about the arrays.
    all_assemblies_dict = {
        k: file_handling.AssemblyCRISPRs() for k in assembly_names}
    for k,v in all_assemblies_dict.items():
        v.accession = k

    all_arrays = blastn_to_arrays(args)

    for array in all_arrays:
        assembly = all_assemblies_dict[array.genome]
        assembly.array_count+=1
        assembly.has_crispr = True
        assembly.arrays[assembly.array_count] = array

    # unpack dict into list for processing
    all_assemblies = [i for i in all_assemblies_dict.values()]

    if args.append:
        # Read in previous CRISPR_spacers.fna and reverse dict
        prev_spacer_id_dict = {
            v:k for k,v in file_handling.fasta_to_dict(
                outdir + "CRISPR_spacers.fna"
                ).items()
            }
        
        # Same for Array_IDs.txt
        prev_array_dict = {
            " ".join(v):k for k,v in file_handling.read_array_file(
                outdir + "Array_seqs.txt"
                ).items()
            }

        (non_red_spacer_dict,
            non_red_spacer_id_dict,
            non_red_array_dict,
            non_red_array_id_dict,
            cluster_reps_dict,
            rev_cluster_reps_dict
        ) = sequence_operations.non_redundant_CR(
            all_assemblies,
            args.snp_thresh,
            prev_spacer_id_dict,
            prev_array_dict,
            outdir
            )

    else:
        (non_red_spacer_dict,
            non_red_spacer_id_dict,
            non_red_array_dict,
            non_red_array_id_dict,
            cluster_reps_dict,
            rev_cluster_reps_dict
        ) = sequence_operations.non_redundant_CR(
            all_assemblies,
            args.snp_thresh,
            outdir=outdir)


    # Output summary info
    sys.stderr.write("Total unique spacers: %i\n" %len(non_red_spacer_id_dict))
    sys.stderr.write("Total unique arrays: %i\n" %len(non_red_array_id_dict))

    # Fill in spacer and array ID info
    sequence_operations.add_ids(
        all_assemblies,
        non_red_spacer_id_dict,
        non_red_array_id_dict,
        rev_cluster_reps_dict)

    # Save array files
    file_handling.write_CRISPR_files(
        all_assemblies,
        non_red_spacer_id_dict,
        non_red_array_id_dict,
        cluster_reps_dict,
        outdir,
        args.append,
        args.min_shared)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Given fasta CRISPR repeats and a blast db of genomes, \
        pulls out spacer arrays of >= 2 spacers")
    parser = build_parser(parser)

    if len(sys.argv) == 1:
        parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

        main(args) 
